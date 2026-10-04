"""R2 + R3 (Part II, E_dual): chuỗi phát xạ đầy đủ có cạnh tranh n/p/α/γ, σ_inv(ε), so tuyệt đối và kiểm định χ² có
sai số thực nghiệm.

Số liệu hạt nhân THẬT (exfor_data/): năng lượng liên kết, S_n, S_p, độ phổ cập (AME/IAEA livechart, nuclides.json, tải
thêm khi thiếu); tiết diện phi đàn hồi đo được σ_non (EXFOR, sigma_non.json): Cu 14.1 MeV 1.486 b, Ni 14.1 MeV 1.397 b.
Mô hình (mỗi sự kiện):
  hợp phần (Z,A) với E* = E_inc·A_t/(A_t+1) + S_proj; giai đoạn exciton n0 = 3 (p = 2, h = 1), mỗi bước n → n+2 do
  (a) va chạm Carnot của động cơ (chuỗi thật của mô hình chặn gãy Γ0 = 2, hiệu chỉnh theo Part I) hoặc (b) Williams K = 400;
  phát n, p (exciton, hệ số điện tích R_ν(p)), α (chỉ khi đã cân bằng, hệ số tạo cụm φ_α = 0.1), γ (GDR Brink–Axel);
  σ_n(ε) Dostrovsky, σ_p, σ_α có hàng rào Coulomb; sau lần phát đầu: bay hơi Weisskopf–Ewing (ρ ∝ e^{2√(aU)}/U^{5/4},
  a = A/8) tới khi U dưới mọi ngưỡng. (c) bay hơi thuần (bắt đầu đã cân bằng).
R3: sai số từng bin thực nghiệm truyền qua tích phân góc; sai số mô hình Poisson × hệ số thiết kế (sự kiện/chuỗi động cơ);
  χ²/bậc tự do của phổ chuẩn hóa và p-value (Wilson–Hilferty); so tuyệt đối ∫dσ/dE trong cửa sổ.
Chạy: python3 giai_so_part2_mo5.py   (≈ 10 phút, 8 lõi)
"""
import csv
import io
import json
import math
import os
import subprocess
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, "exfor_data")
T0 = 2 * math.pi
HBAR, HBARC = 6.582119569e-22, 197.327
TU = HBAR / 10.0
MASS = {"n": 939.565, "p": 938.272, "a": 3727.379}
SPIN = {"n": 2, "p": 2, "a": 1}
ZA = {"n": (0, 1), "p": (1, 1), "a": (2, 4)}
SYM = dict(enumerate("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr".split(), start=1))
_NUC = None


def nuc(Z, A):
    global _NUC
    if _NUC is None:
        _NUC = json.load(open(os.path.join(DIR, "nuclides.json")))
    key = f"{A}{SYM[Z]}"
    if key not in _NUC:
        r = subprocess.run(["curl", "-s", "-m", "60", f"https://nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides={key.lower()}"],
                           capture_output=True, text=True).stdout
        g = list(csv.DictReader(io.StringIO(r)))[0]
        f = lambda k: float(g[k]) if g.get(k) not in (None, "") else 0.0
        _NUC[key] = dict(Z=Z, A=A, BE=f("binding") * A / 1000, sn=f("sn") / 1000, sp=f("sp") / 1000, abundance=None)
        json.dump(_NUC, open(os.path.join(DIR, "nuclides.json"), "w"), indent=1)
    return _NUC[key]


def sep(Z, A, x):
    if x == "n":
        return nuc(Z, A)["sn"]
    if x == "p":
        return nuc(Z, A)["sp"]
    return nuc(Z, A)["BE"] - nuc(Z - 2, A - 4)["BE"] - 28.2957          # S_α = BE(cha) − BE(con) − BE(α)


def sigma_inv(x, eps, Zr, Ar):
    """σ_inv(ε) (fm²) của hạt x trên hạt nhân còn lại (Zr, Ar): Dostrovsky cho n; Coulomb cổ điển cho p, α."""
    R = 1.5 * Ar**(1 / 3)
    sg = math.pi * R * R
    if x == "n":
        if os.environ.get("E_DUAL_OM") == "1":                     # σ_R thế quang học KD 2003 (giai_so_part2_mo6 H2)
            tab = _om_table()
            ia = int(np.argmin(np.abs(tab["A"] - Ar)))
            return np.interp(eps, tab["eps"], tab["sigma"][ia])
        al = 0.76 + 2.2 * Ar**(-1 / 3); be = (2.12 * Ar**(-2 / 3) - 0.050) / al
        return sg * al * (1 + be / np.maximum(eps, 1e-3))
    z = ZA[x][0]
    V = 1.44 * z * Zr / (R + (1.2 if x == "a" else 0.0))
    k, c = (0.70, 0.30) if x == "p" else (0.85, 0.0)
    return np.where(eps > k * V, sg * (1 + c) * (1 - k * V / np.maximum(eps, 1e-6)), 0.0)


_OM = None


def _om_table():
    global _OM
    if _OM is None:
        z = np.load(os.path.join(DIR, "om_sigma_n.npz"))
        _OM = {k: z[k] for k in ("A", "eps", "sigma")}
    return _OM


def gdr(Eg, Z, A):
    EG = 31.2 * A**(-1 / 3) + 20.6 * A**(-1 / 6); G = 5.0
    N = A - Z
    s0 = 2 * (60 * N * Z / A) / (math.pi * G) * 0.1                       # mb → fm²
    return s0 * Eg**2 * G**2 / ((Eg**2 - EG**2)**2 + Eg**2 * G**2)


def kfac(x):
    return SPIN[x] * MASS[x] / (math.pi**2 * HBARC**2 * HBAR) * TU          # nhân với ∫ ε σ (…) dε (MeV², fm²) → 1/đv t


# ----------------------------------------------------------- giai đoạn exciton
def exciton_rates(n, p, E, Z, A, proj, cache):
    key = (n, p, round(E, 3), Z, A)
    if key in cache:
        return cache[key]
    g = A / 13.0
    N = A - Z
    out = {}
    for x in ("n", "p"):
        B = sep(Z, A, x)
        if E <= B:
            out[x] = (0.0, None); continue
        Zr, Ar = Z - ZA[x][0], A - ZA[x][1]
        frac = N / A if x == "n" else Z / A
        own = 1.0 if x == proj else 0.0
        R = (own + (p - 1) * frac) / p
        eps = np.linspace(1e-3, E - B, 300)
        U = E - B - eps
        w = eps * sigma_inv(x, eps, Zr, Ar) * p * R * (n - 1) / g * U**(n - 2) / E**(n - 1)
        out[x] = (kfac(x) * np.trapezoid(w, eps), (eps, w))
    Eg = np.linspace(1e-3, E, 300)
    wg = Eg**2 * gdr(Eg, Z, A) * ((E - Eg) / E)**(n - 1)
    out["g"] = (np.trapezoid(wg, Eg) / (math.pi**2 * HBARC**2 * HBAR) * TU, (Eg, wg))
    cache[key] = out
    return out


def draw(rng, grid):
    e, w = grid
    c = np.cumsum(np.clip(w, 0, None)); c /= c[-1]
    return float(np.interp(rng.random(), c, e))


# ----------------------------------------------------------- bay hơi cân bằng
PAIRING = True
SYS = float(os.environ.get("E_DUAL_SYS", "0.08"))       # sai số hệ thống/bin, lấy từ độ lệch Baba–Takahashi (mo6 H1)


def rho(U, A, Z=None):
    """Fermi gas dịch chuyển cặp: U_eff = U − δ, δ = 12/√A × (2 chẵn–chẵn, 1 A lẻ, 0 lẻ–lẻ)."""
    a = A / 8.0
    if PAIRING and Z is not None:
        N = A - Z
        d = 12 / math.sqrt(A) * ((Z % 2 == 0) + (N % 2 == 0))
        U = U - d
    return np.where(U > 0.05, np.exp(2 * np.sqrt(a * np.maximum(U, 0.05))) / np.maximum(U, 0.05)**1.25, 0.0)


def evaporate(Z, A, U, rng, neut, counts, phi_a=0.1):
    while True:
        W, grids = {}, {}
        for x in ("n", "p", "a"):
            B = sep(Z, A, x)
            if U <= B:
                continue
            Zr, Ar = Z - ZA[x][0], A - ZA[x][1]
            eps = np.linspace(1e-3, U - B, 200)
            w = SPIN[x] * MASS[x] * eps * sigma_inv(x, eps, Zr, Ar) * rho(U - B - eps, Ar, Zr) * (phi_a if x == "a" else 1.0)
            W[x] = np.trapezoid(w, eps); grids[x] = (eps, w)
        Eg = np.linspace(1e-3, U, 200)
        wg = Eg**2 * gdr(Eg, Z, A) * rho(U - Eg, A, Z)
        W["g"] = np.trapezoid(wg, Eg)
        tot = sum(W.values())
        if tot <= 0 or len(W) == 1:
            return                                             # chỉ còn γ: dừng
        r = rng.random() * tot
        for x, wv in W.items():
            r -= wv
            if r <= 0:
                break
        if x == "g":
            return
        e = draw(rng, grids[x])
        counts[x] += 1
        if x == "n":
            neut.append(e)
        U = U - sep(Z, A, x) - e
        Z, A = Z - ZA[x][0], A - ZA[x][1]


# ----------------------------------------------------------- một sự kiện
def event(rng, targets, proj, Einc, mode, seqs, cache, counts, K=400.0):
    ab = np.array([t[1] for t in targets]); i = rng.choice(len(targets), p=ab / ab.sum())
    Zt, At = targets[i][0]
    Zc, Ac = Zt + ZA[proj][0], At + 1
    E = Einc * At / (At + 1) + sep(Zc, Ac, proj)
    neut = []
    if mode == "evap":
        evaporate(Zc, Ac, E, rng, neut, counts); return neut
    g = Ac / 13.0
    neq = int(round(math.sqrt(1.39 * g * E)))
    n, p, t = 3, 2, 0.0
    seq = list(seqs[rng.integers(len(seqs))]) if mode == "engine" else None
    k = 0
    while True:
        if n >= neq:
            evaporate(Zc, Ac, E, rng, neut, counts); return neut
        rates = exciton_rates(n, p, E, Zc, Ac, proj, cache)
        lam = sum(r[0] for r in rates.values())
        if mode == "engine":
            tnext = seq[k] if k < len(seq) else np.inf
        else:
            lp = 2 * math.pi / HBAR * (K / (Ac**3 * E)) * g**3 * E**2 / (2 * (n + 1)) * TU
            tnext = t + rng.exponential(1 / lp)
        te = t + rng.exponential(1 / lam) if lam > 0 else np.inf
        if te < tnext:
            r = rng.random() * lam
            for x, (lv, grid) in rates.items():
                r -= lv
                if r <= 0:
                    break
            e = draw(rng, grid)
            if x == "g":
                E -= e
                if E <= 0:
                    return neut
                t = te; continue
            counts[x] += 1
            if x == "n":
                neut.append(e)
            U = E - sep(Zc, Ac, x) - e
            evaporate(Zc - ZA[x][0], Ac - ZA[x][1], U, rng, neut, counts)
            return neut
        t = tnext; k += 1; n += 2; p += 1


# ----------------------------------------------------------- dữ liệu có sai số
def load(fname):
    d = json.load(open(os.path.join(DIR, fname)))["data_tables"]["data"]
    return {h.strip(): np.array([np.nan if v is None else v for v in c], float) for h, c in zip(d["heads"], d["data"])}


def ang_int_err(ang, emid, val, err):
    es, s, se = [], [], []
    for e in np.unique(emid):
        m = (emid == e) & np.isfinite(val)
        if m.sum() < 3:
            continue
        mu = np.cos(np.radians(ang[m])); o = np.argsort(mu)
        mu, v, er = mu[o], val[m][o], np.nan_to_num(err[m][o], nan=0.0)
        mu2 = np.concatenate([[-1.0], mu, [1.0]]); v2 = np.concatenate([[v[0]], v, [v[-1]]])
        er2 = np.concatenate([[er[0]], er, [er[-1]]])
        wts = np.zeros(len(mu2)); d = np.diff(mu2)
        wts[:-1] += d / 2; wts[1:] += d / 2
        es.append(e); s.append(2 * math.pi * np.sum(wts * v2)); se.append(2 * math.pi * math.sqrt(np.sum((wts * er2)**2)))
    return np.array(es), np.array(s), np.array(se)


def col(D, prefix):
    return next(v for k, v in D.items() if k.replace(" ", "").startswith(prefix.replace(" ", "")))


def datasets():
    out = []
    Cu = [((29, 63), 69.15), ((29, 65), 30.85)]
    Ni = [((28, 58), 68.08), ((28, 60), 26.22), ((28, 61), 1.14), ((28, 62), 3.63), ((28, 64), 0.93)]
    Ge = [((32, 70), 20.52), ((32, 72), 27.45), ((32, 73), 7.76), ((32, 74), 36.52), ((32, 76), 7.75)]
    sn = json.load(open(os.path.join(DIR, "sigma_non.json")))
    D = load("22077_002.json")
    for En in (14.1, 18.0):
        m = np.isclose(D["EN"], En)
        es, s, se = ang_int_err(D["ANG"][m], 0.5 * (D["E-MIN"][m] + D["E-MAX"][m]), 1000 * D["DATA"][m], 1000 * D["DATA-ERR"][m])
        snon = sn.get(f"Cu_{En}", [1.43, 0.05])[0] * 1000 if En == 14.1 else 1430.0
        out.append((f"Cu(n,xn) {En} — Baba 1988", Cu, "n", En, En - 4, es, s, se, snon))
    D = load("22075_068.json")
    es, s, se = ang_int_err(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * col(D, "DATA1"), 1000 * col(D, "DATA-ERR1"))
    out.append(("Cu(n,xn) 14.1 — Takahashi 1988", Cu, "n", 14.1, 10.1, es, s, se, sn["Cu_14.1"][0] * 1000))
    D = load("22278_002.json")
    es, s, se = ang_int_err(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * col(D, "DATA1"), 1000 * col(D, "ERR-S1"))
    out.append(("Ge(n,xn) 14.1 — Gotoh 1992", Ge, "n", 14.1, 10.1, es, s, se, None))
    for f, En in (("22025_044.json", 14.1), ("22025_053.json", 18.0)):
        D = load(f)
        es, s, se = ang_int_err(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * D["DATA"], 1000 * D["DATA-ERR"])
        out.append((f"Ni(n,xn) {En} — Yabuta 1986", Ni, "n", En, En - 4, es, s, se,
                    sn["Ni_14.1"][0] * 1000 if En == 14.1 else None))
    D = load("A0068_024.json")
    es, s, se = ang_int_err(D["ANG"], D["E"], D["DATA"], D["ERR-T"])
    out.append(("60Ni(p,xn) 22 — Biryukov 1980", [((28, 60), 1.0)], "p", 22.0, 15.0, es, s, se, None))
    return out


def chi2_p(chi2, k):
    z = ((chi2 / k)**(1 / 3) - (1 - 2 / (9 * k))) / math.sqrt(2 / (9 * k))
    return 0.5 * math.erfc(z / math.sqrt(2))


def compare(args):
    label, targets, proj, Einc, hi, es, s, se, snon, mode, seqs, nev, seed = args
    rng = np.random.default_rng(seed)
    cache, counts, allE = {}, {"n": 0, "p": 0, "a": 0}, []
    for _ in range(nev):
        allE += event(rng, targets, proj, Einc, mode, seqs, cache, counts)
    allE = np.array(allE)
    m = (es >= 1.0) & (es <= hi) & np.isfinite(s) & (se > 0)
    e, d, de = es[m], s[m], se[m]
    edges = np.concatenate([[e[0] - (e[1] - e[0]) / 2], 0.5 * (e[1:] + e[:-1]), [e[-1] + (e[-1] - e[-2]) / 2]])
    width = np.diff(edges)
    cnt, _ = np.histogram(allE, edges)
    design = nev / len(seqs) if mode == "engine" else 1.0
    mod = cnt / (nev * width); mod_e = np.sqrt(np.maximum(cnt, 1) * design) / (nev * width)   # neutron / sự kiện / MeV
    # hình dạng: chuẩn hóa diện tích
    A_d = np.sum(d * width); A_m = np.sum(mod * width)
    dn, den = d / A_d, de / A_d
    mn, men = mod / A_m, mod_e / A_m
    chi2 = float(np.sum((dn - mn)**2 / (den**2 + men**2 + (0.02 * dn)**2)))
    chi2d = float(np.sum((dn - mn)**2 / (den**2 + (0.02 * dn)**2)))          # chỉ sai số thực nghiệm (công bằng)
    chi2s = float(np.sum((dn - mn)**2 / (den**2 + (SYS * dn)**2)))           # + hệ thống giữa các thí nghiệm (H1)
    k = len(e) - 1
    res = dict(label=label, mode=mode, chi2k=chi2 / k, p=chi2_p(chi2, k), chi2dk=chi2d / k, pd=chi2_p(chi2d, k), chi2sk=chi2s / k, ps=chi2_p(chi2s, k), mean_d=np.sum(e * dn * width),
               mean_m=np.sum(e * mn * width), mult_n=counts["n"] / nev, frac_p=counts["p"] / nev, frac_a=counts["a"] / nev)
    res["snon"] = snon
    if snon:
        res["abs_d"] = A_d; res["abs_m"] = snon * A_m
        res["abs_err"] = math.sqrt(np.sum((de * width)**2) + (0.05 * A_d)**2)
    return res


def engine_sequences():
    f = os.path.join(DIR, "engine_seq_G2.npz")
    if os.path.exists(f):
        z = np.load(f, allow_pickle=True); return list(z["seqs"])
    import giai_so_part2_mo3 as M
    with Pool() as pool:
        res = pool.map(M._traj, [(0.3, 0.5, 2.0, 9700 + i, 2000.0, 0.01) for i in range(24)])
    seqs = M.sequences(res, 2000 / 6)
    np.savez(f, seqs=np.array(seqs, dtype=object))
    return seqs


def main(nev=6000):
    seqs = engine_sequences()
    rng = np.random.default_rng(0)                        # nạp trước số liệu hạt nhân (tránh ghi đồng thời)
    for ds in datasets():
        for mode in ("engine", "std", "evap"):
            for _ in range(150):
                event(rng, ds[1], ds[2], ds[3], mode, seqs, {}, {"n": 0, "p": 0, "a": 0})
    jobs = []
    for ds in datasets():
        for mode in ("engine", "std", "evap"):
            jobs.append(ds + (mode, seqs, nev, 101))
    with Pool() as pool:
        res = pool.map(compare, jobs)
    print(f"=== R2 + R3. Chuỗi phát xạ đầy đủ (n/p/α/γ, σ_inv(ε), số liệu AME) vs EXFOR, χ² có sai số ({nev} sự kiện/mô hình) ===")
    print(f"động cơ: {len(seqs)} chuỗi va chạm (mô hình chặn gãy Γ0 = 2, Γc = 0.5, p* = 0.3)")
    names = {"engine": "động cơ (Part I hiệu chỉnh)", "std": "exciton chuẩn K = 400", "evap": "bay hơi thuần"}
    last = None
    for r in res:
        if r["label"] != last:
            print(f"--- {r['label']}: ⟨ε⟩ dữ liệu {r['mean_d']:.2f} MeV" +
                  (f"; ∫dσ/dE dữ liệu {r['abs_d']:.0f} ± {r['abs_err']:.0f} mb" if 'abs_d' in r else ""))
            print(f"{'':>30} {'⟨ε⟩':>6} {'χ²/k':>7} {'p':>9} {'χ²/k (chỉ TN)':>14} {'p':>9} {'χ²/k (+8% HT)':>14} {'p':>9} {'bội n':>6} {'p/sk':>6} {'α/sk':>6} {'∫dσ/dE (mb)':>12} {'σ(n,xp) (mb)':>13}")
            last = r["label"]
        ab = f"{r['abs_m']:.0f}" if "abs_m" in r else "—"
        sp_ = f"{r['frac_p']*r['snon']:.0f}" if r.get("snon") else "—"
        print(f"{names[r['mode']]:>30} {r['mean_m']:6.2f} {r['chi2k']:7.2f} {r['p']:9.2e} {r['chi2dk']:14.2f} {r['pd']:9.2e} {r['chi2sk']:14.2f} {r['ps']:9.2e} "
              f"{r['mult_n']:6.2f} {r['frac_p']:6.3f} {r['frac_a']:6.3f} {ab:>12} {sp_:>13}")
    print("(χ² dùng sai số thực nghiệm từng bin (truyền qua tích phân góc) + 2% không tương quan + sai số MC của mô hình;"
          " ∫ trong cùng cửa sổ năng lượng; σ_non từ EXFOR)\n")


if __name__ == "__main__":
    main()
