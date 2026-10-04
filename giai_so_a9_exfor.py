"""A9: đối chiếu phổ nucleon phát ra của học thuyết E_dual với dữ liệu thực nghiệm EXFOR (IAEA).

Dữ liệu (tải bằng API https://nds.iaea.org/dataexplorer/api/exfor/subentry/…, lưu ở exfor_data/):
  Cu(n,xn)  14.1, 18.0 MeV   M.Baba 1988          EXFOR 22077-002   (bộ chính)
  Cu(n,xn)  14.1 MeV         A.Takahashi 1988     EXFOR 22075-068   (kiểm chéo, con trỏ 1)
  Ge(n,xn)  14.1 MeV         M.Gotoh 1992         EXFOR 22278-002   (A ≈ 72, con trỏ 1)
  Ni(n,xn)  14.1, 18 MeV     N.Yabuta 1986        EXFOR 22025-044, -053
  60Ni(p,xn) 22 MeV          N.S.Biryukov 1980    EXFOR A0068-024
Xử lý: tích phân theo góc dσ/dE = 2π∫ d²σ/dΩdE dμ (hình thang theo μ = cosθ, ngoại suy phẳng tới 0° và 180°);
vùng liên tục ε ∈ [1, E_max] với E_max = E_inc − 4 MeV cho (n,xn) (bỏ đàn hồi + mức rời rạc), 15 MeV cho (p,xn).
Mô hình (cùng cửa sổ, cùng công thức exciton Eq. (22) Part I: ε = (E_exc − S_n)·Beta(2, n−1), n = 2k+1, tại bước k
đầu tiên vượt S_n = 8 MeV):
  cổ điển Part I (đợt va chạm Carnot từ mô phỏng Langevin Bảng I, 1 MeV),
  lượng tử V.C.1 (Davies mọi lúc, MCWF), lượng tử chặn gãy (Q1, p* = 0.3),
  bay hơi Weisskopf ε e^{−ε/T}, T = √(E*/a), a = A/8, E* = E_inc·A/(A+1) + S_n.
Thước đo: ⟨ε⟩, tỉ phần ε > 8 MeV, khoảng cách Kolmogorov–Smirnov D giữa phân bố chuẩn hóa của mô hình và dữ liệu.
Chạy: python3 giai_so_a9_exfor.py   (≈ 5 phút, 8 lõi)
"""
import json
import math
import os
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
HW = 10.0
SN = 8.0
DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "exfor_data")


def load(fname):
    d = json.load(open(os.path.join(DIR, fname)))["data_tables"]["data"]
    return {h.strip(): np.array([np.nan if v is None else v for v in c], float) for h, c in zip(d["heads"], d["data"])}, \
        [u.strip() for u in d["units"]]


def ang_int(ang, emid, val):
    """dσ/dE (mb/MeV) theo từng bin năng lượng phát ra."""
    out = {}
    for e in np.unique(emid):
        m = (emid == e) & np.isfinite(val)
        if m.sum() < 3:
            continue
        a = np.radians(ang[m]); v = val[m]
        o = np.argsort(np.cos(a)); mu = np.cos(a)[o]; v = v[o]
        mu = np.concatenate([[-1.0], mu, [1.0]]); v = np.concatenate([[v[0]], v, [v[-1]]])
        out[e] = 2 * math.pi * np.trapezoid(v, mu)
    es = np.array(sorted(out)); return es, np.array([out[e] for e in es])


def datasets():
    sets = []
    D, u = load("22077_002.json")                                        # Baba: EN, ANG, E-MIN, E-MAX, DATA (b/sr/MeV)
    for En in (14.1, 18.0):
        m = np.isclose(D["EN"], En)
        es, s = ang_int(D["ANG"][m], 0.5 * (D["E-MIN"][m] + D["E-MAX"][m]), 1000 * D["DATA"][m])
        sets.append((f"Cu(n,xn) {En} MeV — Baba 1988", 64.5, En, En - 4, es, s, 7.6))
    D, u = load("22075_068.json")
    es, s = ang_int(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * D["DATA      1"])
    sets.append(("Cu(n,xn) 14.1 MeV — Takahashi 1988", 64.5, 14.1, 10.1, es, s, 7.6))
    D, u = load("22278_002.json")
    es, s = ang_int(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * D["DATA      1"])
    sets.append(("Ge(n,xn) 14.1 MeV — Gotoh 1992", 73.6, 14.1, 10.1, es, s, 7.4))
    for f, En in (("22025_044.json", 14.1), ("22025_053.json", 18.0)):
        D, u = load(f)
        es, s = ang_int(D["ANG"], 0.5 * (D["E-MIN"] + D["E-MAX"]), 1000 * D["DATA"])
        sets.append((f"Ni(n,xn) {En} MeV — Yabuta 1986", 59.7, En, En - 4, es, s, 9.0))
    D, u = load("A0068_024.json")                                        # Biryukov: ANG, E, DATA (mb/sr/MeV)
    es, s = ang_int(D["ANG"], D["E"], D["DATA"])
    sets.append(("60Ni(p,xn) 22 MeV — Biryukov 1980", 61.0, 22.0, 15.0, es, s, 4.8))   # S_p(61Cu) = 4.8
    return sets


def stats_from_spec(es, s, lo, hi):
    m = (es >= lo) & (es <= hi) & np.isfinite(s)
    e, w = es[m], np.clip(s[m], 0, None)
    w = w / np.trapezoid(w, e)
    mean = np.trapezoid(e * w, e)
    tail = np.trapezoid(w[e > 8], e[e > 8]) if np.sum(e > 8) > 1 else 0.0
    cdf = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(e))])
    return mean, tail, e, cdf / cdf[-1]


def stats_from_samples(x, lo, hi, egrid):
    x = x[(x >= lo) & (x <= hi)]
    cdf = np.searchsorted(np.sort(x), egrid, side="right") / max(len(x), 1)
    return x.mean(), np.mean(x > 8), cdf, len(x)


# --------------------------------------------------------------- mẫu mô hình
def exciton_samples(bursts, rng, per=40):
    """bursts: danh sách năng lượng va chạm (MeV) trong từng đợt → mẫu năng lượng nucleon phát ra."""
    out, nb_emit = [], 0
    for b in bursts:
        acc = np.cumsum(b)
        hit = np.nonzero(acc > SN)[0]
        if len(hit) == 0:
            continue
        nb_emit += 1
        k = hit[0] + 1; n = 2 * k + 1
        out.append((acc[hit[0]] - SN) * rng.beta(2, n - 1, size=per))
    return np.concatenate(out) if out else np.array([]), nb_emit / max(len(bursts), 1)


def group(times, energies, gap=T0):
    bursts, cur, last = [], [], None
    for t, e in zip(times, energies):
        if last is not None and t - last >= gap:
            bursts.append(cur); cur = []
        cur.append(e); last = t
    if cur:
        bursts.append(cur)
    return bursts


def _classical(seed):
    from giai_so_dong_ho_photon import _traj
    from giai_so_burst_g import cfg_for
    _, ev = _traj((cfg_for(0.1), 0.1, seed, 200 * T0, 0.05))
    ks = [(e[1], e[3] * HW) for e in ev if e[0] == "khóa" and e[1] > 5 * T0]
    return group([k[0] for k in ks], [k[1] for k in ks])


def _vc1(seed):
    from giai_so_part2_master import _traj
    jumps = _traj((seed, 1500.0, 0.01, False))[0]
    cj = [(j[0], j[2] * HW) for j in jumps if j[1] == "Carnot"]
    return group([c[0] for c in cj], [c[1] for c in cj])


def _gated(seed):
    from giai_so_part2_mo2 import _q1_traj
    jumps = [j for j in _q1_traj((0.3, seed, 2400.0, 0.01)) if j[0] > 400]
    return group([j[0] for j in jumps], [j[1] * HW for j in jumps])


def main():
    rng = np.random.default_rng(7)
    with Pool() as pool:
        cl = [b for r in pool.map(_classical, range(8100, 8116)) for b in r]
        q1 = [b for r in pool.map(_vc1, range(8200, 8224)) for b in r]
        q2 = [b for r in pool.map(_gated, range(8300, 8324)) for b in r]
    models = {}
    for name, bursts in (("cổ điển Part I", cl), ("lượng tử V.C.1", q1), ("lượng tử chặn gãy (Q1)", q2)):
        x, pe = exciton_samples(bursts, rng)
        models[name] = (x, pe, len(bursts), np.mean([sum(b) for b in bursts]))
    print("=== A9. Phổ neutron phát ra: mô hình E_dual vs dữ liệu EXFOR ===")
    for name, (x, pe, nb, eb) in models.items():
        print(f"{name}: {nb} đợt, năng lượng đợt TB {eb:.1f} MeV, xác suất đợt phát nucleon {pe:.1%}, {len(x)} mẫu ε")
    print()
    for label, A, Einc, hi, es, s, S in datasets():
        lo = 1.0
        mean_d, tail_d, eg, cdf_d = stats_from_spec(es, s, lo, hi)
        Estar = Einc * (A - 1) / A + S
        T = math.sqrt(Estar / (A / 8))
        eg2 = np.linspace(lo, hi, 400)
        w = eg2 * np.exp(-eg2 / T); w /= np.trapezoid(w, eg2)
        cw = np.concatenate([[0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(eg2))]); cw /= cw[-1]
        mean_w = np.trapezoid(eg2 * w, eg2); tail_w = np.trapezoid(w[eg2 > 8], eg2[eg2 > 8]) if hi > 8 else 0
        Dw = np.max(np.abs(np.interp(eg, eg2, cw) - cdf_d))
        print(f"--- {label}: E* ≈ {Estar:.1f} MeV, cửa sổ {lo}–{hi:.1f} MeV, {len(eg)} bin")
        print(f"{'':>26} {'⟨ε⟩ (MeV)':>10} {'P(ε>8 MeV)':>11} {'KS D':>7}")
        print(f"{'dữ liệu':>26} {mean_d:10.2f} {tail_d:11.3f} {'—':>7}")
        print(f"{f'Weisskopf (T = {T:.2f})':>26} {mean_w:10.2f} {tail_w:11.3f} {Dw:7.3f}")
        for name, (x, pe, nb, eb) in models.items():
            if len(x) == 0:
                print(f"{name:>26} {'không phát':>10}"); continue
            mm, tt, cm, nn = stats_from_samples(x, lo, hi, eg)
            Dm = np.max(np.abs(cm - cdf_d))
            print(f"{name:>26} {mm:10.2f} {tt:11.3f} {Dm:7.3f}")
        print()


if __name__ == "__main__":
    main()
