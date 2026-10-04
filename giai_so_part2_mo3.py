"""Ba câu hỏi mở sau A9 (Part II, E_dual, HocThuyetE (7).pdf).

C3  Thống kê đợt của mô hình chặn gãy (Q1): quét tốc độ bơm Γ0 và làm nguội Γc (p* = 0.3) → va chạm/đợt, thời gian lưu
    (TB, CV), ΔE/va chạm, ngưỡng ⟨|ℓ|⟩; chọn bộ gần Part I nhất (5.12; 4.5–6 T0, CV ≈ 1; 0.40–0.46 ħω0; p_θ = 3.56).
C1+C2  Exciton do động cơ điều khiển: E_exc = E* của phản ứng (gắn E*), trạng thái đầu n = 3 (p = 2, h = 1); mỗi va chạm
    Carnot của động cơ (thời điểm lấy từ chuỗi va chạm thật của mô hình chặn gãy) ⟹ n → n + 2; giữa các va chạm phát
    nucleon với tốc độ Griffin–Cline tuyệt đối
        λ_e(n,E) = (2s+1) μc² σ_inv / (π² (ħc)² ħ) · ∫ ε · p(n−1)/g · U^{n−2}/E^{n−1} dε,  U = E − S_n − ε,
    σ_inv = π(1.2 A^{1/3})² fm², g = A/13; phổ từ trạng thái n là (E − S_n)·Beta(2, n−1) — đúng Eq. (22) Part I;
    phát từ n nhỏ cho đuôi tiền cân bằng. n bị chặn ở n_eq = √(1.39 g E).
    Đối chứng: exciton chuẩn với tốc độ chuyển nội Williams λ₊ = (2π/ħ)|M|² g³E²/(2(n+1)), |M|² = K/(A³E), K = 135, 400 MeV³;
    và cân bằng tức thì (n = n_eq). So với dữ liệu EXFOR (giai_so_a9_exfor.datasets).
C4  Bơm mạnh (Γ0 = 2, 3) + phát xạ neutron thứ hai (n,2n) (bay hơi từ năng lượng còn lại): so với dữ liệu.
Chạy: python3 giai_so_part2_mo3.py [C3|C12|C4]
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

import giai_so_part2_mo as P

T0 = 2 * math.pi
HBAR = 6.582119569e-22          # MeV s
TU = HBAR / 10.0                # 1 đơn vị thời gian của mã = ħ/ħω0 (s)
HBARC, MN = 197.327, 939.565


# ======================================================================= C3
def gated(pstar, Gc, G0):
    m, H, Pout = P.a1_model()
    E, V = np.linalg.eigh(H)
    po = np.array([float(np.real(V[:, j].conj() @ Pout @ V[:, j])) for j in range(len(E))])
    ue = V.conj().T @ m["u"] @ V
    A = math.sqrt(Gc) * np.triu(ue, 1) * (po > pstar)[None, :]
    return m, H, [V @ A @ V.conj().T, P.pump(m, G0)], V


def _traj(args):
    pstar, Gc, G0, seed, Tmax, dt = args
    m, H, Ls, V = gated(pstar, Gc, G0)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = V[:, 0].astype(complex)
    t, r, jumps = 0.0, rng.random(), []
    while t < Tmax:
        psi = U @ psi; t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            if k == 0:
                Eb = np.real(np.vdot(psi, H @ psi)); lb = np.real(np.vdot(psi, Labs @ psi))
            psi = Ls[k] @ psi; psi /= np.linalg.norm(psi)
            if k == 0:
                jumps.append((t, Eb - np.real(np.vdot(psi, H @ psi)), lb))
            r = rng.random()
    return jumps


def burst_stats(res, Tskip):
    sizes, dws, eJ, lon = [], [], [], []
    for jl in res:
        jl = [j for j in jl if j[0] > Tskip]
        cur, last = 0, None
        for tj, e, lb in jl:
            eJ.append(e)
            if last is None or tj - last >= T0:
                if cur:
                    sizes.append(cur)
                if last is not None:
                    dws.append((tj - last) / T0)
                cur = 0; lon.append(lb)
            cur += 1; last = tj
    s, d_, e_ = np.array(sizes), np.array(dws), np.array(eJ)
    return dict(n=len(s), size=s.mean() if len(s) else np.nan, dwell=d_.mean() if len(d_) else np.nan,
                cv=d_.std() / d_.mean() if len(d_) > 1 else np.nan, eJ=e_.mean() if len(e_) else np.nan,
                lon=np.mean(lon) if lon else np.nan)


def part_c3(N=16, Tmax=1600.0, dt=0.01):
    print("=== C3. Quét (Γ0, Γc) của mô hình chặn gãy p* = 0.3 → thống kê đợt so với Part I ===")
    print(f"{'Γ0':>5} {'Γc':>5} {'đợt':>5} {'va chạm/đợt':>12} {'lưu TB (T0)':>12} {'CV':>5} {'ΔE/va chạm':>11} "
          f"{'⟨|ℓ|⟩ khởi phát':>16} {'điểm lệch':>9}")
    best = None
    with Pool() as pool:
        for G0 in (0.2, 0.5, 1.0):
            for Gc in (0.5, 2.0, 5.0):
                st = burst_stats(pool.map(_traj, [(0.3, Gc, G0, 7000 + i, Tmax, dt) for i in range(N)]), Tmax / 6)
                score = (abs(st["size"] - 5.12) / 5.12 + abs(st["dwell"] - 5.25) / 5.25 + abs(st["cv"] - 1.05)
                         + abs(st["eJ"] - 0.43) / 0.43 + abs(st["lon"] - 3.56) / 3.56) if st["n"] > 5 else np.inf
                print(f"{G0:5.1f} {Gc:5.1f} {st['n']:5d} {st['size']:12.2f} {st['dwell']:12.2f} {st['cv']:5.2f} "
                      f"{st['eJ']:11.3f} {st['lon']:16.2f} {score:9.2f}", flush=True)
                if best is None or score < best[0]:
                    best = (score, G0, Gc, st)
    print(f"bộ gần Part I nhất: Γ0 = {best[1]}, Γc = {best[2]} (điểm lệch {best[0]:.2f})")
    print("(Part I: 5.12 va chạm/đợt; lưu 4.5–6 T0, CV 1.0–1.1; 0.40–0.46 ħω0; p_θ ngưỡng = M(0)θ̇_br = 3.56)\n")
    return best


# ======================================================================= C1 + C2
def lam_emit(n, p, E, Sn, A):
    """Tốc độ phát nucleon (1/đơn vị thời gian của mã) từ trạng thái exciton n ở năng lượng kích thích E."""
    if E <= Sn:
        return 0.0
    g = A / 13.0
    sig = math.pi * (1.2 * A**(1 / 3))**2
    K0 = 2 * MN * sig / (math.pi**2 * HBARC**2 * HBAR)                 # 1/(MeV² s)
    integ = p * (E - Sn)**n / (g * n * E**(n - 1))                       # ∫ ε p(n−1)/g U^{n−2}/E^{n−1} dε (MeV²)
    return K0 * integ * TU


def lam_plus(n, E, A, K):
    g = A / 13.0
    M2 = K / (A**3 * E)
    return 2 * math.pi / HBAR * M2 * g**3 * E**2 / (2 * (n + 1)) * TU


def cascade(seq, E, Sn, A, rng, nmax):
    """seq: thời điểm va chạm tương đối (đơn vị mã) từ lúc khởi phát; trả năng lượng nucleon phát ra (hoặc None)."""
    n, p, t = 3, 2, 0.0
    for tnext in list(seq) + [np.inf]:
        lam = lam_emit(n, p, E, Sn, A)
        if lam <= 0:
            return None
        te = t + rng.exponential(1 / lam)
        if te < tnext:
            return (E - Sn) * rng.beta(2, n - 1), n
        t = tnext
        if n < nmax:
            n += 2; p += 1
    return None


def cascade_std(E, Sn, A, K, rng, nmax):
    n, p = 3, 2
    while True:
        le = lam_emit(n, p, E, Sn, A)
        lp = lam_plus(n, E, A, K) if n < nmax else 0.0
        if rng.random() < le / (le + lp):
            return (E - Sn) * rng.beta(2, n - 1), n
        n += 2; p += 1


def second_chance(eps1, E, Sn, Sn2, A, rng):
    """Neutron thứ hai (bay hơi Weisskopf ε e^{−ε/T}, T = √(U/a)) nếu U = E − Sn − ε1 > Sn2."""
    U = E - Sn - eps1
    if U <= Sn2:
        return None
    T = math.sqrt(U / (A / 8))
    x = np.linspace(0, U - Sn2, 200)                     # lấy mẫu ngược ε e^{−ε/T} cắt ở U − Sn2
    c = np.cumsum(x * np.exp(-x / T)); c /= c[-1]
    return float(np.interp(rng.random(), c, x))


def with_second(first, E, Sn, Sn2, A, rng):
    out = []
    for o in first:
        if not o:
            continue
        out.append(o[0])
        e2 = second_chance(o[0], E, Sn, Sn2, A, rng)
        if e2 is not None:
            out.append(e2)
    return np.array(out)


def sequences(res, Tskip, horizon=300 * T0):
    seqs = []
    for jl in res:
        ts = np.array([j[0] for j in jl if j[0] > Tskip])
        if len(ts) == 0:
            continue
        onset = [ts[0]] + [ts[i] for i in range(1, len(ts)) if ts[i] - ts[i - 1] >= T0]
        for t0 in onset:
            rel = ts[(ts >= t0) & (ts < t0 + horizon)] - t0
            seqs.append(rel[1:])                                         # va chạm đầu = lúc nucleon vào (n = 3)
    return seqs


def part_c12(G0=None, Gc=None, N=24, Tmax=2400.0, dt=0.01):
    import giai_so_a9_exfor as X
    rng = np.random.default_rng(13)
    cfgs = [("động cơ chặn gãy (Q1: Γ0 0.2, Γc 0.5)", 0.2, 0.5)]
    if G0 is not None and (G0, Gc) != (0.2, 0.5):
        cfgs.append((f"động cơ hiệu chỉnh C3 (Γ0 {G0}, Γc {Gc})", G0, Gc))
    engines = {}
    with Pool() as pool:
        for name, g0, gc in cfgs:
            res = pool.map(_traj, [(0.3, gc, g0, 9500 + i, Tmax, dt) for i in range(N)])
            engines[name] = sequences(res, Tmax / 6)
    print("=== C1 + C2. Exciton do động cơ điều khiển (E_exc = E*, n0 = 3, va chạm Carnot ⟹ n + 2) vs dữ liệu EXFOR ===")
    for name, sq in engines.items():
        print(f"{name}: {len(sq)} chuỗi khởi phát; va chạm trong 1 T0 đầu TB "
              f"{np.mean([np.sum(s < T0) for s in sq]):.2f}, trong 10 T0 TB {np.mean([np.sum(s < 10*T0) for s in sq]):.2f}")
    Sn_emit = {"Cu": 7.6, "Ge": 7.4, "Ni(n": 8.6, "60Ni": 11.7}
    for label, A, Einc, hi, es, s, S in X.datasets():
        Estar = Einc * (A - 1) / A + S
        Sn = next(v for k, v in Sn_emit.items() if label.startswith(k))
        nmax = int(round(math.sqrt(1.39 * (A / 13) * Estar)))
        md, td, eg, cd = X.stats_from_spec(es, s, 1.0, hi)
        print(f"--- {label}: E* = {Estar:.1f} MeV, S_n(phát) = {Sn}, n_eq = {nmax}; λ_e(n=3) = "
              f"{lam_emit(3, 2, Estar, Sn, A)*T0:.3f}/T0, λ₊(n=3, K=400) = {lam_plus(3, Estar, A, 400)*T0:.2f}/T0")
        print(f"{'':>42} {'⟨ε⟩':>6} {'P(>8)':>7} {'KS D':>6} {'n lúc phát TB':>13}")
        print(f"{'dữ liệu':>42} {md:6.2f} {td:7.3f} {'—':>6}")
        rows = []
        for name, sq in engines.items():
            out = [cascade(sq[i % len(sq)], Estar, Sn, A, rng, nmax) for i in range(20000)]
            rows.append((name, [o for o in out if o]))
        for K in (135, 400):
            rows.append((f"exciton chuẩn (Williams, K = {K} MeV³)", [cascade_std(Estar, Sn, A, K, rng, nmax) for _ in range(20000)]))
        rows.append(("cân bằng tức thì (n = n_eq)", [((Estar - Sn) * rng.beta(2, nmax - 1), nmax) for _ in range(20000)]))
        for name, out in rows:
            x = np.array([o[0] for o in out]); nn = np.array([o[1] for o in out])
            mm, tt, cm, _ = X.stats_from_samples(x, 1.0, hi, eg)
            print(f"{name:>42} {mm:6.2f} {tt:7.3f} {np.max(np.abs(cm - cd)):6.3f} {nn.mean():13.2f}")
        print()


SN2 = {"Cu": 10.56, "Ge": 10.0, "Ni(n": 11.9, "60Ni": 10.0}     # ngưỡng tách neutron thứ hai (trung bình đồng vị)


def part_c4(N=24, Tmax=2000.0, dt=0.01):
    import giai_so_a9_exfor as X
    rng = np.random.default_rng(31)
    Sn_emit = {"Cu": 7.6, "Ge": 7.4, "Ni(n": 8.6, "60Ni": 11.7}
    engines = {}
    with Pool() as pool:
        for G0 in (2.0, 3.0):
            res = pool.map(_traj, [(0.3, 0.5, G0, 9700 + i, Tmax, dt) for i in range(N)])
            st = burst_stats(res, Tmax / 6)
            engines[G0] = sequences(res, Tmax / 6)
            print(f"Γ0 = {G0}, Γc = 0.5: {st['n']} đợt, va chạm/đợt {st['size']:.2f}, lưu {st['dwell']:.2f} T0 "
                  f"(CV {st['cv']:.2f}), ΔE/va chạm {st['eJ']:.3f}, ⟨|ℓ|⟩ khởi phát {st['lon']:.2f}", flush=True)
    print("=== C4. Phổ neutron có phát xạ lần hai (n,2n) — động cơ bơm mạnh vs exciton chuẩn vs dữ liệu ===")
    print(f"{'bộ dữ liệu':>36} {'mô hình':>26} {'⟨ε⟩':>6} {'P(>8)':>7} {'KS D':>6} {'bội neutron':>11}")
    for label, A, Einc, hi, es, s_, S in X.datasets():
        Estar = Einc * (A - 1) / A + S
        Sn = next(v for k, v in Sn_emit.items() if label.startswith(k))
        Sn2 = next(v for k, v in SN2.items() if label.startswith(k))
        nmax = int(round(math.sqrt(1.39 * (A / 13) * Estar)))
        md, td, eg, cd = X.stats_from_spec(es, s_, 1.0, hi)
        print(f"{label:>36} {'dữ liệu':>26} {md:6.2f} {td:7.3f} {'—':>6}")
        rows = []
        for G0, sq in engines.items():
            first = [cascade(sq[i % len(sq)], Estar, Sn, A, rng, nmax) for i in range(20000)]
            rows.append((f"động cơ Γ0 = {G0} + (n,2n)", first))
        rows.append(("exciton chuẩn K=400 + (n,2n)", [cascade_std(Estar, Sn, A, 400, rng, nmax) for _ in range(20000)]))
        for name, first in rows:
            x = with_second(first, Estar, Sn, Sn2, A, rng)
            mult = len(x) / sum(1 for o in first if o)
            mm, tt, cm, _ = X.stats_from_samples(x, 1.0, hi, eg)
            print(f"{'':>36} {name:>26} {mm:6.2f} {tt:7.3f} {np.max(np.abs(cm - cd)):6.3f} {mult:11.2f}", flush=True)
    print()


if __name__ == "__main__":
    args = sys.argv[1:] or ["C3", "C12"]
    best = None
    for a in args:
        if a.upper() == "C3":
            best = part_c3()
        elif a.upper() == "C4":
            part_c4()
        else:
            part_c12(*(best[1:3] if best else (None, None)))
