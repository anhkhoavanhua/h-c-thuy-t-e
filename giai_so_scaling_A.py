"""Giải số ba đề xuất:
(A) Thứ nguyên của γθ = (2/3)α ħ²/(m_p c²) và dạng đúng (2/3)α ħ²ω0/(m_p c²); τ = I/γθ với I = ħ/ω0.
(B) Động cơ có duy trì được chu trình với τ ~ 10³–10⁴ khi có cản hạt γx = 0.789 không?
    Cân bằng công suất trong pha liên kết: lái γθ(1 − θ̇²/Ω0²)θ̇² so với tiêu tán ⟨γx ẋ²⟩ = ½γxκ²θ̇²
    ⟹ τ tới hạn τ_c ≈ 2I/(γx κ²). Mô phỏng tất định 50T0 với nhiều τ.
(C) "Dự đoán không tham số": ΔE = η·ħω0 với η lấy từ γx = 0.789 (tức ħω0 = 10 MeV). Nhưng γx = m/τ_relax
    tính theo đơn vị 1/ω0 phụ thuộc A: γx(A) = ħ/(τ_relax·ħω0(A)). Chạy Monte Carlo η cho từng hạt nhân
    với γx(A) của nó, rồi so với mức 3⁻ thấp nhất thực nghiệm (¹⁶O 6.130, ⁴⁰Ca 3.737, ²⁰⁸Pb 2.614 MeV).

Chạy: python3 giai_so_scaling_A.py [A|B|C]   (B ≈ 3–5 phút; C ≈ 5 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L

ALPHA = 1 / 137.035999
HBAR = 1.054571817e-34
HBAR_MEV_S = 6.582119569e-22
MEV = 1.602176634e-13
C = 2.99792458e8
MP_KG = 1.67262192e-27
MPC2 = 938.27208816
TAU_RELAX_S = 5e-15 / (0.2 * C)
HW0 = 10.0


def part_a():
    print("=== A. Thứ nguyên γθ ===")
    w0 = HW0 * MEV / HBAR
    I = HBAR / w0
    bad = (2 / 3) * ALPHA * HBAR**2 / (MP_KG * C**2)            # J·s² (sai)
    good = (2 / 3) * ALPHA * HBAR**2 * w0 / (MP_KG * C**2)      # J·s (đúng)
    print(f"Đề xuất (2/3)αħ²/(m_p c²) = {bad:.3e}  [đơn vị J·s² — không phải J·s]")
    print(f"Dạng đúng (2/3)αħ²ω0/(m_p c²) = {good:.3e} J·s;  I = ħ/ω0 = {I:.3e} kg·m²")
    print(f"τ = I/γθ (đơn vị 1/ω0) = I·ω0/γθ = {I*w0/good:,.0f}   (khớp 19 287 ở tài liệu trước)\n")


def run_det(args):
    tau, gx, T = args
    r = L.run(1.0 / tau, 1.05 * L.WC, 0.0, T, 2e-3, 0, gx=gx)
    br = [e[1] for e in r['events'] if e[0] == 'gãy']
    lk = [e for e in r['events'] if e[0] == 'khóa']
    th, w = r['s'][0], r['s'][1]
    rot = 0.5 * w * w - math.cos(th) > 1.0
    return tau, gx, len(br), sum(1 for b in br if b > T / 2), (br[-1] if br else float('nan')), w, rot, r['carnot'], len(lk)


def part_b():
    gx = 1.0 / 1.267
    kap = L.BASE['kappa']
    tau_c = 2 * L.BASE['I'] / (gx * kap * kap)
    print(f"=== B. Động cơ duy trì được không? γx = {gx:.3f}, κ = {kap:.4f} ⟹ τ_c ≈ 2I/(γxκ²) = {tau_c:.1f} ===")
    T = 50 * L.T0
    jobs = [(t, gx, T) for t in (10, 15, 20, 25, 30, 50, 100, 1000)] + [(1000, 0.0, T)]
    with Pool() as pool:
        res = pool.map(run_det, jobs)
    print(f"{'τ':>6} {'γx':>6} {'gãy tổng':>9} {'gãy nửa sau':>12} {'t gãy cuối':>11} {'θ̇(50T0)':>9} {'động cơ':>9} {'Carnot/khóa MeV':>16}")
    for tau, g, nb, nb2, tl, w, rot, car, nk in res:
        print(f"{tau:6g} {g:6.3f} {nb:9d} {nb2:12d} {tl:11.2f} {w:9.4f} {'còn quay' if rot else 'TẮT':>9} "
              f"{(car/nk*HW0 if nk else float('nan')):16.3f}")
    print()


def _mc(a):
    gx, seed, kT = a
    r = L.run(0.1, 1.05 * L.WC, kT, 10 * L.T0, 2e-3, seed, gx=gx)
    return r['carnot'], sum(e[0] == 'khóa' for e in r['events'])


def part_c(N=48):
    print(f"=== C. η(A) với γx(A) = ħ/(τ_relax·ħω0(A)), τ = 10, k_BT = 1 MeV, N = {N} ===")
    nuclei = [("¹⁶O", 16, 6.130), ("⁴⁰Ca", 40, 3.737), ("⁵⁶Fe", 56, None), ("²⁰⁸Pb", 208, 2.614)]
    print(f"{'hạt nhân':>8} {'ħω0':>7} {'γx(A)':>7} {'η đo':>7} {'ΔE mô hình':>11} {'ΔE công thức 15–17A^-1/3':>25} {'3⁻ thấp nhất':>13}")
    with Pool() as pool:
        for name, A, e3 in nuclei:
            hw = 41 * A ** (-1 / 3)
            gx = (HBAR_MEV_S / hw) / TAU_RELAX_S            # m/τ_relax theo đơn vị ω0 (m=1)
            res = pool.map(_mc, [(gx, 7000 + i, 0.1) for i in range(N)])
            eta = sum(r[0] for r in res) / max(sum(r[1] for r in res), 1)
            lo, hi = 15.2 * A ** (-1 / 3), 16.8 * A ** (-1 / 3)
            print(f"{name:>8} {hw:7.3f} {gx:7.3f} {eta:7.4f} {eta*hw:9.3f} MeV {lo:10.2f} – {hi:5.2f} MeV "
                  f"{(f'{e3:.3f} MeV' if e3 else '—'):>13}", flush=True)
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
