"""Giải số ba đề xuất:
(A) γθ từ cản bức xạ: kiểm thứ nguyên của γθ = (2/3)α_eff ħω0²/c², tính tốc độ tắt bức xạ chuẩn
    Γ = (2/3)α ħω0²/(mc²) và τ = ω0/Γ (đơn vị 1/ω0) cho α_eff = 1/137, 0.01, 0.1.
(B) Langevin đầy đủ 2 bậc tự do: ξθ (γθ = I/τ) + ξx (γx = m/τ_relax = 0.789), k_BT = 1 hoặc 5 MeV.
    Monte Carlo 64 quỹ đạo: tái khóa, thời gian tái khóa, năng lượng phát xạ mỗi lần khóa, τ_deph đo lại.
(C) "7.6–7.9 MeV/xung = B/A": năng lượng xung tỉ lệ ħω0 = 41A^(-1/3) MeV; so với B/A thực nghiệm.

Chạy: python3 giai_so_langevin_day_du.py [A|B|C]   (B ≈ 5–8 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L

ALPHA = 1 / 137.035999
MPC2 = 938.27208816
HW0 = 10.0


def part_a():
    print("=== A. Cản bức xạ ===")
    print("Thứ nguyên ħω0²/c²: [J·s]·[s⁻²]/[m²·s⁻²] = kg/s; γθ cần [I]/[t] = kg·m²/s ⟹ THIẾU [m²].")
    for a in (ALPHA, 0.01, 0.1):
        G_over_w = (2 / 3) * a * HW0 / MPC2        # Γ/ω0 = (2/3) α ħω0/(mc²)
        print(f"α_eff = {a:.5f}: Γ/ω0 = {G_over_w:.3e} ⟹ τ = ω0/Γ = {1/G_over_w:,.0f} đơn vị 1/ω0 "
              f"= {1/G_over_w/(2*math.pi):,.0f} T0   (đề xuất: τ ≈ 10–50)")
    print()


def _worker(a):
    gam, Om, seed, T, dt, kT, gx = a
    r = L.run(gam, Om, kT, T, dt, seed, gx=gx)
    ev = r['events']
    b = next((e[1] for e in ev if e[0] == 'gãy'), None)
    k = next((e[1] for e in ev if e[0] == 'khóa' and b is not None and e[1] > b), None)
    nk = sum(e[0] == 'khóa' for e in ev)
    return b, k, sum(e[0] == 'gãy' for e in ev), nk, r['carnot']


def _phase(a):
    gam, Om, seed, T, dt, every, kT, gx = a
    return L.run(gam, Om, kT, T, dt, seed, breakable=False, sample_every=every, gx=gx)['samples']


def tau_deph(pool, gam, kT, gx, N=64):
    T, dt, every = 10 * L.T0, 2e-3, 0.5
    Om = 0.9 * L.WC
    ref = L.run(gam, Om, 0.0, T, dt, 0, breakable=False, sample_every=every, gx=gx)['samples']
    sims = pool.map(_phase, [(gam, Om, 1000 + i, T, dt, every, kT, gx) for i in range(N)])
    n = min(len(ref), *(len(x) for x in sims))
    var = [sum((x[j] - ref[j])**2 for x in sims) / N for j in range(n)]
    return next((j * every for j in range(n) if var[j] >= 1.0), float('inf'))


def part_b(N=64):
    gx_relax = 1.0 / 1.267
    print(f"=== B. Langevin đầy đủ (ξθ + ξx), Ω0 = 1.05θ̇_gãy, N = {N}, cửa sổ 10T0 ===")
    print(f"{'τ':>6} {'γx':>6} {'k_BT':>6} {'gãy≥1':>6} {'tái khóa':>9} {'t tái khóa min/TB/max':>24} "
          f"{'τ_deph đo':>10} {'kịp':>7} {'gãy/qđ':>7} {'Carnot/khóa (MeV)':>18}")
    T, dt = 10 * L.T0, 2e-3
    with Pool() as pool:
        for tau, gx, kT in [(10.0, 0.0, 0.5), (10.0, gx_relax, 0.1), (10.0, gx_relax, 0.5), (1000.0, gx_relax, 0.5)]:
            gam = 1.0 / tau
            Om = 1.05 * L.WC
            res = pool.map(_worker, [(gam, Om, 5000 + i, T, dt, kT, gx) for i in range(N)])
            broke = [r for r in res if r[0] is not None]
            rel = [r[1] - r[0] for r in broke if r[1] is not None]
            td = tau_deph(pool, gam, kT, gx)
            ok = sum(1 for d in rel if d <= td)
            locks = sum(r[3] for r in res)
            emis = sum(r[4] for r in res) / max(locks, 1) * HW0
            nb = sum(r[2] for r in res) / N
            trel = f"{min(rel):.3f}/{sum(rel)/len(rel):.3f}/{max(rel):.3f}" if rel else "—"
            print(f"{tau:6g} {gx:6.3f} {kT*10:4g}MeV {len(broke):3d}/{N} {len(rel):4d}/{len(broke):<3d} {trel:>24} "
                  f"{td:10.2f} {ok:3d}/{len(broke):<3d} {nb:7.1f} {emis:18.3f}", flush=True)
    print()


def part_c(frac=(0.764, 0.793)):
    print("=== C. Năng lượng xung Carnot = (0.764…0.793)·ħω0 so với B/A thực nghiệm ===")
    BA = {16: ("¹⁶O", 7.976), 40: ("⁴⁰Ca", 8.551), 56: ("⁵⁶Fe", 8.790), 120: ("¹²⁰Sn", 8.505), 208: ("²⁰⁸Pb", 7.867)}
    print(f"{'hạt nhân':>8} {'ħω0 = 41A^-1/3':>15} {'xung mô hình (MeV)':>20} {'B/A (MeV)':>10}")
    for A, (name, ba) in BA.items():
        hw = 41 * A ** (-1 / 3)
        print(f"{name:>8} {hw:15.3f} {frac[0]*hw:9.2f} – {frac[1]*hw:5.2f} {ba:10.3f}")
    A_match = (41 / 10.0) ** 3
    print(f"Chỉ khớp 7.6–7.9 MeV khi ħω0 ≈ 10 MeV ⟺ A ≈ (41/10)³ = {A_match:.0f}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "C", "B"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
