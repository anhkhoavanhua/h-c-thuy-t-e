"""Giải số đề xuất:
(A) τ_step = 2R/v_F ≈ 1.2e-22 s "trùng" khoảng cách xung 1.47e-22 s; so theo A; năng lượng xung so với
    năng lượng tách nucleon (điều kiện phát nucleon tiền cân bằng).
(B) Gợn tất định: lề ứng suất M_det = max|λ|/R mỗi vòng ở trạng thái dừng KHÔNG nhiễu, theo f;
    so với đề xuất 1 − Δλ/R ≈ (θ̇_d/θ̇_gãy)² = 0.673 tại f = 1.05.
(C) "Bão hòa hình học khi k_BT ≥ 1 MeV": giảm nhiễu xuống 0.01–0.3 MeV, đo tốc độ gãy và τ_lưu.

Đơn vị m=ω0=ħ=ℓ0=1, ħω0 = 10 MeV, τ=10, γx=0.70, g=0.1.
Chạy: python3 giai_so_saturation.py [A|B|C]   (C ≈ 6 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import analyse, cfg_for
from giai_so_kramers_A import setup
from giai_so_lien_ket_e_gay import rk4_step

TAU, GX, HW0 = 10.0, 0.70, 10.0
T0 = L.T0
HB = 6.582119569e-22
C_FM_S = 2.99792458e23


def part_a():
    print("=== A. Thang thời gian và năng lượng ===")
    ts = 10.0 / (0.27 * C_FM_S)
    print(f"τ_step = 10 fm/(0.27c) = {ts:.3e} s; khoảng cách xung mô hình 1.47e-22 s ⟹ tỉ số {1.47e-22/ts:.2f}")
    vF = 0.275
    for A in (16, 56, 208):
        hw = 41 * A ** (-1 / 3)
        model = 0.355 * 2 * math.pi * HB / hw
        cross = 2 * 1.2 * A ** (1 / 3) / (vF * C_FM_S)
        print(f"A={A:3d}: khoảng cách xung mô hình (0.355 T0) = {model:.2e} s; 2R/v_F = {cross:.2e} s; tỉ số {model/cross:.2f}")
    print("Năng lượng xung TB 4.0–4.6 MeV so với năng lượng tách nucleon S_n ≈ 6–10 MeV (hạt nhân bền): "
          "một xung KHÔNG đủ đưa nucleon ra khỏi hạt nhân")
    print("⟨N⟩ = 1/p: p = 0.195 được tính TỪ ⟨N⟩ = 5.12 (p = 1/5.12), nên 1/p = 5.13 là đồng nhất thức\n")


def _det_margin(a):
    cfg, f = a
    setup(cfg)
    wc = cfg[2]
    gam, Om = 1 / TAU, f * wc
    s = [0.0, 0.82 * wc, 0.0, L.BASE['kappa'] * 0.82 * wc]
    dt, t, Mmax, ws = 2e-3, 0.0, 0.0, []
    while t < 30 * T0:
        s = rk4_step(lambda z: L.f_con(z, gam, Om, GX), s, dt)
        t += dt
        if t > 25 * T0:
            Mmax = max(Mmax, abs(L.acc_con(s, gam, Om, GX)[2]) / L.R)
            ws.append(s[1])
    return f, Mmax, min(ws), max(ws), sum(ws) / len(ws)


def part_b():
    cfg = cfg_for(0.1)
    wc = cfg[2]
    x = 1 - GX * (1 / 3)**2 / (2 / TAU)
    print(f"=== B. Lề ứng suất tất định (không nhiễu, không cho gãy), θ̇_gãy = {wc:.4f} ===")
    print(f"{'f':>5} {'M_det = max|λ|/R':>17} {'đề xuất (θ̇_d/θ̇_gãy)²':>22} {'θ̇ min–max':>15} {'θ̇ TB':>7}")
    with Pool() as pool:
        for f, M, wmin, wmax, wm in pool.map(_det_margin, [(cfg, f) for f in (0.95, 1.0, 1.05, 1.1, 1.2)]):
            print(f"{f:5.2f} {M:17.4f} {(f*math.sqrt(x))**2:22.4f} {wmin:7.3f}–{wmax:<7.3f} {wm:7.3f}")
    print(f"Tiệm cận R ≈ mκθ̇_gãy² = {(1/3)*wc**2:.3f} so với R thật = 3 (lệch {((1/3)*wc**2-3)/3:+.1%})\n")


def _run(a):
    cfg, kT, seed, T = a
    setup(cfg)
    wc = cfg[2]
    return L.run(1 / TAU, 1.05 * wc, kT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc)['events']


def part_c(N=8):
    cfg = cfg_for(0.1)
    T, warm = 70 * T0, 10 * T0
    print(f"=== C. Tốc độ gãy khi giảm nhiễu (f=1.05, {N} quỹ đạo × 60 T0 sau khởi động) ===")
    print(f"{'k_BT':>9} {'tổng gãy':>9} {'pha dài':>8} {'đợt/T0':>8} {'τ_lưu (T0)':>11}")
    with Pool() as pool:
        for kT in (0.001, 0.003, 0.01, 0.03, 0.1):
            evs = pool.map(_run, [(cfg, kT, 25000 + i, T) for i in range(N)])
            nb = sum(sum(e[0] == 'gãy' and e[1] > warm for e in ev) for ev in evs)
            dw = [w for ev in evs for w in analyse(ev, warm)[0]]
            rate = len(dw) / (N * (T - warm) / T0)
            tl = f"{sum(dw)/len(dw)/T0:11.3f}" if dw else f"{'—':>11}"
            print(f"{kT*HW0:6.2f}MeV {nb:9d} {len(dw):8d} {rate:8.4f} {tl}", flush=True)
    print("(tham chiếu f=1.05 ở tài liệu trước, 1 MeV: 77 pha dài/480 T0, τ_lưu 4.57 T0)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
