"""Giải số đề xuất:
(1) τ_lưu ≈ (3–4)·τ_rel do khoảng cách (θ̇_gãy − θ̇_d)/θ̇_gãy = 1 − 0.782f; dự đoán: f tăng ⟹ τ_lưu giảm mạnh,
    f < 1 ⟹ τ_lưu → ∞. Đo τ_lưu theo f = Ω0/θ̇_gãy.
(2) Đợt vi chấn = chuỗi gamma cascade: bội số hình học/Poisson, E_γ ≈ 3–5 MeV/photon.
    So thời khoảng giữa các xung với thời gian sống Weisskopf của chuyển dời gamma; kiểm phân bố bội số.

Đơn vị m=ω0=ħ=ℓ0=1, ħω0 = 10 MeV, τ=10, γx=0.70, g=0.1, k_BT = 1 MeV.
Chạy: python3 giai_so_f_cascade.py [A|B|C]   (B ≈ 6 phút, C ≈ 3 phút; 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import analyse, cfg_for
from giai_so_kramers_A import setup

TAU, GX, KT, HW0 = 10.0, 0.70, 0.1, 10.0
T0 = L.T0
T0_S = 2 * math.pi * 6.582119569e-22 / HW0
KAP = 1 / 3


def part_a():
    print("=== A. Thời gian sống Weisskopf so với khoảng cách xung trong đợt ===")
    hb = 6.582119569e-22
    for A in (16, 56, 208):
        for E in (3.0, 4.0, 5.0):
            e1 = 1.0e14 * A ** (2 / 3) * E**3        # s⁻¹ (E theo MeV)
            m1 = 5.6e13 * E**3
            e2 = 7.3e7 * A ** (4 / 3) * E**5
            print(f"A={A:3d}, E_γ={E} MeV: τ(E1) = {1/e1:.1e} s, τ(M1) = {1/m1:.1e} s, τ(E2) = {1/e2:.1e} s")
    print(f"Khoảng cách xung trong đợt của mô hình: < 0.2 T0 = {0.2*T0_S:.1e} s  (T0 = {T0_S:.3e} s)\n")


def _run(a):
    cfg, f, seed, T = a
    setup(cfg)
    wc = cfg[2]
    return L.run(1 / TAU, f * wc, KT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc)['events']


def part_b(N=8):
    cfg = cfg_for(0.1)
    wc = cfg[2]
    I_eff = 1 + 0.5 * KAP**2
    x = 1 - GX * KAP**2 / (2 / TAU)
    t_rel = I_eff / (2 * (1 / TAU) * x) / T0            # |F'| = 2γθ θ̇_d²/Ω0² = 2γθ·x — không phụ thuộc f
    print(f"=== B. τ_lưu theo f = Ω0/θ̇_gãy (θ̇_gãy = {wc:.4f}, τ_rel = {t_rel:.3f} T0 cho MỌI f) ===")
    print(f"{'f':>5} {'θ̇_d/θ̇_gãy':>11} {'khoảng cách':>11} {'gãy':>5} {'pha dài':>8} {'τ_lưu (T0)':>11} {'τ_lưu/τ_rel':>12} {'CV dư':>6}")
    T, warm = 70 * T0, 10 * T0
    with Pool() as pool:
        for f in (0.90, 0.95, 1.00, 1.05, 1.10, 1.20):
            evs = pool.map(_run, [(cfg, f, 21000 + i, T) for i in range(N)])
            dw = [w for ev in evs for w in analyse(ev, warm)[0]]
            nb = sum(sum(e[0] == 'gãy' and e[1] > warm for e in ev) for ev in evs)
            ratio = f * math.sqrt(x)
            if len(dw) > 1:
                m = sum(dw) / len(dw) / T0
                res = [w - T0 for w in dw]
                mr = sum(res) / len(res)
                cv = math.sqrt(sum((r - mr)**2 for r in res) / (len(res) - 1)) / mr
                print(f"{f:5.2f} {ratio:11.4f} {1-ratio:11.4f} {nb:5d} {len(dw):8d} {m:11.3f} {m/t_rel:12.2f} {cv:6.2f}", flush=True)
            else:
                print(f"{f:5.2f} {ratio:11.4f} {1-ratio:11.4f} {nb:5d} {len(dw):8d} {'—':>11} {'—':>12} {'—':>6}", flush=True)
    print()


def part_c(N=16):
    cfg = cfg_for(0.1)
    T, warm = 70 * T0, 10 * T0
    print(f"=== C. Bội số xung mỗi đợt (f = 1.05, {N} quỹ đạo × 60 T0) so với phân bố hình học ===")
    with Pool() as pool:
        evs = pool.map(_run, [(cfg, 1.05, 23000 + i, T) for i in range(N)])
    bs = [b for ev in evs for b in analyse(ev, warm)[1]]
    ns = [b['n'] for b in bs]
    n = len(ns)
    mean = sum(ns) / n
    p = 1 / mean
    print(f"{n} đợt, bội số TB {mean:.2f} ⟹ hình học p = {p:.3f}")
    print(f"{'N':>4} {'đếm':>5} {'kỳ vọng hình học':>17} {'kỳ vọng Poisson(TB−1)+1':>24}")
    chi_g = chi_p = 0.0
    lam = mean - 1
    for k in range(1, 13):
        obs = sum(1 for v in ns if v == k) if k < 12 else sum(1 for v in ns if v >= 12)
        eg = n * (p * (1 - p)**(k - 1) if k < 12 else (1 - p)**11)
        ep = n * (math.exp(-lam) * lam**(k - 1) / math.factorial(k - 1) if k < 12
                  else 1 - sum(math.exp(-lam) * lam**j / math.factorial(j) for j in range(11)))
        chi_g += (obs - eg)**2 / eg
        chi_p += (obs - ep)**2 / max(ep, 1e-9)
        print(f"{('≥' if k == 12 else '')+str(k):>4} {obs:5d} {eg:17.1f} {ep:24.1f}")
    print(f"χ² (12 ô): hình học {chi_g:.1f}, Poisson dịch {chi_p:.1f}  (ngưỡng 5% với 10–11 bậc tự do ≈ 18–20)")
    gaps = []
    for b in bs:
        if b['n'] >= 2:
            gaps.append(b['dur'] / (b['n'] - 1) / T0)
    print(f"Khoảng cách TB giữa hai lần khóa trong đợt (đợt ≥2 khóa): {sum(gaps)/len(gaps):.3f} T0 = "
          f"{sum(gaps)/len(gaps)*T0_S:.2e} s\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
