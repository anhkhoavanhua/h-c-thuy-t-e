"""Giải số đề xuất "tốc độ gãy ≈ 0.15/T0 là chu kỳ nạp Van der Pol (dao động tích thoát)".

(A) Phương trình trung bình đề xuất: I_eff·dθ̇/dt = F(θ̇) = γθ(1 − θ̇²/Ω0²)θ̇ − ½γxκ²θ̇.
    Điểm dừng θ̇_d, hệ số hồi phục F'(θ̇_d) = −2γθθ̇_d²/Ω0², thời gian hồi phục τ_rel = I_eff/|F'|;
    kiểm dấu F giữa θ̇_d và θ̇_gãy (tích phân nạp tới ngưỡng có tồn tại không); số học T_cycle của đề xuất.
(B) Đo trực tiếp (Langevin, k_BT = 1 MeV) theo lưới τ, γx: tốc độ đợt thoát r, thời gian từ khóa tới
    lần gãy kế tiếp, θ̇ ngay sau khóa — so với τ_rel và với công thức đề xuất.
(C) Kiểm "Carnot là cận trên bao trùm các mức 2⁺, 3⁻ thấp nhất" với tiên đề tự đồng dạng E_vt/ħω0 = 3.

Chạy: python3 giai_so_recharge.py [A|B|C]   (B ≈ 5 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_kramers_A import setup, std_cfg

KAP = 1 / 3
I_EFF = 1 + 0.5 * KAP**2


def averaged(tau, gx, wc):
    gth = 1.0 / tau
    Om = 1.05 * wc
    x = 1 - gx * KAP**2 / (2 * gth)
    wd = Om * math.sqrt(x) if x > 0 else 0.0
    F = lambda w: gth * (1 - w * w / Om**2) * w - 0.5 * gx * KAP**2 * w
    fp = -2 * gth * wd**2 / Om**2
    return Om, wd, F, fp


def part_a():
    cfg = std_cfg()
    wc = cfg[2]
    print(f"=== A. Phương trình trung bình (θ̇_gãy = {wc:.4f}, I_eff = {I_EFF:.4f}) ===")
    for tau, gx in ((10, 0.70), (10, 0.789)):
        Om, wd, F, fp = averaged(tau, gx, wc)
        mid = 0.5 * (wd + wc)
        print(f"τ={tau}, γx={gx}: Ω0={Om:.4f}, θ̇_d={wd:.4f}, F(θ̇ giữa θ̇_d và θ̇_gãy = {mid:.3f}) = {F(mid):+.5f} "
              f"(<0 ⟹ không thể nạp tất định lên ngưỡng); F'(θ̇_d) = {fp:.5f} ⟹ τ_rel = {I_EFF/abs(fp):.3f} "
              f"= {I_EFF/abs(fp)/L.T0:.3f} T0")
    Om, wd = 1.05 * wc, 2.63
    lnr = math.log((Om - wd) / (Om - wc))
    print(f"Đề xuất: ln((Ω0−θ̇_d)/(Ω0−θ̇_gãy)) = ln({Om-wd:.3f}/{Om-wc:.3f}) = {lnr:.3f}; "
          f"τ_eff không được cho — muốn ra 6–7 đv cần τ_eff = {6/lnr:.2f}–{7/lnr:.2f}")
    print(f"Đề xuất: T_recharge 6–7 đv + t_ổn định 4–5 đv = 10–12 đv = {10/L.T0:.2f}–{12/L.T0:.2f} T0 "
          f"⟹ r = {L.T0/12:.3f}–{L.T0/10:.3f} /T0   (đề xuất viết '6.5–7.0 đv ≈ 6.7 T0'; 6.7 T0 = {6.7*L.T0:.1f} đv)\n")


def _worker(a):
    cfg, tau, gx, kT, seed, T, warm = a
    setup(cfg)
    wc = cfg[2]
    r = L.run(1.0 / tau, 1.05 * wc, kT, T, 2e-3, seed, gx=gx, w_start=0.82 * wc)
    ev = r['events']
    episodes, waits, w_after, last_lock = 0, [], [], None
    for e in ev:
        if e[0] == 'khóa':
            last_lock = e[1]
            if e[1] > warm:
                w_after.append(e[2])
        elif e[0] == 'gãy' and e[1] > warm and last_lock is not None and e[1] - last_lock >= L.T0:
            episodes += 1
            waits.append(e[1] - last_lock)
    return episodes, waits, w_after


def part_b(N=8):
    cfg = std_cfg()
    wc = cfg[2]
    T, warm, kT = 70 * L.T0, 10 * L.T0, 0.1
    print(f"=== B. Đo chu kỳ gãy theo τ, γx (k_BT = 1 MeV, {N} quỹ đạo × 60 T0 sau khởi động) ===")
    print(f"{'τ':>4} {'γx':>5} {'θ̇_d':>7} {'τ_rel (T0)':>11} {'r (/T0)':>8} {'1/r (T0)':>9} {'1/r ÷ τ_rel':>11} "
          f"{'chờ khóa→gãy TB (T0)':>21} {'θ̇ sau khóa TB':>15}")
    with Pool() as pool:
        for tau in (5.0, 10.0, 15.0):
            for gx in (0.5, 0.7):
                Om, wd, F, fp = averaged(tau, gx, wc)
                trel = I_EFF / abs(fp) / L.T0 if fp else float('inf')
                res = pool.map(_worker, [(cfg, tau, gx, kT, 11000 + i, T, warm) for i in range(N)])
                ne = sum(r[0] for r in res)
                waits = [w for r in res for w in r[1]]
                wa = [w for r in res for w in r[2]]
                rate = ne / (N * (T - warm) / L.T0)
                print(f"{tau:4g} {gx:5.2f} {wd:7.4f} {trel:11.3f} {rate:8.4f} {(1/rate if rate else float('inf')):9.2f} "
                      f"{((1/rate)/trel if rate else float('inf')):11.2f} "
                      f"{(sum(waits)/len(waits)/L.T0 if waits else float('nan')):21.2f} "
                      f"{(sum(wa)/len(wa) if wa else float('nan')):15.4f}", flush=True)
    print()


def part_c():
    print("=== C. 'Cận trên Carnot' (0.42–0.48)ħω0 so với mức 2⁺ và 3⁻ thấp nhất ===")
    data = [("¹⁶O", 16, 6.917, 6.130), ("⁴⁰Ca", 40, 3.904, 3.737), ("²⁰⁸Pb", 208, 4.085, 2.614)]
    print(f"{'hạt nhân':>8} {'ħω0':>7} {'E_vt=3ħω0':>10} {'L=ħc/E_vt':>10} {'1.2A^1/3':>9} {'Carnot':>14} {'2⁺':>7} {'3⁻':>7} {'bao trùm?':>10}")
    for name, A, e2, e3 in data:
        hw = 41 * A ** (-1 / 3)
        lo, hi = 0.42 * hw, 0.48 * hw
        cover = "có" if hi >= max(e2, e3) else "KHÔNG"
        print(f"{name:>8} {hw:7.3f} {3*hw:10.2f} {197.327/(3*hw):9.2f}f {1.2*A**(1/3):8.2f}f {lo:6.2f}–{hi:5.2f} "
              f"{e2:7.3f} {e3:7.3f} {cover:>10}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "C", "B"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
