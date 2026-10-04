"""Giải số hai đề xuất:
(A) "Thang giữ pha tự thân" τ_E = ħ/δE với δE = 2%·ΔE. Đo độ gợn THỰC của động cơ ở bộ thông số hạt nhân
    (I=1, ΔE=1, κ=1/3, R=3, g=0.1; 1 đơn vị năng lượng = ħω0 = 10 MeV) khi quay đồng bộ ngay dưới ngưỡng gãy:
    gợn vận tốc δθ̇/θ̇, gợn động năng động cơ δK, và τ = ħ/δK nếu áp đúng công thức đề xuất.
    Kiểm tính "xyclic": p_θ = ∂L/∂θ̇ có bảo toàn không.
(B) Hàm tiêu tán Rayleigh Φ = ½γx·x'² + ½γθ·θ'², γ = I/τ (và m/τ): mô phỏng chu kỳ gãy/khóa có cản,
    đo thời gian đồng bộ cuối và động cơ còn quay (vượt rào 2ΔE) hay đã tắt.

Chạy: python3 giai_so_giu_pha_rayleigh.py [A|B]   (B ≈ 10 phút)
"""
import math
import sys

from giai_so_bo_thong_so_hat_nhan import HBAR_MEV_S, HW0
from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params
from giai_so_lien_ket_e_gay import constrained_acc, rk4_step
from giai_so_va_cham_tai_lien_ket import simulate

T0_S = 2 * math.pi * HBAR_MEV_S / HW0
BASE, R = params(1.0, 3.0, I=1.0, g=0.1)          # I=1, m=1, k=1, g=0.1, ΔE=1, κ=1/3
WC = omega_crit_general(BASE, R)


def part_a():
    print("=== A. Độ gợn thực của động cơ khi quay đồng bộ (bộ thông số hạt nhân, g=0.1) ===")
    print(f"θ̇_gãy = {WC:.5f}")
    f_c = lambda s: (s[1], *constrained_acc(s, **BASE)[:1], s[3], constrained_acc(s, **BASE)[1])
    for frac in (0.95, 0.8, 0.6):
        w0 = frac * WC
        s, t, dt = [0.0, w0, 0.0, BASE['kappa'] * w0], 0.0, 1e-4
        ws, ks, ps = [], [], []
        while t < 4 * 2 * math.pi:
            s = rk4_step(f_c, s, dt)
            t += dt
            th, w, x, v = s
            ws.append(w)
            ks.append(0.5 * BASE['I'] * w * w)                       # động năng riêng động cơ
            ps.append(BASE['I'] * w)                                   # p_θ = ∂L/∂θ̇ = Iθ̇
        wmax, wmin = max(ws), min(ws)
        dK = (max(ks) - min(ks))
        rip = (wmax - wmin) / (wmax + wmin)
        tau = HBAR_MEV_S / (dK * HW0)
        print(f"θ̇0 = {frac:.2f}·θ̇_gãy = {w0:.4f}: θ̇ ∈ [{wmin:.4f}, {wmax:.4f}] → δθ̇/θ̇ = {rip:.2%};"
              f"  p_θ dao động {min(ps):.4f}…{max(ps):.4f} (KHÔNG bảo toàn);"
              f"  δK = {dK:.4f} = {dK*HW0:.2f} MeV ⟹ ħ/δK = {tau:.3e} s = {tau/T0_S:.3f} T0")
    print(f"Đề xuất: δE = 0.02·10 MeV = 0.2 MeV ⟹ ħ/δE = {HBAR_MEV_S/0.2:.3e} s = {HBAR_MEV_S/0.2/T0_S:.2f} T0")
    print(f"Tham chiếu: 2% trong tài liệu đầu là δθ̇/θ̇ của bộ I=1, ΔE=0.5, Ω=5: ΔE/(IΩ²) = {0.5/25:.3f}\n")


def part_b():
    tau_list = [("τ_relax", 1.267), ("ħ/Γ1MeV", 10.0), ("τ_Pauli", 17.344), ("τ_E đề xuất 8T0", 8 * 2 * math.pi)]
    T = 10 * 2 * math.pi
    print(f"=== B. Chu kỳ gãy/khóa có cản Rayleigh (g=0.1, cửa sổ {T:.2f} = 10T0, θ̇_gãy={WC:.5f}) ===")
    print(f"{'chế độ cản':>22} {'τ':>7} {'f':>5} {'e':>4} {'gãy':>4} {'va chạm':>8} {'t đồng bộ cuối':>15} {'/T0':>6} "
          f"{'θ̇(T)':>9} {'động cơ':>10} {'%E mất':>7}")
    modes = [(f"γθ=I/τ, γx=m/τ", t, True, True) for _, t in tau_list] + [("chỉ γx=m/τ", 10.0, False, True)]
    for name, tau, on_th, on_x in modes:
        p = dict(BASE, R=R, gth=(BASE['I'] / tau) if on_th else 0.0, gx=(BASE['m'] / tau) if on_x else 0.0)
        for f in (1.05, 1.5, 2.0):
            for e in (0.0, 0.2):
                w0 = f * WC
                E0 = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * w0**2 - p['dE']
                r = simulate(p, w0, t_end=T, e=e)
                ev = r['events']
                nb = sum(x['kind'] == 'gãy' for x in ev)
                ni = len(ev) - nb
                tf = ev[-1]['t'] if (r['locked'] and ev and ev[-1]['kind'] == 'khóa') else float('nan')
                # động cơ còn quay nếu năng lượng quay vượt đỉnh rào: ½Iθ̇² − ΔEcosθ > ΔE
                th, w = r['s'][0], r['s'][1]
                rot = 0.5 * p['I'] * w * w - p['dE'] * math.cos(th) > p['dE']
                print(f"{name:>22} {tau:7.3f} {f:5.2f} {e:4.1f} {nb:4d} {ni:8d} {tf:15.4f} {tf/(2*math.pi):6.2f} "
                      f"{w:9.4f} {'còn quay' if rot else 'TẮT':>10} {(E0-r['E'])/E0:7.1%}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
