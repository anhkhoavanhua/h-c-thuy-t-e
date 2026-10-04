"""Nguyên lý toàn ảnh AdS/CFT cho học thuyết E_dual (bản thảo HocThuyetE (7).pdf, Bảng I: ħω0 = 10 MeV, A ~ 70,
E_vt = 30 MeV, ℓ0 = 2.037 fm, κ = ℓ0/3).

Cầu nối kiểm được: AdS toàn cục là "bẫy điều hòa" — metric g_00 = −(1 + r²/L²) cho thế Newton c²r²/(2L²), tức
ω0 = c/L. Hạt x của E_dual (dao động tử ħω0) = hạt khối lượng m_N trong AdS với L = ħc/ħω0.
(A) Từ điển thang: L, ℓ0 = √(ƛ_C L), E_vt theo ħc/L.
(B) Trắc địa kiểu thời gian trong AdS toàn cục (RK4): chu kỳ theo thời gian toàn cục = 2πL/c với MỌI biên độ
    (đẳng thời), so với T0 của bản thảo; chu kỳ thời gian riêng cũng = 2πL/c (mọi trắc địa hội tụ lại).
(C) Phổ trường vô hướng khối lượng m_N trong AdS2 (x một chiều như bản thảo) và AdS4 (hạt nhân 3 chiều),
    giải số phương trình Klein–Gordon toàn cục bằng sai phân hữu hạn: ω = Δ + n (AdS2), ω = Δ + 2n + l (AdS4),
    Δ(Δ − d) = m²L² (từ điển AdS/CFT); so E − m c² với dao động tử điều hòa và bậc suy biến vỏ.
(D) Tính tuần hoàn của thời gian: e^{−iΔτ} sau một chu kỳ AdS 2π; bó sóng hạt tự do lặp lại sau 2π.
(E) Cận Bekenstein cho đợt Carnot (E ≤ 20.4 MeV, R = 1.2 A^{1/3} fm) so với số trạng thái exciton Ericson.
Chạy: python3 giai_so_ads_cft.py   (≈ 10 s)
"""
import math

import numpy as np

HBARC = 197.3270          # MeV fm
MN = 938.919              # MeV, khối lượng nucleon trung bình
HW = 10.0                 # ħω0, MeV (Bảng I)
EVT = 30.0                # MeV
A = 70
L = HBARC / HW            # fm
MUL = MN / HW             # m L c/ħ (không thứ nguyên)


def part_a():
    lc = HBARC / MN
    l0 = math.sqrt(HBARC**2 / (MN * HW))
    T0 = 2 * math.pi * L * 1e-15 / 2.99792458e8
    print("=== A. Từ điển thang ===")
    print(f"L_AdS = ħc/ħω0 = {L:.3f} fm;  2πL/c = {T0:.4e} s (Bảng I: T0 = 4.136e-22 s)")
    print(f"ℓ0 = √(ħ/mω0) = {l0:.4f} fm;  √(ƛ_C·L) = √({lc:.4f}·{L:.3f}) = {math.sqrt(lc*L):.4f} fm")
    print(f"ℓ0/L = {l0/L:.4f} (v/c của hạt ở biên độ ℓ0);  κ/L = {l0/3/L:.4f}")
    print(f"E_vt = {EVT} MeV = {EVT*L/HBARC:.3f} ħc/L;  giả thuyết cũ E_vt = ħc/L_AdS cho L_AdS = {HBARC/EVT:.3f} fm "
          f"= L/{L/(HBARC/EVT):.2f}")
    print(f"R_hạt nhân = 1.2·A^(1/3) = {1.2*A**(1/3):.3f} fm = {1.2*A**(1/3)/L:.3f} L\n")


def part_b():
    print("=== B. Trắc địa kiểu thời gian xuyên tâm trong AdS toàn cục (tọa độ ρ, r = L tan ρ, τ = ct/L) ===")
    print(f"{'r_max/L':>8} {'r_max (fm)':>11} {'v_max/c':>8} {'chu kỳ τ (/2π)':>15} {'chu kỳ riêng (/2π)':>19} "
          f"{'E':>8}")
    for rmax in (0.0344, 0.1032, 0.25, 1.0, 10.0):          # κ/L, ℓ0/L, ..., tương đối tính mạnh
        E = math.sqrt(1 + rmax**2)                          # cos ρ_max = 1/E, tan ρ_max = r_max/L
        rho, v, tau, prop, h = 0.0, math.sqrt(1 - 1 / E**2), 0.0, 0.0, 1e-4
        acc = lambda r: -math.tan(r) / (math.cos(r)**2 * E**2)
        cross = []
        while len(cross) < 2:
            r0 = rho
            k1r, k1v = v, acc(rho)
            k2r, k2v = v + h / 2 * k1v, acc(rho + h / 2 * k1r)
            k3r, k3v = v + h / 2 * k2v, acc(rho + h / 2 * k2r)
            k4r, k4v = v + h * k3v, acc(rho + h * k3r)
            rho += h / 6 * (k1r + 2 * k2r + 2 * k3r + k4r)
            v += h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
            prop += h * 1 / (E * math.cos(r0)**2)            # dλ/dτ = 1/(E cos²ρ) (λ theo đơn vị L/c)
            tau += h
            if r0 > 0 and rho <= 0 or r0 < 0 and rho >= 0:
                cross.append(tau - h * rho / (rho - r0))
        vmax = math.sqrt(1 - 1 / E**2)                      # vận tốc riêng cục bộ tại ρ = 0
        print(f"{rmax:8.4f} {rmax*L:11.3f} {vmax:8.4f} {cross[1]/(2*math.pi):15.8f} "
              f"{prop/(2*math.pi)*cross[1]/tau:19.6f} {E:8.4f}")
    print("(cả hai chu kỳ = 2πL/c cho mọi biên độ, kể cả v → c; hai hàng đầu là biên độ κ và ℓ0 của bản thảo)\n")


def spectrum(d, l, N=2400, rmax=0.8):
    nu = math.sqrt(d * d / 4 + MUL**2)
    if d == 1:
        x = np.linspace(-rmax, rmax, N + 2)[1:-1]
        V = (nu**2 - 0.25) / np.cos(x)**2
    else:
        x = np.linspace(0, rmax, N + 2)[1:-1]
        V = (nu**2 - 0.25) / np.cos(x)**2 + ((l + (d - 2) / 2)**2 - 0.25) / np.sin(x)**2
    h = x[1] - x[0]
    M = np.diag(2 / h**2 + V) - np.diag(np.ones(N - 1) / h**2, 1) - np.diag(np.ones(N - 1) / h**2, -1)
    return np.sqrt(np.linalg.eigvalsh(M)[:6]), d / 2 + nu


def part_c():
    print(f"=== C. Phổ Klein–Gordon toàn cục, m_N L c/ħ = {MUL:.3f} ===")
    w1, D1 = spectrum(1, 0)
    w1b, _ = spectrum(1, 0, N=4800)
    w1r = (4 * w1b - w1) / 3                                 # ngoại suy Richardson O(h²)
    print(f"AdS2/CFT1: Δ = 1/2 + √(1/4 + m²L²) = {D1:.6f}")
    print(f"{'n':>3} {'ω số':>12} {'Δ + n':>12} {'sai':>9} {'E − mc² (MeV)':>14} {'ħω0(n+½)':>10}")
    for n in range(5):
        print(f"{n:3d} {w1r[n]:12.6f} {D1+n:12.6f} {w1r[n]-D1-n:9.1e} {(w1r[n]*HW-MN):14.4f} {HW*(n+0.5):10.4f}")
    print("AdS4/CFT3 (hạt nhân 3 chiều):")
    levels = {}
    for l in range(0, 6):
        wa, D3 = spectrum(3, l)
        wb, _ = spectrum(3, l, N=4800)
        wr = (4 * wb - wa) / 3
        for n in range(3):
            levels.setdefault(2 * n + l, []).append((n, l, wr[n]))
    print(f"Δ = 3/2 + √(9/4 + m²L²) = {D3:.6f}; Δ − mL − 3/2 = {D3-MUL-1.5:.6f} "
          f"(dịch tương đối tính {HW*(D3-MUL-1.5)*1000:.1f} keV)")
    print(f"{'N=2n+l':>7} {'(n,l) có mặt':>22} {'max|ω − Δ − N|':>15} {'suy biến (2l+1)':>16} {'HO 3D (N+1)(N+2)/2':>19} "
          f"{'E − mc² (MeV)':>14}")
    for N in range(0, 5):
        lv = levels[N]
        err = max(abs(w - D3 - N) for _, _, w in lv)
        deg = sum(2 * l + 1 for _, l, _ in lv)
        print(f"{N:7d} {' '.join(f'({n},{l})' for n, l, _ in lv):>22} {err:15.1e} {deg:16d} {(N+1)*(N+2)//2:19d} "
              f"{lv[0][2]*HW-MN:14.4f}")
    shells = np.cumsum([2 * (N + 1) * (N + 2) // 2 for N in range(5)])
    print(f"số nucleon lấp đầy vỏ (×2 spin): {list(shells)}  ⟹ A ~ 70 ≈ vỏ N = 4 (bản thảo: A ~ 70)\n")
    return D1


def part_d(D1):
    print("=== D. Tính tuần hoàn của thời gian toàn cục ===")
    print(f"pha e^(−iΔτ) sau τ = 2π: Δ mod 1 = {D1 % 1:.6f}  ⟹ Δ không nguyên, trường chỉ xác định trên PHỦ PHỔ QUÁT "
          "(τ ∈ ℝ); trên hyperboloid (τ ~ τ + 2π) thì bị cấm")
    ns = np.arange(40)
    c = np.exp(-(ns - 6.0)**2 / 8.0) * np.exp(1j * 0.7 * ns)
    c /= np.linalg.norm(c)
    for tau in (math.pi, 2 * math.pi, 4 * math.pi):
        ov = abs(np.sum(np.abs(c)**2 * np.exp(-1j * (D1 + ns) * tau)))**2
        print(f"bó sóng (≈ 7 mức): |⟨ψ(0)|ψ(τ = {tau/math.pi:.0f}π)⟩|² = {ov:.6f}")
    print("⟹ trạng thái vật lý của hạt lặp lại sau mỗi chu kỳ 2π (mọi trắc địa hội tụ lại): hạt trong AdS là đồng hồ "
          "tuần hoàn, như OAM\n")


def part_e():
    R = 1.2 * A**(1 / 3)
    g = A / 13.0
    print("=== E. Cận Bekenstein S ≤ 2πER/(ħc) cho đợt Carnot ===")
    print(f"{'E (MeV)':>8} {'cận, chỉ năng lượng kích thích (nat)':>36} {'ln #exciton trong 1 MeV (nat)':>30} "
          f"{'cận, kể cả khối lượng nghỉ (nat)':>32}")
    for E in (8.0, 13.0, 20.4):
        best = max(n * math.log(g) + (n - 1) * math.log(E) - math.lgamma(p + 1) - math.lgamma(n - p + 1)
                   - math.lgamma(n) for n in range(1, 12) for p in range(0, n + 1))
        print(f"{E:8.1f} {2*math.pi*E*R/HBARC:36.2f} {best:30.2f} {2*math.pi*(E+A*MN)*R/HBARC:32.0f}")
    print(f"(R = {R:.2f} fm, g = A/13 = {g:.2f} MeV⁻¹)\n")


if __name__ == "__main__":
    part_a()
    part_b()
    D1 = part_c()
    part_d(D1)
    part_e()
