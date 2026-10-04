"""Năm loại đồng hồ lượng tử áp vào động cơ E_dual (HocThuyetE (7).pdf; ħω0 = 10 MeV, T0 = 4.136·10⁻²² s, A ~ 70).
Đòi hỏi chung (dong-ho-hat-nhan-page-wootters-kiem-chung.md, N1): xếp đúng 99% ba pha ⟺ σ_t ≤ 0.071 T0; 90% ⟺ 0.27 T0;
một chu trình ≈ 5 T0.

1 Salecker–Wigner: M c² ≳ ħT/δt²; khối lượng tối thiểu; thời gian tối đa một nucleon / một hạt nhân giữ được độ chính xác.
2 WKB: rô-to là biến "nặng"? Tiêu chí |d(1/p_θ)/dθ| ≪ 1 dọc quỹ đạo ràng buộc (Eq. 6) ở θ̇_d = 0.821 θ̇_br; tác dụng J/ħ.
3 Brown–Kuchař (bụi không tương tác): quãng đường tự do trung bình nucleon λ = 1/(ρσ_eff) so với bán kính hạt nhân.
4 Entropy: entropy hệ S(ρ_S) (bão hòa, Part III N1) vs entropy bể tích lũy ΔS = Σ ΔE_Carnot/k_BT.
5 Đồng hồ tự hành (Erker et al. 2017) + TUR (độ chính xác N ≤ Σ_tick/2k_B): entropy và năng lượng cần tỏa mỗi T0.
Chạy: python3 giai_so_nam_dong_ho.py
"""
import math

import numpy as np

HBAR = 6.582119569e-22          # MeV s
T0 = 4.1357e-22
E_T0 = HBAR / T0                # ħ/T0 (MeV)
MN = 939.565


def sw():
    print("=== 1. Salecker–Wigner: M c² ≳ ħ T / δt² ===")
    print(f"ħ/T0 = {E_T0:.3f} MeV")
    print(f"{'δt (T0)':>8} {'T (T0)':>7} {'M c² tối thiểu':>15} {'so nucleon':>11} {'so hạt nhân A=70':>17}")
    for dt in (0.27, 0.071):
        for T in (5, 50, 300):
            M = E_T0 * T / dt**2
            print(f"{dt:8.3f} {T:7d} {M:12.0f} MeV {M/MN:11.1f} {M/(70*MN):17.2f}")
    for dt in (0.27, 0.071):
        print(f"một nucleon (939 MeV) giữ δt = {dt} T0 tối đa T = {MN*dt**2/E_T0:.1f} T0; "
              f"cả hạt nhân A = 70: {70*MN*dt**2/E_T0:.0f} T0")
    print()


def wkb():
    import giai_so_schrodinger_rang_buoc as S
    th = np.linspace(0, 2 * math.pi, 4001)
    wd = 0.821 * S.WBR
    E = 0.5 * S.Mf(0.0) * wd**2 + S.Uf(0.0)
    p = np.sqrt(2 * S.Mf(th) * (E - S.Uf(th)))
    crit = np.abs(np.gradient(1 / p, th))
    lam = 2 * math.pi / p
    print("=== 2. Đồng hồ WKB: rô-to là biến 'nặng'? ===")
    print(f"tại θ̇_d = {wd:.2f}: p_θ dọc quỹ đạo {p.min():.2f}–{p.max():.2f} ħ; bước sóng de Broglie theo góc 2π/p = "
          f"{lam.min():.2f}–{lam.max():.2f} rad (thế biến thiên trên ~1 rad)")
    print(f"tiêu chí WKB |d(1/p)/dθ|: max {crit.max():.3f}, TB {crit.mean():.3f}  (cần ≪ 1)")
    print("PW với đồng hồ rô-to tác dụng J (Part III, N3): 1−F = 0.054 (J = 3ħ) → 0.014 (J = 60ħ), sàn do trộn vòng;"
          " giới hạn WdW 0.31 (J = 3ħ) → 0.05")
    print(f"⟹ động cơ J ≈ {p.mean():.1f} ħ: WKB {'đạt' if crit.max() < 0.1 else 'chỉ gần đúng'}; "
          "thời gian WKB nổi lên với sai số vài % trong một vòng nhưng không đếm vòng\n")


def bk():
    """Quãng đường tự do trung bình nucleon từ thế quang học KD 2003 (giai_so_part2_mo6): λ = ħv/(2W)."""
    from giai_so_part2_mo6 import kd_neutron
    R = 1.2 * 70**(1 / 3)
    print("=== 3. Brown–Kuchař (bụi không tương tác): nucleon trong hạt nhân có phải 'bụi'? ===")
    print(f"{'E (MeV)':>8} {'W_V (MeV)':>10} {'W_V+W_D (MeV)':>14} {'λ thể tích (fm)':>16} {'λ cả bề mặt (fm)':>17} {'λ/R':>10}")
    for E in (5.0, 10.0, 20.0, 40.0):
        p = kd_neutron(70, 31, E)
        v = math.sqrt(2 * E / MN)                       # v/c
        hv = 197.327 * v
        lv = hv / (2 * max(p["WV"], 1e-6)); lt = hv / (2 * (p["WV"] + p["WD"]))
        print(f"{E:8.0f} {p['WV']:10.2f} {p['WV']+p['WD']:14.2f} {lv:16.1f} {lt:17.1f} {lt/R:5.2f}–{lv/R:.1f}")
    print(f"(R = {R:.1f} fm.) Nucleon ở 5–40 MeV: λ ≈ vài fm tới hàng chục fm ~ R ⟹ hạt nhân 'bán trong suốt': nucleon "
          "chưa va chạm qua cả hạt nhân, nhưng KHÔNG đủ tự do để làm bụi Brown–Kuchař suốt nhiều T0;")
    print("đồng hồ Brown–Kuchař trong học thuyết chỉ là đồng hồ lý tưởng hình thức (Ĥ_C = p̂, Part III)\n")


def entropic():
    print("=== 4. Đồng hồ entropy ===")
    print("entropy hệ S(ρ_S) theo Part III (N1): 0 → 1.250 (1 T0) → 1.572 (5 T0) → 1.615 (20 T0) → 1.615 (60 T0) nat:"
          " BÃO HÒA sau ~20 T0 ⟹ không phải đồng hồ dài hạn")
    for rate, dE, lab in ((0.50, 2.05, "V.C.1 (Davies)"), (0.61, 0.41, "chặn gãy + bể rô-to (khớp Part I)")):
        dS = dE * 10 / 1.0
        print(f"   bể, {lab}: {rate:.2f} bản ghi/T0 × ΔE/k_BT = {dS:.1f} k_B ⟹ {rate*dS:.1f} k_B mỗi T0, tăng TUYẾN TÍNH mãi")
    print("⟹ entropy của BỂ (số bản ghi × ΔE/kT) là đồng hồ entropy đúng; entropy của hệ thì không\n")


def erker():
    print("=== 5. Đồng hồ tự hành (Erker 2017) + bất đẳng thức bất định nhiệt động (TUR) ===")
    for lab, cv, dS in (("chỉ bản ghi Carnot", 2.05, 4.1), ("mọi bản ghi", 1.52, 4.1)):
        N = 1 / cv**2
        print(f"   {lab}: độ chính xác N = 1/CV² = {N:.2f}; giới hạn TUR Σ/2 = {dS/2:.2f} ⟹ đạt {N/(dS/2):.0%} giới hạn")
    print("   yêu cầu σ_t ≤ 0.071 T0 trên T = 5 T0 với r tích tắc/T0, độ chính xác N: σ_t = √(T/(rN)) ⟹ rN ≥ T/σ²")
    for dt, T in ((0.071, 5), (0.27, 5)):
        rN = T / dt**2
        Smin = 2 * rN
        print(f"   σ_t = {dt} T0: rN ≥ {rN:.0f} ⟹ entropy ≥ 2rN = {Smin:.0f} k_B mỗi T0 ⟹ ở k_BT = 1 MeV tỏa ≥ {Smin:.0f} MeV/T0 "
              f"(so năng lượng kích thích ~20–30 MeV: thiếu {Smin/30:.0f}–{Smin/20:.0f} lần)")
    print()


if __name__ == "__main__":
    sw(); wkb(); bk(); entropic(); erker()
