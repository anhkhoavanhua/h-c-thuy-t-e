"""Kiểm kiến trúc "đồng hồ hạt nhân gamma" Page–Wootters cho E_dual (HocThuyetE (7).pdf, A ~ 70, ħω0 = 10 MeV).

G1  GDR làm đồng hồ: E_G = 31.2 A^{-1/3} + 20.6 A^{-1/6} MeV, Γ ≈ 5 MeV (hệ thống thực nghiệm A ~ 60–80):
    thời gian sống ħ/Γ, số chu kỳ dao động trước khi tắt E_G/(2πΓ), hệ số phẩm chất.
G2  Chuỗi xung gamma bên trong hạt nhân: để đọc thời gian với σ_t ≤ 0.07 T0 (99% xếp đúng 3 pha) trên khoảng
    t ≈ một chu trình (lưu + đợt): (a) phát Poisson: k ≥ (t/σ)² photon; (b) chuỗi đều lý tưởng: τ ≤ √12 σ, mỗi xung
    cần băng thông ΔE ≥ ħ/(2σ); năng lượng tối thiểu so với E* của hạt nhân (≈ 20–30 MeV).
G3  Dải quay làm đồng hồ: E(J) = J(J+1) ħ²/2𝓘, ħ²/2𝓘 = E(2⁺)/6; spin J cần để chu kỳ quay ≤ T0/25.
G4  "Bắt trọn 100% va chạm Carnot": khoảng cách giữa các va chạm liên tiếp trong đợt (Langevin Part I) ⟹ tỉ lệ cặp va
    chạm được xếp đúng thứ tự với σ_t = 0.07 T0.
G5  Cầu nối ²²⁹Th: δν/ν = K·δE_vt/E_vt với K ≈ |ΔV|/E_c (ΔV = hiệu năng lượng Coulomb/hạt nhân giữa hai trạng thái,
    cỡ 0.1–0.5 MeV); độ nhạy nhỏ nhất với độ chính xác đồng hồ 10⁻¹⁸; hệ số triệt tiêu khi đồng hồ chậm lấy trung
    bình dao động nhanh của động cơ.
Chạy: python3 giai_so_dong_ho_gamma.py   (≈ 1 phút)
"""
import math
from multiprocessing import Pool

import numpy as np

HBAR = 6.582119569e-22          # MeV s
T0 = 4.1357e-22                 # s
HW = 10.0
E_TH = 8.355733e-6              # MeV


def g1(A=70):
    EG = 31.2 * A**(-1 / 3) + 20.6 * A**(-1 / 6)
    print("=== G1. Cộng hưởng lưỡng cực khổng lồ (GDR) làm đồng hồ ===")
    print(f"{'Γ (MeV)':>8} {'E_G (MeV)':>10} {'sống ħ/Γ (s)':>13} {'= T0':>7} {'chu kỳ dao động trước khi tắt':>30} {'Q = E_G/Γ':>10}")
    for G in (4.0, 5.0, 6.0):
        tau = HBAR / G
        cyc = EG / (2 * math.pi * G)
        print(f"{G:8.1f} {EG:10.2f} {tau:13.2e} {tau/T0:7.2f} {cyc:30.2f} {EG/G:10.2f}")
    print(f"cần ≥ 25 tích tắc/T0 trên ≥ 5 T0 ⟹ ≥ 125 tích tắc kết hợp; GDR cho < 1 chu kỳ ⟹ KHÔNG làm đồng hồ được\n")


def g2(sig=0.0708, t_cycle=4.0, Estar=(20.0, 30.0)):
    sig_s = sig * T0
    dE = HBAR / (2 * sig_s)
    k_poisson = (t_cycle / sig)**2
    tau = math.sqrt(12) * sig
    k_reg = t_cycle / tau
    print("=== G2. Chuỗi xung gamma phát từ chính hạt nhân ===")
    print(f"mục tiêu σ_t = {sig} T0 = {sig_s:.2e} s trên một chu trình t ≈ {t_cycle} T0")
    print(f"(a) phát Poisson (phân rã gamma là ngẫu nhiên): cần k ≥ (t/σ)² = {k_poisson:.0f} photon trong {t_cycle} T0")
    print(f"    mỗi photon ≥ 1 MeV ⟹ ≥ {k_poisson*1e-3:.1f} GeV; so E* = {Estar[0]:.0f}–{Estar[1]:.0f} MeV: "
          f"thiếu {k_poisson/Estar[1]:.0f}–{k_poisson/Estar[0]:.0f} lần")
    print(f"(b) chuỗi đều lý tưởng: khoảng cách τ ≤ √12 σ = {tau:.3f} T0 ⟹ {k_reg:.0f} xung; mỗi xung định vị trong σ cần "
          f"băng thông ΔE ≥ ħ/2σ = {dE:.1f} MeV")
    print(f"    năng lượng tối thiểu ≈ {k_reg*dE:.0f} MeV so E* = {Estar[0]:.0f}–{Estar[1]:.0f} MeV: thiếu "
          f"{k_reg*dE/Estar[1]:.1f}–{k_reg*dE/Estar[0]:.1f} lần; và phát đều đòi một nguồn kết hợp bên ngoài (hạt nhân "
          "phát gamma theo Poisson)\n")


def g3(E2s=(0.5, 0.8, 1.0, 1.2)):
    print("=== G3. Dải quay hạt nhân làm đồng hồ (tích tắc = một vòng quay, cần ≤ T0/25) ===")
    need = 2 * math.pi * HBAR / (T0 / 25)                 # ħω_rot cần (MeV)
    print(f"cần ħω_rot = dE/dJ ≥ {need:.0f} MeV")
    print(f"{'E(2+) (MeV)':>11} {'ħ²/2𝓘 (MeV)':>12} {'J cần':>8} {'chu kỳ quay ở J = 20 (T0)':>26}")
    for E2 in E2s:
        a = E2 / 6
        J = (need / a - 1) / 2
        Trot20 = 2 * math.pi * HBAR / ((2 * 20 + 1) * a) / T0
        print(f"{E2:11.2f} {a:12.3f} {J:8.0f} {Trot20:26.2f}")
    print("(hạt nhân A ~ 70 vỡ/phân hạch trước J ~ 40–60; dải quay còn là đồng hồ tuần hoàn — vấn đề đếm vòng)\n")


def _lang(seed):
    from giai_so_dong_ho_photon import _traj
    from giai_so_burst_g import cfg_for
    return _traj((cfg_for(0.1), 0.1, seed, 80 * 2 * math.pi, 0.05))[1]


def g4(sig=0.0708):
    from math import erf
    with Pool() as pool:
        evs = pool.map(_lang, range(22000, 22016))
    Tn = 2 * math.pi
    gaps = []
    for ev in evs:
        ks = [e[1] / Tn for e in ev if e[0] == "khóa" and e[1] > 5 * Tn]
        gaps += [b - a for a, b in zip(ks, ks[1:]) if b - a < 1.0]      # va chạm liên tiếp trong cùng đợt
    g = np.array(gaps)
    p = 0.5 * (1 + np.vectorize(erf)(g / (2 * sig)))                    # P(đúng thứ tự) = Φ(Δ/(√2σ))
    print("=== G4. Đồng hồ σ_t = 0.07 T0 có 'bắt trọn 100%' các va chạm Carnot không? ===")
    print(f"{len(g)} cặp va chạm liên tiếp trong đợt: khoảng cách trung vị {np.median(g):.3f} T0, 10% {np.percentile(g,10):.4f} T0, "
          f"tỉ lệ < 0.01 T0: {np.mean(g < 0.01):.1%}")
    print(f"tỉ lệ cặp xếp đúng thứ tự: TB {p.mean():.1%}; cặp phân giải tốt (> 2σ): {np.mean(g > 2*sig):.1%}")
    print("⟹ ở mức ĐỢT xếp đúng 99%, ở mức TỪNG va chạm thì không: các va chạm rung sát nhau dưới ngưỡng phân giải\n")


def g5(Kmin=1.2e4, Kmax=6e4, prec=1e-18, Tint=1.0):
    print("=== G5. Cầu nối với đồng hồ ²²⁹Th (phép đo tĩnh) ===")
    for dV in (0.1, 0.3, 0.5):
        K = dV / E_TH
        print(f"|ΔV| = {dV} MeV (hiệu Coulomb/hạt nhân giữa đồng phân và trạng thái cơ bản) ⟹ K = |ΔV|/E_c = {K:.1e}; "
              f"δE_vt/E_vt nhỏ nhất đo được (δν/ν = {prec:g}) = {prec/K:.1e}")
    f_eng = 1 / T0
    supp = 1 / (math.pi * f_eng * Tint)
    print(f"dao động của E_vt ở tần số động cơ f = 1/T0 = {f_eng:.2e} Hz bị lấy trung bình trong thời gian hỏi đồng hồ "
          f"T = {Tint} s: hệ số còn lại ≲ 1/(π f T) = {supp:.1e}")
    print("⟹ đồng hồ Th chỉ thấy (i) độ dịch HẰNG của E_vt — bị hấp thụ vào định nghĩa ν0, không quan sát được tuyệt đối; "
          "(ii) biến thiên CHẬM (≲ kHz) của E_vt theo thời gian/môi trường. Học thuyết cần dự đoán loại (ii) mới kiểm được\n")


if __name__ == "__main__":
    g1(); g2(); g3(); g4(); g5()
