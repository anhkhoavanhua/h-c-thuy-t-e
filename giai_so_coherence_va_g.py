"""Giải số hai đề xuất: (1) thời gian giữ pha τ_coherence (chặn Pauli, độ rộng mode tập thể);
(2) hằng số ghép cặp g = 0.1·g_max, g_max = √(k·ΔE).

(A) Số học SI: τ_Pauli, ħ/Γ_giant, g_max, g, độ dốc dλ/dθ|0, tách mode ω± tại cộng hưởng.
(B) Mô phỏng lại chu kỳ gãy/khóa với g = 0.1 (các kết quả 12/12 trước đó là với g = 0),
    so thời gian tái đồng bộ với các cửa sổ τ_relax, ħ/Γ (Γ = 1 và 3 MeV), τ_Pauli.

Chạy: python3 giai_so_coherence_va_g.py [A|B]   (B ≈ 5–6 phút)
"""
import math
import sys

from giai_so_bo_thong_so_hat_nhan import HBAR_C, HBAR_MEV_S, HW0, MEV, MPC2, C, part_b

EPS_F = 37.0     # MeV, năng lượng Fermi hạt nhân
T0 = 2 * math.pi * HBAR_MEV_S / HW0
TAU_RELAX = 5e-15 / (0.2 * C)
TAU_PAULI = (EPS_F / HW0)**2 * TAU_RELAX


def part_a():
    print("=== A. Số học SI (ħω0 = 10 MeV, ΔE = 10 MeV) ===")
    print(f"T0 = {T0:.4e} s;  τ_relax = {TAU_RELAX:.4e} s = {TAU_RELAX/T0:.3f} T0")
    print(f"τ_Pauli = (ε_F/(E−ε_F))² τ_relax = ({EPS_F}/{HW0})² = {(EPS_F/HW0)**2:.3f} × τ_relax "
          f"= {TAU_PAULI:.4e} s = {TAU_PAULI/T0:.3f} T0   (đề xuất: 1.2e-21 s ≈ 3 T0)")
    for G in (1.0, 3.0, 4.5):
        tg = HBAR_MEV_S / G
        print(f"ħ/Γ, Γ = {G} MeV: {tg:.4e} s = {tg/T0:.3f} T0   (đề xuất cho Γ=1÷3 MeV: (2÷6)e-21 s ≈ (5÷15) T0)")
    k = MPC2 * HW0**2 / HBAR_C**2
    gmax = math.sqrt(k * HW0)
    l0 = HBAR_C / math.sqrt(MPC2 * HW0)
    unit_force = HW0 / l0
    print(f"k = mω0² = {k:.4f} MeV/fm²;  g_max = √(kΔE) = {gmax:.4f} MeV/fm = {gmax*MEV/1e-15:.1f} N")
    print(f"đơn vị lực không thứ nguyên ħω0/ℓ0 = {unit_force:.4f} MeV/fm ⟹ g_max = {gmax/unit_force:.4f} đơn vị")
    g = 0.1 * gmax
    print(f"g = 0.1 g_max = {g:.4f} MeV/fm = {g*MEV/1e-15:.2f} N = {g/unit_force:.3f} đơn vị")
    I, m, kk, dE, ka, gd = 1.0, 1.0, 1.0, 1.0, 1 / 3, 0.1
    slope = (ka * (I * kk - m * dE) - gd * I + gd * m * ka**2) / (I + m * ka**2)
    print(f"dλ/dθ|0 (đủ số hạng) = ({ka*(I*kk-m*dE):.4f} − {gd*I:.4f} + {gd*m*ka**2:.4f}) / {I+m*ka**2:.4f} "
          f"= {slope:.4f}   (đề xuất ≈ −0.09, bỏ số hạng +gmκ²)")
    wp, wm = math.sqrt(1 + gd), math.sqrt(1 - gd)       # hệ tự do (1)–(2) tuyến tính hóa, ωθ = ω0 = 1
    print(f"Tách mode tự do tại cộng hưởng: ω± = √(1±g) = {wp:.5f}, {wm:.5f};  ħΔω = {(wp-wm)*HW0:.4f} MeV;"
          f"  chu kỳ phách = {2*math.pi/(wp-wm):.3f} = {1/(wp-wm):.2f} T0\n")


def part_b_g():
    windows = [("τ_relax", TAU_RELAX / (HBAR_MEV_S / HW0)),
               ("ħ/Γ3MeV", (HBAR_MEV_S / 3.0) / (HBAR_MEV_S / HW0)),
               ("ħ/Γ1MeV", (HBAR_MEV_S / 1.0) / (HBAR_MEV_S / HW0)),
               ("τ_Pauli", TAU_PAULI / (HBAR_MEV_S / HW0)),
               ("10T0", 10 * 2 * math.pi)]
    print("Cửa sổ (đơn vị 1/ω0): " + ", ".join(f"{n}={v:.3f} ({v/(2*math.pi):.2f} T0)" for n, v in windows))
    part_b(g=0.1, windows=windows)


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b_g}[name.upper()]()
