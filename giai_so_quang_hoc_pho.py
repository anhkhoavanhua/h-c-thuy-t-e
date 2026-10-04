"""Giải số hai khẳng định:
(A) Tần số cắt ħω_c ≈ 80–100 MeV "đến từ phần ảo thế quang học": dùng dạng thể tích Koning–Delaroche (2003)
    W_V(ΔE) = w1 ΔE²/(ΔE² + w2²), w1 = 12.195 + 0.0167A MeV, w2 = 73.55 + 0.0795A MeV, ΔE = E − ε_F.
    Tính W_V ở năng lượng của mô hình, tốc độ tắt tương đương γ = 2W/ħω0, và ΔE cần để γ = γx = 0.7.
(B) "Không tồn tại chuyển dời 0.4 ħω0 giữa các trạng thái riêng của H": khe phổ trong cửa sổ năng lượng của
    rô-to quay (ℓ ≈ 2.4), và phân bố năng lượng chuyển dời của û₋ có trọng số |⟨j|û|i⟩|².
Chạy: python3 giai_so_quang_hoc_pho.py [A|B]
"""
import math
import sys

import numpy as np

from giai_so_u_hermite_drude import basis

HW0 = 10.0


def part_a():
    print("=== A. Phần ảo thế quang học (Koning–Delaroche, số hạng thể tích) ===")
    for A in (56, 208):
        w1, w2 = 12.195 + 0.0167 * A, 73.55 + 0.0795 * A
        print(f"A={A}: w1 = {w1:.2f} MeV, w2 = {w2:.2f} MeV (thang bão hòa)")
        for dE in (5, 10, 20, 26.3, 50, 100):
            W = w1 * dE**2 / (dE**2 + w2**2)
            print(f"   ΔE = {dE:5.1f} MeV: W_V = {W:6.3f} MeV ⟹ γ = 2W/ħω0 = {2*W/HW0:.3f} (mô hình dùng γx = 0.70)")
        Wt = 0.7 * HW0 / 2
        x = math.sqrt(Wt * w2**2 / (w1 - Wt))
        print(f"   γ = 0.70 (W = {Wt} MeV) cần ΔE = {x:.0f} MeV trên mặt Fermi")
    print("(Năng lượng hạt trong mô hình: ~0.57 ħω0 ≈ 6 MeV; xung Carnot ~4 MeV; rô-to tới hạn ~51 MeV nhưng rô-to không phải nucleon)\n")


def part_b():
    ls, _, _, _, H, u = basis()
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    print("=== B. Phổ H và chuyển dời û₋ ===")
    win = (E > 2.0) & (E < 6.0)
    Ew = E[win]
    gaps = np.diff(Ew)
    print(f"{win.sum()} mức riêng trong cửa sổ E ∈ (2, 6) ħω0 (vùng rô-to quay ℓ≈2–3); khe liền kề: min {gaps.min():.3f}, "
          f"trung vị {np.median(gaps):.3f}, max {gaps.max():.3f}; số khe < 0.5: {int(np.sum(gaps < 0.5))}/{len(gaps)}")
    idx = np.where(win)[0]
    dEs, wts = [], []
    for i in idx:
        for j in range(i):
            w = abs(ue[j, i])**2
            if w > 1e-12:
                dEs.append(E[i] - E[j]); wts.append(w)
    dEs, wts = np.array(dEs), np.array(wts)
    wts = wts / wts.sum()
    mean = float(np.sum(dEs * wts))
    print(f"chuyển dời đi xuống từ cửa sổ này, trọng số |⟨j|û|i⟩|²: TB {mean:.3f} ħω0; "
          f"trọng số dành cho ΔE < 0.5: {float(wts[dEs < 0.5].sum()):.1%}; ΔE ∈ [0.3, 0.5]: {float(wts[(dEs>=0.3)&(dEs<0.5)].sum()):.1%}; "
          f"ΔE > 1.5: {float(wts[dEs > 1.5].sum()):.1%}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
