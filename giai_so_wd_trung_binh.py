"""Trung bình hấp thụ mặt ⟨W_D⟩ trên hàm sóng nucleon, và γ_x hiệu dụng.

Thừa số dạng mặt (đạo hàm Woods–Saxon, đỉnh = 1 tại r = R):  g(r) = 4 e^{(r−R)/a} / (1 + e^{(r−R)/a})²
Thừa số thể tích:  f(r) = 1 / (1 + e^{(r−R_V)/a_V})
Hàm sóng: dao động tử 3D, ℓ0 = 2.04 fm, R_nl(r) ∝ r^l e^{−r²/2ℓ0²} L_n^{l+1/2}(r²/ℓ0²).
Hình học: (i) R = 1.2 A^{1/3}, a = 0.54 fm (đề xuất); (ii) Koning–Delaroche neutron:
  r_D = 1.3424 − 0.01585 A^{1/3}, a_D = 0.5446 − 1.656e-4 A; r_V = 1.3039 − 0.4054 A^{−1/3}, a_V = 0.6778 − 1.487e-4 A.
γ_eff = 2(W_V⟨f⟩ + W_D⟨g⟩)/ħω0 tại ΔE = 7 MeV (W_V = 0.105, W_D = 3.448 MeV cho A = 56), so với đỉnh 0.711.
Chạy: python3 giai_so_wd_trung_binh.py
"""
import math

import numpy as np

A, L0, HW0 = 56, 2.04, 10.0
WV, WD = 0.105, 3.448          # MeV tại ΔE = 7 MeV (A = 56), từ giai_so_phan_bien.py


def laguerre(n, alpha, x):
    if n == 0:
        return np.ones_like(x)
    lm, l = np.ones_like(x), 1 + alpha - x
    for k in range(1, n):
        lm, l = l, ((2 * k + 1 + alpha - x) * l - (k + alpha) * lm) / (k + 1)
    return l


def density(n, l, r):
    R = r**l * np.exp(-r**2 / (2 * L0**2)) * laguerre(n, l + 0.5, r**2 / L0**2)
    p = r**2 * R**2
    return p / np.trapezoid(p, r)


def main():
    r = np.linspace(1e-4, 15.0, 30001)
    A13 = A ** (1 / 3)
    geoms = {
        "R = 1.2A^1/3, a = 0.54": (1.2 * A13, 0.54, 1.2 * A13, 0.54),
        "Koning–Delaroche": ((1.3424 - 0.01585 * A13) * A13, 0.5446 - 1.656e-4 * A,
                             (1.3039 - 0.4054 / A13) * A13, 0.6778 - 1.487e-4 * A),
    }
    orbits = [("0s", 0, 0, 0), ("0p", 0, 1, 1), ("1s", 1, 0, 2), ("0d", 0, 2, 2),
              ("1p", 1, 1, 3), ("0f", 0, 3, 3), ("2s", 2, 0, 4), ("1d", 1, 2, 4), ("0g", 0, 4, 4)]
    peak = 2 * (WV + WD) / HW0
    print(f"A = {A}, ℓ0 = {L0} fm; đỉnh γ_total (ΔE = 7 MeV) = {peak:.3f}")
    for gname, (RD, aD, RV, aV) in geoms.items():
        g = 4 * np.exp((r - RD) / aD) / (1 + np.exp((r - RD) / aD))**2
        f = 1 / (1 + np.exp((r - RV) / aV))
        print(f"\n== {gname}: R_D = {RD:.3f} fm, a_D = {aD:.3f}; R_V = {RV:.3f}, a_V = {aV:.3f} ==")
        print(f"{'quỹ đạo':>7} {'N':>2} {'r_rms (fm)':>10} {'⟨g⟩':>7} {'⟨f⟩':>7} {'γ_eff':>7} {'so với đỉnh':>12}")
        for name, n, l, N in orbits:
            p = density(n, l, r)
            rr = math.sqrt(np.trapezoid(r**2 * p, r))
            mg, mf = np.trapezoid(g * p, r), np.trapezoid(f * p, r)
            ge = 2 * (WV * mf + WD * mg) / HW0
            tag = "  ← gần mức Fermi (A=56)" if N in (3, 4) else ""
            print(f"{name:>7} {N:2d} {rr:10.3f} {mg:7.3f} {mf:7.3f} {ge:7.3f} {(ge-peak)/peak:+11.1%}{tag}")


if __name__ == "__main__":
    main()
