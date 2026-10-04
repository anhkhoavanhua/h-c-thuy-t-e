"""Giải số đề xuất Part II: toán tử nhảy chung L = √Γ e^{−iθ} a, bộ tiêu tán Caldeira–Leggett, ranh giới k_BT ≲ 1 MeV.

(A) Chiều quay ở 1 MeV: tỉ lệ thời gian dΦ/dt > 0, tốc độ ròng từng quỹ đạo (đề xuất: "θ̇>0 với xác suất > 91.6%,
    tốc độ trôi ròng luôn dương").
(B) Năng lượng và mô-men mỗi bước nhảy của L = e^{−iθ}a: |ℓ,1⟩ → |ℓ−1,0⟩ so với va chạm Carnot cổ điển
    (ΔE ≈ 0.40 ħω0, Δℓ = Λκ|cosθ| ≤ κ√(2m_eff ΔE)).
(C) Γ_eff = ⟨ΔE⟩/ħ · r_burst · ⟨N⟩: thứ nguyên; tốc độ va chạm cổ điển cần tái lập.
(D) Caldeira–Leggett ở k_BT = 1 và 5 MeV: tính dương của ρ, ⟨x²⟩⟨p²⟩ ≥ 1/4; so với dạng Lindblad L = Ax + iBp.
Chạy: python3 giai_so_lindblad_joint.py [A|B|C|D]   (A ≈ 1 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

from giai_so_burst_g import cfg_for
from giai_so_winding_lindblad import _traj, ops

T0 = 2 * math.pi


def part_a(N=8):
    cfg = cfg_for(0.1)
    T, warm, every = 70 * T0, 10 * T0, 0.05
    with Pool() as pool:
        res = pool.map(_traj, [(cfg, 0.1, 33000 + i, T, warm, every) for i in range(N)])
    pos, tot, rates = 0, 0, []
    for th, _ in res:
        phi = th[int(warm / every):]
        d = np.diff(phi)
        pos += int(np.sum(d > 0)); tot += len(d)
        rates.append((phi[-1] - phi[0]) / (2 * math.pi) / ((T - warm) / T0))
    print("=== A. Chiều quay ở 1 MeV (cùng hạt giống với kiểm chứng số cuộn) ===")
    print(f"tỉ lệ thời gian dΦ/dt > 0: {pos/tot:.1%};  tốc độ ròng từng quỹ đạo (vòng/T0): "
          + ", ".join(f"{r:+.2f}" for r in rates) + f";  số quỹ đạo âm: {sum(1 for r in rates if r < 0)}/{N}\n")


def part_b():
    kap, dE_cl = 1 / 3, 0.40
    print("=== B. Bước nhảy L = e^{−iθ}a: |ℓ,1⟩ → |ℓ−1,0⟩ (I = 1, ħ = ω0 = 1) ===")
    print(f"{'ℓ':>3} {'ΔE_hạt':>7} {'ΔE_rô-to = ((ℓ−1)²−ℓ²)/2I':>27} {'ΔE tổng':>8} {'so với Carnot 0.40':>19}")
    for l in (-2, -1, 0, 1, 2, 3, 4, 5):
        drot = ((l - 1)**2 - l**2) / 2
        tot = -1 + drot
        print(f"{l:3d} {-1:7d} {drot:27.2f} {tot:8.2f} {abs(tot)/dE_cl:18.1f}×" + ("   (rô-to NHẬN năng lượng)" if drot > 0 else ""))
    print(f"Cổ điển: Δℓ mỗi va chạm = Λκ|cosθ| ≤ κ√(2 m_eff ΔE) ≤ {kap*math.sqrt(2*dE_cl):.2f} ħ "
          f"(với ΔE = {dE_cl}, m_eff ≤ 1), so với 1 ħ của e^(−iθ)\n")


def part_c():
    rb, Nb = 0.152, 5.12
    rate = rb * Nb / T0
    print("=== C. Γ_eff ===")
    print("Thứ nguyên: [⟨ΔE⟩/ħ] = 1/s, [r_burst] = 1/s, [⟨N⟩] = 1 ⟹ [Γ_eff] = 1/s² (không phải 1/s)")
    print(f"Tốc độ va chạm cổ điển cần tái lập: r_burst·⟨N⟩ = {rb}×{Nb}/T0 = {rate:.4f} /đv t "
          f"(= {rate*T0:.3f} /T0); với ⟨n⟩ = 0.071, Γ = tốc độ/⟨a†a⟩ = {rate/0.071:.2f} /đv t")
    print(f"Công suất nếu mỗi nhảy mất 3.5 ħω0 (ℓ=3): {rate*3.5:.3f} ħω0/đv t, so với Carnot cổ điển {rate*0.40:.3f} "
          f"(gấp {3.5/0.40:.1f} lần)\n")


def part_d(nmax=40):
    a, ad = ops(nmax)
    x = (a + ad) / math.sqrt(2)
    p = 1j * (ad - a) / math.sqrt(2)
    H = ad @ a + 0.5 * np.eye(nmax + 1)
    gam = 0.7
    print("=== D. Caldeira–Leggett (ma sát −(iγ/2)[x,{p,ρ}], khuếch tán γ m k_BT — đúng hệ số đề xuất), γ = 0.7 ===")

    def run(rhs, T=15.0, dt=0.002):
        rho = np.zeros((nmax + 1, nmax + 1), complex); rho[0, 0] = 1
        t, mineig = 0.0, 1.0
        while t < T:
            k1 = rhs(rho); k2 = rhs(rho + dt / 2 * k1); k3 = rhs(rho + dt / 2 * k2); k4 = rhs(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            mineig = min(mineig, float(np.min(np.linalg.eigvalsh((rho + rho.conj().T) / 2))))
        return float(np.real(np.trace(x @ x @ rho))), float(np.real(np.trace(p @ p @ rho))), mineig

    for kT_mev in (1.0, 5.0):
        kT = kT_mev / 10
        q = 0.5 / math.tanh(1 / (2 * kT))

        def cl(rho):
            d = -1j * (H @ rho - rho @ H)
            d += -0.5j * gam * (x @ (p @ rho + rho @ p) - (p @ rho + rho @ p) @ x)
            xx = x @ rho - rho @ x
            return d - gam * kT * (x @ xx - xx @ x)

        A, B = math.sqrt(2 * gam * kT), math.sqrt(gam / (8 * kT))
        Lk = A * x + 1j * B * p

        def lind(rho):
            d = -1j * (H @ rho - rho @ H)
            LdL = Lk.conj().T @ Lk
            return d + Lk @ rho @ Lk.conj().T - 0.5 * (LdL @ rho + rho @ LdL)

        for label, f in (("Caldeira–Leggett", cl), ("Lindblad L = Ax + iBp (dạng Diósi)", lind)):
            x2, p2, me = run(f)
            print(f"k_BT = {kT_mev:g} MeV: {label}: ⟨x²⟩ = {x2:.4f} (lượng tử cân bằng {q:.4f}), ⟨p²⟩ = {p2:.4f}, "
                  f"⟨x²⟩⟨p²⟩ = {x2*p2:.4f} (≥ 0.25), trị riêng nhỏ nhất của ρ = {me:.2e}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["B", "C", "D", "A"]:
        {"A": part_a, "B": part_b, "C": part_c, "D": part_d}[name.upper()]()
