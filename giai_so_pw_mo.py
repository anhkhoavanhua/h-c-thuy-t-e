"""Page–Wootters có tiêu tán cho E_dual: dρ_CS/dt = −i[H, ρ] + D_Carnot + D_pump, rồi ρ_S(θ0) = ⟨θ0|ρ_CS|θ0⟩/Tr.

Kiểm các khẳng định: (i) θ0 nhỏ ⟹ ρ_S thuần (Tr ρ_S² ≈ 1); (ii) θ0 "quét qua ngưỡng" ⟹ độ tinh khiết sụp;
(iii) θ0 lớn ⟹ trạng thái mới lệch pha vĩnh viễn.
Mô hình: H (rô-to + dao động tử + ghép −g x sinθ, dấu đúng), D_Carnot với L_c = √Γc û₋ (Davies T=0),
D_pump với L_p = √Γ0 e^{iθ} f(ℓ̂), f = √max(0, 1 − ℓ²/ℓ0²); Γ0 = 0.2, Γc = 0.5, ℓ0 = 4 (chế độ ⟨ℓ⟩ ≈ 2.4).
Tham chiếu PW lý tưởng (Moreva): trạng thái thuần một chiều (cặp ±ℓ) — ρ_S(θ0) thuần với mọi θ0.
Đơn vị ħ = m = ω0 = 1, I = 1, ΔE = 1, g = 0.1. Chạy: python3 giai_so_pw_mo.py
"""
import math

import numpy as np

from giai_so_u_hermite_drude import basis

NMAX = 6


def cond_state(rho, ls, th, dC, dS):
    c = np.exp(1j * ls * th)
    R = rho.reshape(dC, dS, dC, dS)
    rs = np.einsum('a,anbm,b->nm', c, R, c.conj())
    return rs / np.trace(rs).real, np.trace(rs).real


def lindblad(rho, H, Ls):
    d = -1j * (H @ rho - rho @ H)
    for L in Ls:
        LdL = L.conj().T @ L
        d += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return d


def profile(rho, ls, x, dC, dS, nth=64):
    ths = np.linspace(0, 2 * math.pi, nth, endpoint=False)
    pur, xm, pw = [], [], []
    for th in ths:
        rs, w = cond_state(rho, ls, th, dC, dS)
        pur.append(float(np.real(np.trace(rs @ rs))))
        xm.append(float(np.real(np.trace(x @ rs))))
        pw.append(w)
    pur, xm, pw = np.array(pur), np.array(xm), np.array(pw)
    A = 2 * np.sum(pw * xm * np.sin(ths)) / np.sum(pw)          # hệ số sinθ (xấp xỉ)
    return pur, xm, pw, A


def main():
    ls, cos_t, sin_t, pth, H, u = basis(lmax=8, nmax=NMAX)
    dC, dS = len(ls), NMAX + 1
    a = np.diag(np.sqrt(np.arange(1, dS)), 1)
    x = (a + a.T) / math.sqrt(2)
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(dS))
    E, V = np.linalg.eigh(H)
    Lc = math.sqrt(0.5) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    up = np.diag(np.ones(dC - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(dS))
    Ls = [Lc, Lp]

    # tham chiếu PW lý tưởng: cặp ±ℓ gần suy biến quanh E ≈ 5, chồng chập cho ⟨ℓ⟩ lớn nhất
    k = next(k for k in range(len(E) - 1) if 4.5 < E[k] < 6 and abs(E[k + 1] - E[k]) < 1e-3)
    S2 = V[:, k:k + 2]
    w2, u2 = np.linalg.eigh(S2.conj().T @ Lop @ S2)
    psi = S2 @ u2[:, -1]
    pur, xm, pw, A = profile(np.outer(psi, psi.conj()), ls, x, dC, dS)
    print("=== Tham chiếu PW đơn nhất (trạng thái thuần, đồng hồ quay một chiều, ⟨ℓ⟩ = "
          f"{float(np.real(psi.conj() @ Lop @ psi)):.2f}) ===")
    print(f"Tr ρ_S²(θ0): min {pur.min():.6f}, max {pur.max():.6f}  ⟹ thuần với MỌI θ0; hệ số sinθ của ⟨x⟩ ≈ {A:+.4f}\n")

    # Lindblad theo thời gian ngoài t, xuất phát |ℓ=0, n=0⟩
    rho = np.zeros((dC * dS, dC * dS), complex)
    i0 = int(np.where(ls == 0)[0][0]) * dS
    rho[i0, i0] = 1
    dt, t = 0.01, 0.0
    marks = [0, 2, 5, 10, 20, 40, 80, 120]
    print("=== PW có tiêu tán (Lindblad bơm + làm nguội), trạng thái điều kiện ρ_S(θ0) theo thời gian ngoài t ===")
    print(f"{'t':>6} {'⟨ℓ⟩':>7} {'Tr ρ_CS²':>9} {'Tr ρ_S²(θ0): min':>17} {'TB':>7} {'max':>7} "
          f"{'biến thiên theo θ0':>18} {'A sinθ':>8} {'p(θ0) max/min':>14}")
    snaps = {}
    for mt in marks:
        while t < mt - 1e-9:
            k1 = lindblad(rho, H, Ls); k2 = lindblad(rho + dt / 2 * k1, H, Ls)
            k3 = lindblad(rho + dt / 2 * k2, H, Ls); k4 = lindblad(rho + dt * k3, H, Ls)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
        pur, xm, pw, A = profile(rho, ls, x, dC, dS)
        snaps[mt] = (pur, xm)
        lm = float(np.real(np.trace(Lop @ rho)))
        print(f"{mt:6g} {lm:7.3f} {float(np.real(np.trace(rho @ rho))):9.4f} {pur.min():17.4f} {pur.mean():7.4f} "
              f"{pur.max():7.4f} {pur.max()-pur.min():18.4f} {A:+8.4f} {pw.max()/max(pw.min(),1e-300):14.2f}")
    d1 = np.max(np.abs(snaps[120][0] - snaps[80][0]))
    d2 = np.max(np.abs(snaps[120][1] - snaps[80][1]))
    print(f"\nthay đổi giữa t = 80 và t = 120: max|ΔTr ρ_S²| = {d1:.2e}, max|Δ⟨x⟩| = {d2:.2e} (đã tới trạng thái dừng)")
    print("ρ_S(θ0) tính trên θ0 ∈ [0, 2π); θ0 + 2π cho cùng ρ_S (tuần hoàn theo định nghĩa)")


if __name__ == "__main__":
    main()
