"""Giải số hai câu hỏi mở của Part II.

Câu 1 — toán tử vận tốc tương đối û = (κ/2I){cosθ, p_θ} − p_x/m:
  (A) kiểm dấu dạng khai triển: {cosθ,p_θ}/2 = cosθ p_θ + (iħ/2) sinθ (đề xuất ghi dấu −).
  (B) L = √Γ û (Hermite): năng lượng ⟨H⟩ theo thời gian (tăng = nung nóng).
  (C) L = √Γ û₋ (dạng Davies T=0: chỉ giữ thành phần hạ năng lượng của û trong cơ sở riêng của H):
      quỹ đạo lượng tử MCWF — số bước nhảy tới khi gần trạng thái cơ bản, năng lượng mỗi bước nhảy,
      thời gian chờ giữa các bước nhảy (CV), kiểm phân bố hình học/mũ.
Câu 2 — bể Drude–Lorentz, ħω_c = 4.5 MeV:
  (D) ma sát hiệu dụng tại ω0 và tại Ω = 2.63ω0; cân bằng CHÍNH XÁC của dao động tử trong bể Drude
      (chuyển động Brown lượng tử tuyến tính, định lý thăng giáng–tiêu tán) so với Caldeira–Leggett,
      Lindblad–Diósi và Lindblad nhiệt a/a†.
Đơn vị ħ = m = ω0 = 1 (năng lượng theo ħω0 = 10 MeV), I = 1, κ = 1/3, ΔE = 1, g = 0.1, γ = 0.7.
Chạy: python3 giai_so_u_hermite_drude.py [A|B|C|D]   (C ≈ 1–2 phút)
"""
import math
import random
import sys

import numpy as np

KAP, I_, DE, G, GAM = 1 / 3, 1.0, 1.0, 0.1, 0.7


def basis(lmax=8, nmax=6):
    ls = np.arange(-lmax, lmax + 1)
    d = len(ls)
    up = np.diag(np.ones(d - 1), -1).astype(complex)          # e^{iθ}
    cos_t, sin_t = (up + up.conj().T) / 2, (up - up.conj().T) / (2j)
    pth = np.diag(ls).astype(complex)
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
    x = (a + a.conj().T) / math.sqrt(2)
    px = 1j * (a.conj().T - a) / math.sqrt(2)
    Ir, Ip = np.eye(d), np.eye(nmax + 1)
    K = lambda A, B: np.kron(A, B)
    H = (K(pth @ pth / (2 * I_) - DE * cos_t, Ip) + K(Ir, a.conj().T @ a + 0.5 * Ip) - G * K(sin_t, x))
    sym = 0.5 * (cos_t @ pth + pth @ cos_t)
    u = KAP / I_ * K(sym, Ip) - K(Ir, px)
    return ls, cos_t, sin_t, pth, H, u


def part_a():
    ls, cos_t, sin_t, pth, H, u = basis()
    sym = 0.5 * (cos_t @ pth + pth @ cos_t)
    plus = cos_t @ pth + 0.5j * sin_t
    minus = cos_t @ pth - 0.5j * sin_t
    inner = slice(1, -1)   # bỏ hàng/cột biên của cơ sở cắt
    print("=== A. Dạng khai triển của ½{cosθ, p_θ} ===")
    print(f"‖½{{cos,p}} − (cos p + (i/2) sin)‖ = {np.linalg.norm((sym - plus)[inner, inner]):.2e}")
    print(f"‖½{{cos,p}} − (cos p − (i/2) sin)‖ (đề xuất) = {np.linalg.norm((sym - minus)[inner, inner]):.2e}")
    print(f"‖û − û†‖ = {np.linalg.norm(u - u.conj().T):.2e};  ‖(cos p − (i/2)sin) − (…)†‖ = "
          f"{np.linalg.norm((minus - minus.conj().T)[inner, inner]):.2e} (không Hermite)\n")


def lind_evolve(rho, H, Ls, T, dt):
    def rhs(r):
        d = -1j * (H @ r - r @ H)
        for L in Ls:
            LdL = L.conj().T @ L
            d += L @ r @ L.conj().T - 0.5 * (LdL @ r + r @ LdL)
        return d
    out, t = [], 0.0
    while t < T:
        k1 = rhs(rho); k2 = rhs(rho + dt / 2 * k1); k3 = rhs(rho + dt / 2 * k2); k4 = rhs(rho + dt * k3)
        rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        out.append((t, float(np.real(np.trace(H @ rho)))))
    return rho, out


def start_state(ls, nmax=6, l0=3):
    v = np.zeros(len(ls) * (nmax + 1), complex)
    v[int(np.where(ls == l0)[0][0]) * (nmax + 1)] = 1
    return v


def part_b():
    ls, _, _, _, H, u = basis()
    psi = start_state(ls)
    rho = np.outer(psi, psi.conj())
    _, tr = lind_evolve(rho, H, [math.sqrt(0.2) * u], 10.0, 0.002)
    E = [e for _, e in tr]
    print("=== B. L = √Γ û (Hermite), Γ = 0.2, xuất phát |ℓ=3, n=0⟩ ===")
    print(f"⟨H⟩: t=0 → {E[0]:.3f}; t=2.5 → {E[1249]:.3f}; t=5 → {E[2499]:.3f}; t=10 → {E[-1]:.3f}; "
          f"tăng đơn điệu: {'CÓ' if all(E[i+1] >= E[i]-1e-9 for i in range(len(E)-1)) else 'KHÔNG'}")
    print(f"(trần của cơ sở cắt: năng lượng TB của trạng thái hỗn hợp tối đa = {np.trace(H).real/H.shape[0]:.3f})\n")


def part_c(ntraj=300, Gam=0.5, dt=0.01, Tmax=300.0):
    ls, _, _, _, H, u = basis()
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    Ldown_e = np.triu(ue, 1)                                   # chỉ giữ ⟨j|û|i⟩ với E_j < E_i
    L = math.sqrt(Gam) * (V @ Ldown_e @ V.conj().T)
    Heff = H - 0.5j * (L.conj().T @ L)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    psi0 = start_state(ls)
    E0 = float(np.real(psi0.conj() @ H @ psi0))
    Eg = E[0]
    rng = random.Random(7)
    njumps, dEs, waits_all, first_wait = [], [], [], []
    for _ in range(ntraj):
        psi, t, r, tj, n = psi0.copy(), 0.0, rng.random(), 0.0, 0
        while t < Tmax:
            psi = U @ psi
            t += dt
            if np.vdot(psi, psi).real < r:
                Eb = float(np.real(np.vdot(psi, H @ psi) / np.vdot(psi, psi)))
                psi = L @ psi
                psi /= np.linalg.norm(psi)
                Ea = float(np.real(np.vdot(psi, H @ psi)))
                dEs.append(Eb - Ea)
                waits_all.append(t - tj)
                if n == 0:
                    first_wait.append(t - tj)
                tj, n, r = t, n + 1, rng.random()
                if Ea - Eg < 0.05:
                    break
        njumps.append(n)
    nj = np.array(njumps)
    dE = np.array(dEs)
    wt = np.array(waits_all)
    print(f"=== C. MCWF với L = √Γ û₋ (Davies T=0), Γ = {Gam}, {ntraj} quỹ đạo, xuất phát |ℓ=3, n=0⟩ ===")
    print(f"E đầu = {E0:.3f}, E cơ bản = {Eg:.3f} ⟹ cần tiêu tán {E0-Eg:.3f} ħω0")
    print(f"số bước nhảy tới trạng thái cơ bản: TB {nj.mean():.2f}, trung vị {np.median(nj):.0f}, min {nj.min()}, max {nj.max()}; "
          f"CV {nj.std()/nj.mean():.3f} (hình học có CV = √(1−p) ≈ {math.sqrt(1-1/nj.mean()):.3f})")
    counts = {k: int(np.sum(nj == k)) for k in range(int(nj.min()), int(nj.max()) + 1)}
    print("   phân bố số bước:", ", ".join(f"{k}:{v}" for k, v in counts.items()))
    print(f"năng lượng mỗi bước nhảy: TB {dE.mean():.3f} ħω0, trung vị {np.median(dE):.3f}, 10–90%: "
          f"{np.percentile(dE,10):.3f}–{np.percentile(dE,90):.3f}; P(ΔE < 0.5) = {np.mean(dE < 0.5):.1%}")
    print(f"thời gian chờ giữa các bước nhảy: TB {wt.mean():.3f}, CV {wt.std()/wt.mean():.3f} (mũ: CV = 1)\n")


def part_d():
    wc, kT, w0, gam = 0.45, 0.1, 1.0, GAM
    print("=== D. Bể Drude–Lorentz, ħω_c = 4.5 MeV (0.45 ħω0), k_BT = 1 MeV, γ = 0.7 ===")
    for Om in (1.0, 2.63):
        print(f"ma sát hiệu dụng Re γ̂(ω) = γω_c²/(ω_c²+ω²) tại ω={Om}: {gam*wc*wc/(wc*wc+Om*Om):.4f} "
              f"(= {wc*wc/(wc*wc+Om*Om):.1%} của γ)")
    wgrid = np.linspace(1e-6, 60.0, 600001)

    def moments(wc_):
        gh = gam * wc_ / (wc_ - 1j * wgrid)
        chi = 1 / (w0**2 - wgrid**2 - 1j * wgrid * gh)
        cth = 1 / np.tanh(wgrid / (2 * kT))
        x2 = np.trapezoid(cth * chi.imag, wgrid) / math.pi
        p2 = np.trapezoid(wgrid**2 * cth * chi.imag, wgrid) / math.pi
        return x2, p2

    for label, wc_ in (("Drude ω_c = 0.45 (4.5 MeV)", 0.45), ("Drude ω_c = 5 (50 MeV)", 5.0), ("Drude ω_c = 20", 20.0)):
        x2, p2 = moments(wc_)
        print(f"{label}: ⟨x²⟩ = {x2:.4f}, ⟨p²⟩ = {p2:.4f}, ⟨x²⟩⟨p²⟩ = {x2*p2:.4f}")
    print("so sánh (từ các kiểm chứng trước, cùng γ = 0.7, k_BT = 1 MeV): Caldeira–Leggett ⟨x²⟩ = 0.100 (ρ âm); "
          "Lindblad–Diósi 0.782; Lindblad nhiệt a/a† = ½coth(1/2k_BT) = "
          f"{0.5/math.tanh(1/(2*kT)):.4f}; không ghép = 0.5000\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "D", "C"]:
        {"A": part_a, "B": part_b, "C": part_c, "D": part_d}[name.upper()]()
