"""Kim K̂ là đồng hồ NGOÀI hạt nhân: kiểm số đề xuất "K̂ = chuỗi xung attosecond laser/tia X".

Cùng hệ S (36 chiều) và Lindblad E_dual (Carnot + bơm) của giai_so_hai_vai_tro / giai_so_kim_chi_gio.
  A  Đơn vị vật lý: σ_t ≤ 0.071 T0 và nhịp τ = T0/25 đổi ra giây và MeV (băng thông ħ/2σ, photon h/τ);
     so với xung attosecond ngắn nhất đã đo (43 as) và chu kỳ chuỗi xung attosecond (~1.3 fs).
  B  Rung không tích lũy (mỗi lần đọc lệch ε ~ N(0, σ_K²), như độ dài/độ trễ xung):
     ρ̃(t̂) = E_ε[ρ(max(0, t̂ + ε))]. Đo sai số so với ρ(t̂), độ tương phản động lực còn lại,
     và phần dư Lindblad ‖dρ̃/dt − ℒρ̃‖ (rung không tích lũy chỉ làm nhòe điều kiện đầu?).
  C  Rung tích lũy (đồng hồ tự chạy, phương sai Dt): ρ̃(t) = E[e^{(t+W_t)ℒ}]ρ0 = e^{t(ℒ + Dℒ²/2)}ρ0
     (Gambini–Porto–Pullin). Kiểm đẳng thức bằng cầu phương, đo D tối đa cho sai số ≤ nhiễu thống kê 0.016.
Chạy: python3 giai_so_kim_ngoai.py   (≈ 1–2 phút)
"""
import math

import numpy as np

from giai_so_wdw_pw_e7 import T0, channels
from giai_so_kim_chi_gio import psi0, tdist

T0_S = 4.136e-22          # s (Bảng I: ħω0 = 10 MeV)
HBAR_MEV_S = 6.582119569e-22
STAT = 0.016              # nhiễu thống kê 1/√4000 của phép kiểm vai trò 1


def liouville(H, Ls):
    d = H.shape[0]; I = np.eye(d); S = sum(L.conj().T @ L for L in Ls)
    M = -1j * (np.kron(H, I) - np.kron(I, H.T)) - 0.5 * (np.kron(S, I) + np.kron(I, S.T))
    for L in Ls:
        M = M + np.kron(L, L.conj())
    w, V = np.linalg.eig(M)
    return M, w, V, np.linalg.inv(V)


def states(w, V, c, times, d):
    """ρ(t) = V e^{wt} V⁻¹ρ0 cho cả mảng thời gian (nt, d, d)."""
    return ((np.exp(np.outer(times, w)) * c) @ V.T).reshape(len(times), d, d)


def part_a():
    print("=== A. Yêu cầu của K̂ theo đơn vị vật lý ===")
    s99 = 0.071 * T0_S; tick = T0_S / 25
    rows = [("T0", T0_S), ("σ_t cho 99% xếp đúng ba pha = 0.071 T0", s99), ("nhịp τ = T0/25", tick)]
    for lab, v in rows:
        print(f"  {lab:42s} = {v:.3e} s = {v * 1e24:.1f} ys (yoctogiây)")
    bw = HBAR_MEV_S / (2 * s99)
    print(f"  băng thông tối thiểu của xung dài σ_t:  ħ/(2σ_t) = {bw:.1f} MeV")
    print(f"  photon có ≥ 1 chu kỳ trong nhịp τ:      h/τ = {2 * math.pi * HBAR_MEV_S / tick:.0f} MeV"
          f"  (ngưỡng sinh pion ≈ 140 MeV; ngân sách kích thích 20–30 MeV)")
    for lab, v in (("xung attosecond cô lập ngắn nhất đã đo (43 as)", 43e-18),
                   ("chu kỳ chuỗi xung attosecond (nửa chu kỳ laser 800 nm ≈ 1.33 fs)", 1.33e-15),
                   ("photon XUV 100 eV — năng lượng", None)):
        if v is None:
            print(f"  {lab:62s}: 1e−4 MeV = {100e-6 / bw:.1e} lần băng thông cần")
        else:
            print(f"  {lab:62s}: {v:.2e} s = {v / T0_S:7.0f} T0 = {v / s99:7.0f} × σ_t yêu cầu")
    return s99


def smeared(w, V, c, d, th, sig, xg, wg):
    """ρ̃(t̂) = E_ε[ρ(max(0, t̂ + ε))], ε ~ N(0, sig²): xung đến trước khi kích thì hệ còn ở ρ0."""
    if sig == 0:
        return states(w, V, c, th, d)
    tt = np.maximum(th[:, None] + sig * xg[None, :], 0.0)
    S = states(w, V, c, tt.ravel(), d).reshape(len(th), len(xg), d, d)
    return np.einsum('j,tjab->tab', wg, S)


def part_b(w, V, Vi, c, d):
    print("\n=== B. Rung không tích lũy σ_K (mỗi lần đọc kim lệch ε ~ N(0, σ_K²)) ===")
    th = np.linspace(0, 5 * T0, 126)                      # đọc kim mỗi 0.04 T0 trên 5 T0
    ideal = states(w, V, c, th, d)
    xg, wg = np.polynomial.hermite_e.hermegauss(120)      # cầu phương Gauss–Hermite (xác suất)
    wg = wg / wg.sum()
    span = lambda R: max(tdist(R[i], R[j]) for i in range(0, len(th), 5) for j in range(i + 5, len(th), 5))
    span0 = span(ideal)
    print("  σ_K/T0    σ_K (s)    max‖ρ̃−ρ‖   biên độ động lực max‖ρ̃(t̂₁)−ρ̃(t̂₂)‖   nửa nhóm ‖ρ̃(t̂+1T0) − e^{T0ℒ}ρ̃(t̂)‖, t̂ = 3σ+0.5T0")
    out = {}
    for sK in (0.0, 0.001, 0.002, 0.003, 0.005, 0.01, 0.071, 0.27, 1.0, 5.0, 104.0):
        sig = sK * T0
        sm = smeared(w, V, c, d, th, sig, xg, wg)
        err = max(tdist(sm[i], ideal[i]) for i in range(len(th)))
        t1 = 3 * sig + 0.5 * T0
        pair = smeared(w, V, c, d, np.array([t1, t1 + T0]), sig, xg, wg)
        prop = (V @ (np.exp(w * T0) * (Vi @ pair[0].reshape(-1)))).reshape(d, d)
        semi = tdist(pair[1], prop)
        out[sK] = err
        print(f"  {sK:8.3f}  {sK * T0_S:.2e}  {err:9.4f}   {span0:.3f} → {span(sm):.3f}"
              f"                          {semi:.1e}")
    return out


def part_c(w, V, c, d):
    print("\n=== C. Rung tích lũy (đồng hồ tự chạy, Var[t̂ − t] = D t): ℒ_eff = ℒ + (D/2)ℒ² ===")
    th = np.linspace(0, 5 * T0, 51)
    ideal = states(w, V, c, th, d)
    xg, wg = np.polynomial.hermite_e.hermegauss(80); wg = wg / wg.sum()
    print("  D/T0      σ_t tại 5 T0   ‖quadrature − e^{t(ℒ+Dℒ²/2)}‖   max‖ρ̃ − ρ‖ (sai số do đồng hồ)")
    for D in (1e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1):
        Dt = D * T0
        tt = np.maximum(th[:, None] + np.sqrt(Dt * th)[:, None] * xg[None, :], 0.0)
        quad = np.einsum('j,tjab->tab', wg, states(w, V, c, tt.ravel(), d).reshape(len(th), len(xg), d, d))
        gpp = ((np.exp(np.outer(th, w + 0.5 * Dt * w * w)) * c) @ V.T).reshape(len(th), d, d)
        diff = max(tdist(quad[i], gpp[i]) for i in range(len(th)))
        err = max(tdist(quad[i], ideal[i]) for i in range(len(th)))
        print(f"  {D:7.0e}   {math.sqrt(Dt * 5 * T0) / T0:8.3f} T0      {diff:.2e}                    {err:.4f}")


if __name__ == '__main__':
    ls, H, Ls = channels()
    d = H.shape[0]
    M, w, V, Vi = liouville(H, Ls)
    rho0 = np.outer(psi0(), psi0().conj())
    c = Vi @ rho0.reshape(-1)
    iss = int(np.argmin(np.abs(w)))
    rho_ss = (V[:, iss] * c[iss]).reshape(d, d)
    rho_ss = rho_ss / np.trace(rho_ss)
    LIOU = lambda r: (M @ r.reshape(-1)).reshape(d, d)
    print(f"S {d} chiều, Liouvillian {d * d}×{d * d}; khe phổ (tốc độ hồi phục chậm nhất) "
          f"= {sorted(-w.real)[1]:.4f} ω0 ⟹ thời gian hồi phục {1 / sorted(-w.real)[1] / T0:.2f} T0\n")
    part_a()
    part_b(w, V, Vi, c, d)
    part_c(w, V, c, d)
