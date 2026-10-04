"""Kiểm số đề xuất "thời gian tự sinh từ góc rô-to mở cuộn" (Dirac–Kuchař/WKB) cho E_dual.

  A  [T_geom, p_θ] = iħ với T_geom = (I/⟨p_θ⟩) Φ̂?  Trên S¹ (ℓ nguyên): chéo của [Φ, L] trong cơ sở ℓ luôn bằng 0,
     nên min_Φ ‖[Φ, L] − i·1‖_F = √d (không có toán tử góc nào). Trên ℝ: [T_geom, p] = iħ·I/⟨p⟩, không phải iħ.
  B  Mở cuộn Φ ∈ ℝ có xóa sàn 1 − F ≈ 0.05 không?  Mở rộng N3 của giai_so_part3.clock_fidelity:
     rô-to trên vòng chu vi 2π/δ (mômen là bội của δ; δ = 1 là S¹ của học thuyết, δ → 0 là đường thẳng = "mở cuộn").
     Độ trùng Page–Wootters của hạt điều kiện theo góc so với Schrödinger, theo δ và độ rộng cửa sổ năng lượng
     (tw: 1 vòng … WdW tw = 50 T_rot).
  C  "Biên độ x vượt bán kính giam r_br ≈ 6.11 fm khi gãy": so x trên đa tạp với r_br.
Chạy: python3 giai_so_dong_ho_tu_sinh.py   (≈ 3–6 phút)
"""
import math

import numpy as np


def part_a():
    print("=== A. Toán tử thời gian góc: [Φ, L] = i có nghiệm không? ===")
    for lmax in (3, 10, 30):
        ls = np.arange(-lmax, lmax + 1).astype(float); d = ls.size
        # [Φ, L]_{ab} = Φ_ab (ℓ_b − ℓ_a): bình phương tối thiểu từng phần tử (chéo không điều khiển được)
        diff = ls[None, :] - ls[:, None]
        Phi = np.where(diff != 0, 0.0, 0.0).astype(complex)          # tối ưu: Φ_ab = 0 khử mọi phần tử ngoài chéo
        C = Phi * diff
        res = np.linalg.norm(C - 1j * np.eye(d))
        print(f"  ℓ ∈ [−{lmax}, {lmax}] (d = {d}): min‖[Φ, L] − i·1‖_F = {res:.4f} = √d = {math.sqrt(d):.4f}"
              f"  (sai số tương đối {res / math.sqrt(d):.0%}); Tr[Φ, L] = 0 ≠ Tr(i·1) = {d}i")
    for ell in (2.4, 3.0, 3.2):
        print(f"  trên ℝ: [T_geom, p_θ] = iħ·I/⟨p_θ⟩ = iħ/{ell} (I = 1, ⟨p_θ⟩ = {ell}ħ) — không phải iħ;"
              f" cặp liên hợp đúng: [T_geom, H_rot] = iħ·p̂/⟨p⟩ ≈ iħ chỉ gần ⟨p⟩")
    print("  Pauli: H_total bị chặn dưới ⟹ không có T tự liên hợp với [T, H] = iħ chính xác\n")


def clock_fid_delta(J, delta, tws, g=1.0, Om=0.3, nmax=10, half=16, sl=None):
    """Như giai_so_part3.clock_fidelity nhưng mômen rô-to ℓ ∈ J + δ·ℤ (vòng chu vi 2π/δ).
    e^{iθ} dời ℓ đi 1 = 1/δ bước lưới; chu kỳ của đồng hồ là T_rot/δ."""
    I = J / Om
    Trot = 2 * math.pi / Om
    sl = sl or max(1.0, math.sqrt(J / 3))
    step = int(round(1 / delta))
    ls = J + delta * np.arange(-half * step, half * step + 1)
    n = ls.size
    sh = np.zeros((n, n)); sh[np.arange(step, n), np.arange(n - step)] = 1.0      # |ℓ⟩ → |ℓ+1⟩
    sin_t = (sh - sh.T) / 2j
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    x = (a + a.T) / math.sqrt(2)
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    H = np.kron(np.diag(ls**2 / (2 * I)), np.eye(nmax + 1)) + np.kron(np.eye(n), Hx) - g * np.kron(sin_t, x)
    E, Uv = np.linalg.eigh(H)
    c = np.exp(-(ls - J)**2 / (4 * sl * sl)).astype(complex); c /= np.linalg.norm(c)
    chi0 = np.zeros(nmax + 1, complex); chi0[0] = 1
    amp = Uv.conj().T @ np.kron(c, chi0)
    E0 = float(E[np.argmax(np.abs(amp))])
    refs = {}
    for f in (0.1, 0.2, 0.3, 0.4):                          # Schrödinger của hạt với θ = Ωt
        t = f * Trot; psi = chi0.copy(); n_sub = 1500; dt = t / n_sub
        for i in range(n_sub):
            w_, v_ = np.linalg.eigh(Hx - g * math.sin(Om * (i + 0.5) * dt) * x)
            psi = v_ @ (np.exp(-1j * w_ * dt) * (v_.conj().T @ psi))
        refs[f] = psi
    out = {}
    for tw in tws:
        Psi = (Uv @ (np.exp(-((E - E0) * tw * Trot)**2 / 2) * amp)).reshape(n, nmax + 1)
        Fs = []
        for f, psi in refs.items():
            chi = np.exp(1j * ls * Om * f * Trot) @ Psi; chi /= np.linalg.norm(chi)
            Fs.append(abs(np.vdot(psi, chi))**2)
        out[tw] = 1 - min(Fs)
    return out, n * (nmax + 1)


def part_b():
    print("=== B. 'Mở cuộn' Φ ∈ ℝ: sàn 1 − F theo chu vi đồng hồ 2π/δ và cửa sổ năng lượng tw ===")
    tws = (0.35, 5.0, 50.0)
    print("  kiểm lại với N3 (δ = 1, nmax 10, ±16): ", end="")
    for J in (3, 30):
        r, dim = clock_fid_delta(J, 1.0, tws)
        print(f"J={J}: 1−F(1 vòng) {r[0.35]:.3f}, (WdW) {r[50.0]:.3f};  ", end="")
    print("(N3 gốc nmax 15, ±24: J=3 0.054/0.307, J=30 0.015/0.048)")
    print("   J   δ      chu kỳ đồng hồ   chiều    1−F (tw=0.35 T_rot)   (tw=5)   (tw=50, WdW)")
    for J in (3, 30):
        for delta in (1.0, 0.5, 0.25, 0.125):
            r, dim = clock_fid_delta(J, delta, tws)
            print(f"  {J:3d}  {delta:5.3f}   {1 / delta:5.0f} T_rot     {dim:5d}    {r[0.35]:.4f}               "
                  f"{r[5.0]:.4f}   {r[50.0]:.4f}", flush=True)
    print()


def part_c():
    print("=== C. Biên độ x trên đa tạp khi gãy so với r_br ===")
    l0 = 2.037                                  # fm (Bảng I)
    kappa = l0 / 3
    print(f"  r_br = 3ℓ0 = {3 * l0:.2f} fm;  |x| = κ|sinθ| ≤ κ = ℓ0/3 = {kappa:.3f} fm  (= {kappa / (3 * l0):.3f} r_br)")
    print(f"  dao động điểm không: x_rms = ℓ0/√2 = {l0 / math.sqrt(2):.2f} fm;  gãy là ngưỡng LỰC |λ| ≥ R = k·r_br,"
          f" không phải ngưỡng vị trí\n")


if __name__ == '__main__':
    part_a()
    part_c()
    part_b()
