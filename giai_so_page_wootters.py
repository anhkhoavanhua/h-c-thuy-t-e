"""Page–Wootters cho hệ (θ, x): chạy script đề xuất và kiểm tra.

(A) Script đề xuất (nguyên văn, chỉ thay scipy.linalg.eigh bằng numpy.linalg.eigh vì máy không có scipy).
(B) Kiểm tra:
    B1 dấu: ma trận sinθ của script so với (e^{iθ} − e^{−iθ})/2i; phép chiếu ⟨θ0| của script so với ⟨θ|ℓ⟩ = e^{iℓθ}.
    B2 trị riêng "gần 0 nhất" là trạng thái nào (thứ tự trong phổ, ⟨ℓ²⟩, phân bố p(θ) của đồng hồ).
    B3 chọn các nghiệm HΨ = EΨ ứng với đồng hồ QUAY (⟨ℓ⟩ ≈ 2–3): biên độ, pha của ⟨x⟩(θ0); so với κ sinθ (κ = 1/3),
       g sinθ/k (tĩnh) và g sinθ/(k − Ω²) (cưỡng bức cổ điển với Ω = ⟨ℓ⟩/I).
Đơn vị ħ = m = ω0 = 1, I = 1, ΔE = 1, g = 0.1.
Chạy: python3 giai_so_page_wootters.py
"""
import numpy as np


def solve_page_wootters():
    # ---- (A) nguyên văn đề xuất ----
    I = 1.0
    Delta_E = 1.0
    g = 0.1
    kappa = 1.0 / 3.0
    l_max = 8
    n_max = 6
    dim_C = 2 * l_max + 1
    dim_S = n_max + 1
    dim_total = dim_C * dim_S
    print(f"Kích thước Hilbert: C={dim_C}, S={dim_S} => Tổng={dim_total}")
    l_vals = np.arange(-l_max, l_max + 1)
    p_theta = np.diag(l_vals)
    H_rot = (p_theta ** 2) / (2.0 * I)
    cos_theta = 0.5 * (np.diag(np.ones(dim_C - 1), 1) + np.diag(np.ones(dim_C - 1), -1))
    sin_theta = -0.5j * (np.diag(np.ones(dim_C - 1), 1) - np.diag(np.ones(dim_C - 1), -1))
    H_C = H_rot - Delta_E * cos_theta
    H_S = np.diag(np.arange(dim_S) + 0.5)
    a = np.diag(np.sqrt(np.arange(1, dim_S)), 1)
    a_dag = a.T
    x_op = (a + a_dag) / np.sqrt(2.0)
    H_total = np.kron(H_C, np.eye(dim_S)) + np.kron(np.eye(dim_C), H_S) - g * np.kron(sin_theta, x_op)
    evals, evecs = np.linalg.eigh(H_total)
    idx_zero = np.argmin(np.abs(evals))
    E_zero = evals[idx_zero]
    psi_total = evecs[:, idx_zero]
    print(f"Trị riêng gần 0 nhất: E = {E_zero:.6e}")
    theta_steps = 100
    theta_arr = np.linspace(0, 4 * np.pi, theta_steps)
    x_mean = np.zeros(theta_steps)
    E_s_mean = np.zeros(theta_steps)
    norm_s = np.zeros(theta_steps)
    psi_mat = psi_total.reshape((dim_C, dim_S))
    for i, th in enumerate(theta_arr):
        clock_state = np.exp(1j * l_vals * th) / np.sqrt(dim_C)
        psi_S = np.dot(clock_state.conj(), psi_mat)
        norm = np.linalg.norm(psi_S)
        norm_s[i] = norm
        if norm > 1e-12:
            psi_S_norm = psi_S / norm
            x_mean[i] = np.real(np.vdot(psi_S_norm, x_op @ psi_S_norm))
            E_s_mean[i] = np.real(np.vdot(psi_S_norm, H_S @ psi_S_norm))
    print("\n--- KẾT QUẢ ĐO THEO PHA ĐỒNG HỒ ---")
    print(f"Biên độ dao động <x> đo được: [{np.min(x_mean):.4f}, {np.max(x_mean):.4f}]")
    print(f"Kỳ vọng lý thuyết liên kết cổ điển (kappa = 1/3): [-0.3333, +0.3333]")
    print(f"Năng lượng hạt trung bình: {np.mean(E_s_mean):.4f} hbar*omega_0")
    return dict(evals=evals, evecs=evecs, idx=idx_zero, H=H_total, sin=sin_theta, cos=cos_theta,
                l=l_vals, x=x_op, HS=H_S, dimC=dim_C, dimS=dim_S, theta=theta_arr, xm=x_mean, norm=norm_s, g=g)


def conditional(psi_mat, l_vals, th, x_op):
    """Trạng thái điều kiện đúng quy ước ⟨θ|ℓ⟩ = e^{iℓθ}: ψ_S(θ) = Σ_ℓ e^{iℓθ} Ψ_ℓ."""
    v = np.exp(1j * l_vals * th) @ psi_mat
    p = np.vdot(v, v).real
    return p, float(np.real(np.vdot(v, x_op @ v) / p))


def fit_sin(th, xm, w):
    """Khớp ⟨x⟩(θ) ≈ A sinθ + B cosθ + C (bình phương tối thiểu có trọng số p(θ))."""
    M = np.stack([np.sin(th), np.cos(th), np.ones_like(th)], 1)
    W = np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(M * W, xm * np.sqrt(w), rcond=None)
    return c


def checks(r):
    l_vals, dC, dS = r['l'], r['dimC'], r['dimS']
    up = np.diag(np.ones(dC - 1), -1)             # e^{iθ}: |ℓ⟩ → |ℓ+1⟩
    sin_true = (up - up.T) / (2j)
    print("\n=== B1. Kiểm dấu ===")
    print(f"‖sinθ(script) + sinθ(đúng)‖ = {np.linalg.norm(r['sin'] + sin_true):.2e}  ⟹ script dùng −sinθ")
    print("Phép chiếu của script: Σ e^{−iℓθ}Ψ_ℓ = ψ(−θ). Hai lỗi dấu triệt tiêu trong tương quan ⟨x⟩(θ) "
          "(vì đổi θ → −θ cũng đổi dấu sinθ). Mô phỏng dưới đây dùng dấu đúng.")
    up2 = np.kron(up, np.eye(dS)); Lsq = np.kron(np.diag(l_vals.astype(float)**2), np.eye(dS))
    Lop = np.kron(np.diag(l_vals.astype(float)), np.eye(dS))
    H = np.kron(np.diag(l_vals**2 / 2.0) - (up + up.T) / 2, np.eye(dS)) + np.kron(np.eye(dC), r['HS']) \
        - r['g'] * np.kron(sin_true, r['x'])
    E, V = np.linalg.eigh(H)
    k0 = int(np.argmin(np.abs(E)))
    print("\n=== B2. Trị riêng 'gần 0 nhất' ===")
    print(f"E = {E[k0]:.4f}, là mức thứ {k0} trong phổ (0 = cơ bản); mức cơ bản E0 = {E[0]:.4f}")
    th = np.linspace(0, 2 * np.pi, 400, endpoint=False)
    print(f"{'nghiệm':>10} {'E':>8} {'⟨ℓ⟩':>7} {'√⟨ℓ²⟩':>7} {'p(θ) max/min':>13} {'A (sin)':>9} {'B (cos)':>9} "
          f"{'κ':>6} {'g/k':>6} {'g/(k−Ω²)':>9}")
    # ứng viên: (i) nghiệm "gần 0"; (ii) các nghiệm riêng đơn lẻ trong E ∈ (3, 6);
    # (iii) chồng chập trong cặp gần suy biến (|ΔE| < 1e-3) chọn ⟨ℓ⟩ lớn nhất ⟹ đồng hồ quay MỘT chiều
    states = [("gần 0", E[k0], V[:, k0])]
    win = [k for k in range(len(E)) if 3.0 < E[k] < 6.0]
    for k in win[:3]:
        states.append(("riêng đơn", E[k], V[:, k]))
    for k in win:
        if k + 1 < len(E) and abs(E[k + 1] - E[k]) < 1e-3 and E[k + 1] < 6.0:
            S2 = V[:, k:k + 2]
            w2, u2 = np.linalg.eigh(S2.conj().T @ Lop @ S2)
            states.append(("cặp ±ℓ", 0.5 * (E[k] + E[k + 1]), S2 @ u2[:, -1]))
    print(f"số nghiệm riêng trong E∈(3,6): {len(win)}; ⟨ℓ⟩ của chúng: max |⟨ℓ⟩| = "
          f"{max(abs(float(np.real(V[:, k].conj() @ Lop @ V[:, k]))) for k in win):.2e}")
    for tag, Ek, vec in states:
        pm = vec.reshape(dC, dS)
        ps, xs = zip(*[conditional(pm, l_vals, t, r['x']) for t in th])
        ps, xs = np.array(ps), np.array(xs)
        A, B, C0 = fit_sin(th, xs, ps)
        lm = float(np.real(vec.conj() @ Lop @ vec))
        l2 = float(np.real(vec.conj() @ Lsq @ vec))
        Om = max(abs(lm), np.sqrt(l2), 1e-9)
        print(f"{tag:>10} {Ek:8.3f} {lm:7.3f} {np.sqrt(l2):7.3f} {ps.max()/max(ps.min(),1e-300):13.2e} "
              f"{A:+9.4f} {B:+9.4f} {1/3:6.3f} {r['g']:6.3f} {r['g']/(1-Om**2):+9.4f}")
    print("(A, B: hệ số khớp ⟨x⟩(θ) ≈ A sinθ + B cosθ + C theo trọng số p(θ); p(θ) max/min đo độ đều của đồng hồ)")
    print("Lưu ý: nghiệm riêng của H với ⟨ℓ⟩ ≈ 0 là chồng chập ℓ và −ℓ (đồng hồ quay cả hai chiều).")


if __name__ == "__main__":
    res = solve_page_wootters()
    checks(res)
