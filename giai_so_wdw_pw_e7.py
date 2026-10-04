"""Trạng thái Wheeler–DeWitt / Page–Wootters trọn vẹn cho học thuyết E_dual theo HocThuyetE (7).pdf.

Vũ trụ = Đồng hồ C ⊗ Hệ S ⊗ Môi trường E:
  C: đồng hồ chuỗi xung (V.C.4), mỗi tích tắc τ = T0/25 (25 xung mỗi chu kỳ T0);
  S: E_dual lượng tử (rô-to ℓ ≤ 4 ⊗ hạt n ≤ 3), H Bảng I với −g x sinθ (giai_so_u_hermite_drude.basis);
  E: mỗi tích tắc một ancilla 3 mức {không, nhận lượng tử Carnot, nhận lượng tử bơm} — mô hình va chạm đơn nhất
     U = exp(−i[τ H ⊗ 1 + √τ Σ_j (L_j ⊗ |j⟩⟨0| + h.c.)]),  L_c = √0.5 û₋ (Davies), L_p = √0.2 e^{iθ} f(ℓ̂) (V.C.1).
Trạng thái lịch sử: |Ψ⟩⟩ = N^{-1/2} Σ_k |k⟩_C ⊗ U_k ⋯ U_1 |ψ0⟩_S |0…0⟩_E.
Kiểm:
 (1) Ràng buộc WdW rời rạc (Giovannetti–Lloyd–Maccone): Ĵ|Ψ⟩⟩ = |Ψ⟩⟩, Ĵ = Σ_k |k+1⟩⟨k| ⊗ U_{k+1} (đóng vòng).
     Người quan sát ngoài: mọi phân bố trên C là đều, vũ trụ tĩnh.
 (2) Người quan sát trong: ρ_S(k) = Tr_E ⟨k|Ψ⟩⟩ so với nghiệm Lindblad liên tục tại t = kτ (sai số mô hình va chạm O(τ)).
 (3) Mũi tên thời gian và bản ghi: S toàn cục = 0; S(ρ_S(k)) và số bản ghi Carnot tăng theo k; điều kiện hóa trên
     một bản ghi E cụ thể ⟹ S thuần.
Chạy: python3 giai_so_wdw_pw_e7.py   (≈ 10 s)
"""
import math

import numpy as np

from giai_so_u_hermite_drude import basis

T0 = 2 * math.pi
LMAX, NMAX, NT = 4, 3, 7          # NT tích tắc ⟹ E có 3^NT trạng thái
TAU = T0 / 25                      # chuỗi xung V.C.4: 25 xung mỗi chu kỳ


def channels():
    ls, cos_t, sin_t, pth, H, u = basis(lmax=LMAX, nmax=NMAX)
    E, V = np.linalg.eigh(H)
    Lc = math.sqrt(0.5) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    up = np.diag(np.ones(len(ls) - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(NMAX + 1))
    return ls, H, [Lc, Lp]


def collision_unitary(H, Ls):
    d = H.shape[0]
    G = np.kron(H, np.eye(3)) * TAU
    for j, L in enumerate(Ls, start=1):
        P = np.zeros((3, 3)); P[j, 0] = 1
        G = G + math.sqrt(TAU) * (np.kron(L, P) + np.kron(L.conj().T, P.T))
    w, Q = np.linalg.eigh(G)
    return (Q * np.exp(-1j * w)) @ Q.conj().T, d


def apply_tick(psi, U, d, k):
    """Áp U lên trục S (0) và ancilla thứ k (trục k+1) của tensor psi."""
    psi = np.moveaxis(psi, k + 1, 1)
    sh = psi.shape
    psi = (U @ psi.reshape(d * 3, -1)).reshape(sh)
    return np.moveaxis(psi, 1, k + 1)


def lindblad_ref(H, Ls, rho, T, dt=0.002):
    S = sum(L.conj().T @ L for L in Ls)
    f = lambda r: -1j * (H @ r - r @ H) - 0.5 * (S @ r + r @ S) + sum(L @ r @ L.conj().T for L in Ls)
    n = int(round(T / dt))
    for _ in range(n):
        k1 = f(rho); k2 = f(rho + dt / 2 * k1); k3 = f(rho + dt / 2 * k2); k4 = f(rho + dt * k3)
        rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return rho


def vn(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


def main():
    ls, H, Ls = channels()
    U, d = collision_unitary(H, Ls)
    print(f"=== Vũ trụ C ⊗ S ⊗ E: S {d} chiều (ℓ ≤ {LMAX}, n ≤ {NMAX}), {NT} tích tắc τ = T0/25, E 3^{NT} = {3**NT} chiều ===")
    print(f"‖U†U − 1‖ = {np.linalg.norm(U.conj().T @ U - np.eye(3 * d)):.1e} (đơn nhất)")
    # trạng thái đầu: bó sóng thỏa ràng buộc (pha Quá khứ, con trỏ của bể Carnot), ℓ0 = 3, θ0 = 0
    cr = np.exp(-(ls - 3)**2 / (4 * 0.7**2)); cr = cr / np.linalg.norm(cr)
    al = 1j * (1 / 3) * 3 / math.sqrt(2)
    n = np.arange(NMAX + 1)
    cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in n]) * math.exp(-abs(al)**2 / 2)
    cx /= np.linalg.norm(cx)
    psi0 = np.kron(cr, cx)
    psi = np.zeros((d,) + (3,) * NT, complex)
    psi[(slice(None),) + (0,) * NT] = psi0
    hist = [psi.copy()]
    for k in range(NT):
        psi = apply_tick(psi, U, d, k)
        hist.append(psi.copy())
    N = len(hist)                                   # các số đọc đồng hồ k = 0..NT
    Psi = np.array(hist) / math.sqrt(N)             # |Ψ⟩⟩ = Σ_k |k⟩ ⊗ ψ(k)

    # (1) ràng buộc WdW rời rạc: Ĵ|Ψ⟩⟩ với đóng vòng |0⟩ ← W |NT⟩, W = (U_NT ⋯ U_1)^{-1}
    JPsi = np.zeros_like(Psi)
    for k in range(N - 1):
        JPsi[k + 1] = apply_tick(Psi[k], U, d, k)
    back = Psi[N - 1].copy()
    for k in reversed(range(NT)):                   # áp U^† theo thứ tự ngược
        back = apply_tick(back, U.conj().T, d, k)
    JPsi[0] = back
    print("\n(1) Wheeler–DeWitt / Page–Wootters")
    print(f"‖Ĵ|Ψ⟩⟩ − |Ψ⟩⟩‖ = {np.linalg.norm(JPsi - Psi):.1e}   (ràng buộc thỏa chính xác ⟹ vũ trụ tĩnh)")
    pC = [float(np.sum(np.abs(Psi[k])**2)) for k in range(N)]
    print(f"phân bố số đọc đồng hồ P(k): min {min(pC):.6f}, max {max(pC):.6f} (= 1/{N} đều; người quan sát ngoài "
          "không thấy gì thay đổi)")
    print(f"S(toàn vũ trụ) = 0 (trạng thái thuần), ‖Ψ‖ = {np.linalg.norm(Psi):.12f}")

    # (2) + (3) người quan sát trong
    print("\n(2) Người quan sát trong: trạng thái điều kiện ρ_S(k) = Tr_E ⟨k|Ψ⟩⟩, so với Lindblad liên tục t = kτ")
    print(f"{'k':>3} {'t/T0':>6} {'‖ρ_S − ρ_Lindblad‖₁/2':>21} {'Tr ρ_S²':>8} {'S(ρ_S) = S(E)':>14} "
          f"{'P(≥1 bản ghi Carnot)':>21} {'P(≥1 bản ghi bơm)':>18}")
    rho_ref = np.outer(psi0, psi0.conj())
    for k in range(N):
        M = hist[k].reshape(d, -1)
        rS = M @ M.conj().T
        if k > 0:
            rho_ref = lindblad_ref(H, Ls, rho_ref, TAU)
        td = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(rS - rho_ref)))
        prob = np.abs(hist[k])**2
        tot = prob.sum(axis=0)                      # phân bố bản ghi E (3^NT)
        idx = np.indices(tot.shape).reshape(NT, -1).T
        tv = tot.reshape(-1)
        pcool = float(tv[np.any(idx == 1, axis=1)].sum())
        ppump = float(tv[np.any(idx == 2, axis=1)].sum())
        print(f"{k:3d} {k*TAU/T0:6.3f} {td:21.2e} {np.real(np.trace(rS @ rS)):8.4f} {vn(rS):14.4f} {pcool:21.4f} {ppump:18.4f}")

    # điều kiện hóa trên một bản ghi E cụ thể
    M = hist[-1].reshape(d, -1)
    pr = np.sum(np.abs(M)**2, axis=0)
    order = np.argsort(pr)[::-1]
    print(f"\n(3) Điều kiện hóa trên bản ghi môi trường cụ thể tại k = {NT}:")
    recs = np.array(np.unravel_index(np.arange(pr.size), (3,) * NT)).T
    rC = int(np.argmax(np.where(np.any(recs == 1, axis=1), pr, 0)))      # bản ghi có Carnot hay gặp nhất
    for r in dict.fromkeys(list(order[:3]) + [rC]):
        rec = np.unravel_index(r, (3,) * NT)
        v = M[:, r] / math.sqrt(pr[r])
        print(f"   bản ghi {''.join('·CB'[x] for x in rec)} (P = {pr[r]:.4f}): Tr ρ_S² = "
              f"{abs(np.vdot(v, v))**2:.6f} (thuần), ⟨H⟩ = {np.real(v.conj() @ H @ v):.3f} ħω0")
    print("   (· = không, C = nhận lượng tử Carnot, B = nhận lượng tử bơm)")


if __name__ == "__main__":
    main()
