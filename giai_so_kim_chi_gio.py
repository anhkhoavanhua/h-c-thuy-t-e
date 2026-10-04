"""Kiểm chứng số chứng minh "kim chỉ giờ" Page–Wootters của E_dual (HocThuyetE (7).pdf, Part III).
Dùng lại vũ trụ C ⊗ S ⊗ E của giai_so_wdw_pw_e7.py: S = rô-to ℓ ≤ 4 ⊗ hạt n ≤ 3 (20 chiều), mỗi tích tắc một ancilla
3 mức, U = exp(−i[τH ⊗ 1 + √τ Σ_j (L_j ⊗ |j⟩⟨0| + h.c.)]), τ = T0/25.

Bước 1  Ràng buộc Ĵ|Ψ⟩⟩ = |Ψ⟩⟩: đóng vòng |0⟩ ← W|N−1⟩ chọn W = (U_N⋯U_1)^{-1} ⟹ thỏa theo cách dựng (đồng nhất thức).
Bước 3a Khai triển Kraus: K0 = 1 − iτH − (τ/2)ΣL†L + O(τ^p), K_j = −i√τ L_j + O(τ^q): đo p, q.
Bước 3b Hội tụ Λ^k → e^{tℒ}: hai cách lấy giới hạn (τ → 0 với k cố định, và với t = kτ cố định).
Bước 4  Mũi tên: S(E) toàn bộ = S(ρ_S) (trạng thái SE thuần) — có đơn điệu không? Phản ví dụ làm lạnh thuần.
        So với entropy bản ghi thô Σ_i S(E_i) (từng ancilla) và entropy sinh Spohn.
Bước 5  Kim = bộ đếm Carnot N_C thay cho k: điều kiện hóa ⟨N_C = n|, so với Lindblad tại t̄(n) và độ nhòe σ_t(n).
Chạy: python3 giai_so_kim_chi_gio.py
"""
import math
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")             # ma trận 36×36: BLAS đa luồng chậm hơn 4–14 lần trên máy này

import numpy as np

from giai_so_wdw_pw_e7 import TAU, T0, LMAX, NMAX, channels, vn
import giai_so_wdw_pw_e7 as W


def unitary(H, Ls, tau):
    d = H.shape[0]; na = len(Ls) + 1
    G = np.kron(H, np.eye(na)) * tau
    for j, L in enumerate(Ls, start=1):
        P = np.zeros((na, na)); P[j, 0] = 1
        G = G + math.sqrt(tau) * (np.kron(L, P) + np.kron(L.conj().T, P.T))
    w, Q = np.linalg.eigh(G)
    return (Q * np.exp(-1j * w)) @ Q.conj().T


def kraus(H, Ls, tau):
    d = H.shape[0]; na = len(Ls) + 1
    U4 = unitary(H, Ls, tau).reshape(d, na, d, na)
    return [U4[:, j, :, 0] for j in range(na)]


_EIG = {}


def lind_at(H, Ls, rho, times):
    """Nghiệm Lindblad chính xác e^{tℒ}ρ tại các thời điểm: chéo hóa Liouvillian (d² × d²) một lần rồi lưu lại."""
    d = H.shape[0]
    key = (id(H),) + tuple(id(L) for L in Ls)
    if key not in _EIG:
        I = np.eye(d); S = sum(L.conj().T @ L for L in Ls)
        M = -1j * (np.kron(H, I) - np.kron(I, H.T)) - 0.5 * (np.kron(S, I) + np.kron(I, S.T))
        for L in Ls:
            M = M + np.kron(L, L.conj())
        w, V = np.linalg.eig(M)
        _EIG[key] = (w, V, np.linalg.inv(V))
    w, V, Vi = _EIG[key]
    c = Vi @ rho.reshape(-1)
    return [(V @ (np.exp(w * t) * c)).reshape(d, d) for t in times]


def lind_deriv(H, Ls, r):
    S = sum(L.conj().T @ L for L in Ls)
    return -1j * (H @ r - r @ H) - 0.5 * (S @ r + r @ S) + sum(L @ r @ L.conj().T for L in Ls)


def psi0():
    ls = np.arange(-LMAX, LMAX + 1)
    cr = np.exp(-(ls - 3)**2 / (4 * 0.7**2)); cr /= np.linalg.norm(cr)
    al = 1j * (1 / 3) * 3 / math.sqrt(2)
    cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in range(NMAX + 1)]); cx /= np.linalg.norm(cx)
    return np.kron(cr, cx)


def apply(Ks, r):
    r = sum(K @ r @ K.conj().T for K in Ks)
    r[np.abs(r) < 1e-200] = 0                    # tránh số subnormal (chậm ~100 lần)
    return r


def tdist(a, b):
    return 0.5 * np.sum(np.abs(np.linalg.eigvalsh((a - b + (a - b).conj().T) / 2)))


def step1(H, Ls):
    print("=== Bước 1. Ràng buộc WdW rời rạc ===")
    print("chạy lại giai_so_wdw_pw_e7 (7 tích tắc, E 3^7 chiều):")
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        W.main()
    print("   " + [l for l in buf.getvalue().splitlines() if "Ĵ|Ψ" in l][0])
    # ràng buộc có tính phân biệt không? thử một "lịch sử" SAI (dùng U với τ khác) nhưng vẫn đóng vòng bằng W nghịch đảo
    print("   kiểm độ phân biệt: thay U bằng U' (τ' = 2τ, hoặc H' = 0) nhưng vẫn đóng vòng bằng W' = (U'_N⋯U'_1)^{-1}:")
    for lab, Hx, tx in (("τ' = 2τ", H, 2 * TAU), ("H' = 0", 0 * H, TAU)):
        U = unitary(Hx, Ls, tx); d = H.shape[0]; NT = 5
        psi = np.zeros((d,) + (3,) * NT, complex); psi[(slice(None),) + (0,) * NT] = psi0()
        hist = [psi.copy()]
        for k in range(NT):
            psi = W.apply_tick(psi, U, d, k); hist.append(psi.copy())
        Psi = np.array(hist) / math.sqrt(NT + 1)
        J = np.zeros_like(Psi)
        for k in range(NT):
            J[k + 1] = W.apply_tick(Psi[k], U, d, k)
        b = Psi[-1].copy()
        for k in reversed(range(NT)):
            b = W.apply_tick(b, U.conj().T, d, k)
        J[0] = b
        print(f"      {lab}: ‖Ĵ'|Ψ'⟩⟩ − |Ψ'⟩⟩‖ = {np.linalg.norm(J - Psi):.1e}")
    print("   ⟹ mọi dãy unita đóng vòng bằng nghịch đảo của chính nó đều thỏa tới sai số máy: ràng buộc là ĐỒNG NHẤT THỨC "
          "của cách dựng, không chọn ra động lực E_dual (động lực nằm ở việc chọn U, không ở Ĵ = 1)\n")


def step3(H, Ls, rho0):
    print("=== Bước 3a. Khai triển Kraus theo τ ===")
    S = sum(L.conj().T @ L for L in Ls); d = H.shape[0]; I = np.eye(d)
    taus = [TAU / 2**m for m in range(5)]
    e0, ej = [], []
    for t in taus:
        K = kraus(H, Ls, t)
        e0.append(np.linalg.norm(K[0] - (I - 1j * t * H - 0.5 * t * S), 2))
        ej.append(max(np.linalg.norm(K[j] - (-1j) * math.sqrt(t) * Ls[j - 1], 2) for j in (1, 2)))
    p = np.polyfit(np.log(taus), np.log(e0), 1)[0]; q = np.polyfit(np.log(taus), np.log(ej), 1)[0]
    for t, a, b in zip(taus, e0, ej):
        print(f"   τ = T0/{T0/t:5.0f}: ‖K0 − (1 − iτH − τ/2 ΣL†L)‖ = {a:.2e}, ‖K_j + i√τ L_j‖ = {b:.2e}")
    print(f"   bậc đo được: K0 sai O(τ^{p:.2f}), K_j sai O(τ^{q:.2f})  (lý thuyết 2 và 3/2; pha −i của K_j không ảnh hưởng ρ)")
    Ks = kraus(H, Ls, TAU)
    print(f"   Σ K†K − 1 = {np.linalg.norm(sum(K.conj().T @ K for K in Ks) - I):.1e} (bảo toàn vết chính xác)\n")

    print("=== Bước 3b. Hội tụ Λ^k → e^{tℒ} ===")
    print("   (i) t = kτ cố định = 1 T0 và 5 T0, τ → 0:")
    for tf in (1.0, 5.0):
        errs = []
        for m in range(5):
            t = TAU / 2**m; k = int(round(tf * T0 / t)); Ks = kraus(H, Ls, t)
            r = rho0.copy()
            for _ in range(k):
                r = apply(Ks, r)
            errs.append(tdist(r, lind_at(H, Ls, rho0, [tf * T0])[0]))
        rat = [errs[i] / errs[i + 1] for i in range(4)]
        print(f"      t = {tf:g} T0: sai số vết " + ", ".join(f"{e:.2e}" for e in errs)
              + "  tỉ số khi chia đôi τ: " + ", ".join(f"{x:.2f}" for x in rat))
    print("   (ii) k cố định = 1, τ → 0 (cách viết trong chứng minh): ‖Λ(ρ) − ρ‖/τ vs ‖ℒρ‖:")
    Lr = lind_deriv(H, Ls, rho0)
    for m in range(0, 5, 2):
        t = TAU / 2**m
        diff = (apply(kraus(H, Ls, t), rho0) - rho0) / t
        print(f"      τ = T0/{T0/t:4.0f}: ‖(Λ−1)ρ/τ − ℒρ‖₁ = {np.sum(np.abs(np.linalg.eigvalsh((diff - Lr + (diff-Lr).conj().T)/2))):.2e}")
    print("   ⟹ đạo hàm khớp ℒ ở bậc O(τ); nhưng phương trình chủ chỉ nổi lên khi k → ∞ cùng τ → 0 với t = kτ cố định "
          "(k cố định thì t → 0, chỉ được đạo hàm tại 0)\n")


def marginal_entropy(Ks, r):
    """Entropy của ancilla vừa va chạm: ρ_E[i, j] = Tr(K_i ρ K_j†)."""
    na = len(Ks)
    rE = np.array([[np.trace(Ks[i] @ r @ Ks[j].conj().T) for j in range(na)] for i in range(na)])
    return vn(rE)


def step4(H, Ls, rho0):
    print("=== Bước 4. Mũi tên thời gian: S(E) có đơn điệu theo k? ===")
    Ks = kraus(H, Ls, TAU); r = rho0.copy()
    nk = 1500
    S = [0.0]; Sm = [0.0]; Pc = [0.0]
    pnoC = rho0.copy()                    # nhánh chưa có bản ghi Carnot (không chuẩn hóa)
    for k in range(1, nk + 1):
        Sm.append(Sm[-1] + marginal_entropy(Ks, r))
        r = apply(Ks, r)
        pnoC = Ks[0] @ pnoC @ Ks[0].conj().T + Ks[2] @ pnoC @ Ks[2].conj().T
        S.append(vn(r)); Pc.append(1 - np.real(np.trace(pnoC)))
    S = np.array(S); Sm = np.array(Sm); dS = np.diff(S)
    print(f"{'k':>5} {'t/T0':>6} {'S(E)=S(ρ_S)':>12} {'P(≥1 Carnot)':>13} {'Σ_i S(E_i)':>11}")
    for k in (0, 1, 3, 7, 25, 50, 125, 250, 500, 1000, 1500):
        print(f"{k:5d} {k*TAU/T0:6.2f} {S[k]:12.4f} {Pc[k]:13.4f} {Sm[k]:11.2f}")
    print(f"   số bước ΔS(E) < 0: {int(np.sum(dS < -1e-12))}/{nk}; ΔS nhỏ nhất {dS.min():.2e}; "
          f"S max {S.max():.4f} tại k = {int(S.argmax())}, cuối {S[-1]:.4f}; trần ln d = {math.log(H.shape[0]):.3f}")
    rate = (Sm[-1] - Sm[-501]) / (500 * TAU / T0)
    print(f"   entropy bản ghi thô Σ_i S(E_i): tăng tuyến tính {rate:.2f} nat/T0, không bão hòa")

    # phản ví dụ: chỉ kênh Carnot (làm lạnh Davies, nhiệt độ 0), trạng thái đầu thuần ở mức kích thích cao
    Hd, Vd = np.linalg.eigh(H)
    for lab, e in (("trạng thái riêng thứ 6", 5), ("trạng thái riêng thứ 12", 11)):
        v = Vd[:, e]; rc = np.outer(v, v.conj()); Kc = kraus(H, [Ls[0]], TAU)
        Sc = [0.0]
        for k in range(1, 2501):
            rc = apply(Kc, rc); Sc.append(vn(rc))
        Sc = np.array(Sc); km = int(Sc.argmax())
        print(f"   phản ví dụ (chỉ Carnot, T = 0), xuất phát {lab}: S(E) 0 → {Sc.max():.3f} (k = {km}, {km*TAU/T0:.1f} T0)"
              f" → {Sc[-1]:.2e} (k = 2500); bước ΔS < 0: {int(np.sum(np.diff(Sc) < -1e-12))}")
    print("   ⟹ S(E) của toàn môi trường = S(ρ_S) (Schmidt), bị chặn bởi ln d và KHÔNG đơn điệu nói chung;"
          " 'ΔS(E)/Δk ≥ 0' không phải định luật hai.")
    print("   Đại lượng đơn điệu đúng: entropy bản ghi thô Σ_i S(E_i) (mỗi tích tắc +S(E_k) ≥ 0, tăng tuyến tính) — "
          "= ' đồng hồ entropy của bể' ở nam-dong-ho-luong-tu-kiem-chung.md\n")
    return S


def _cmul(X, Y):
    """Nhân ma trận phức bằng 4 gemm thực (zgemm trên máy này chậm hơn ~10 lần)."""
    Xr, Xi = np.ascontiguousarray(X.real), np.ascontiguousarray(X.imag)
    Yr, Yi = np.ascontiguousarray(Y.real), np.ascontiguousarray(Y.imag)
    return (Xr @ Yr - Xi @ Yi) + 1j * (Xr @ Yi + Xi @ Yr)


def sandwich(K, rn):
    """K ρ_n K† cho cả chồng ρ_n."""
    n, d, _ = rn.shape
    A = _cmul(rn.reshape(n * d, d), K.conj().T).reshape(n, d, d)                       # ρ K†
    B = _cmul(A.transpose(0, 2, 1).reshape(n * d, d), K.T).reshape(n, d, d)           # (K ρ K†)ᵀ
    return B.transpose(0, 2, 1)


def step5(H, Ls, rho0):
    print("=== Bước 5. Kim = bộ đếm Carnot N_C thay cho số tích tắc k ===")
    Ks = kraus(H, Ls, TAU); d = H.shape[0]
    NMX, nk = 60, 1500
    rn = np.zeros((NMX + 1, d, d), complex); rn[0] = rho0
    Pn = np.zeros((nk + 1, NMX + 1)); Pn[0, 0] = 1
    acc = np.zeros((NMX + 1, d, d), complex); acc += rn
    for k in range(1, nk + 1):
        new = sandwich(Ks[0], rn) + sandwich(Ks[2], rn)
        new[1:] += sandwich(Ks[1], rn[:-1])
        new[np.abs(new) < 1e-200] = 0          # số subnormal làm phép nhân chậm ~100 lần
        rn = new; acc += rn
        Pn[k] = np.real(np.einsum('naa->n', rn))
    t = np.arange(nk + 1) * TAU / T0
    mean_n = Pn @ np.arange(NMX + 1)
    r = (mean_n[-1] - mean_n[-501]) / (500 * TAU / T0)
    print(f"   tốc độ bản ghi Carnot dừng r = {r:.3f}/T0 (ideal k: 25/T0); Fano(N_C) tại 60 T0 = "
          f"{(Pn[-1] @ np.arange(NMX+1)**2 - mean_n[-1]**2)/mean_n[-1]:.2f}")
    print(f"{'n':>4} {'⟨t|n⟩ (T0)':>11} {'σ_t|n (T0)':>11} {'σ_t/T0':>7} {'‖ρ_S|n − ρ_Lind(⟨t⟩)‖':>22} {'‖ρ_S(k) − ρ_Lind(kτ)‖ cùng ⟨t⟩':>30}")
    ns = (0, 1, 2, 3, 5, 10, 20, 30); rows = []
    for n in ns:
        w = Pn[:, n]; w = w / w.sum()
        tm = w @ t; ts = math.sqrt(max(w @ t**2 - tm**2, 0)); rows.append((n, tm, ts))
    # các thời điểm cần nghiệm Lindblad: ⟨t|n⟩ và kτ gần nhất
    need = sorted({tm * T0 for _, tm, _ in rows} | {int(round(tm * 25)) * TAU for _, tm, _ in rows})
    ref = dict(zip(need, lind_at(H, Ls, rho0, need)))
    rk, kdone = rho0.copy(), 0
    for n, tm, ts in rows:
        rc = acc[n] / np.real(np.trace(acc[n]))
        rl = ref[tm * T0]
        k = int(round(tm * 25))
        if k < kdone:
            rk, kdone = rho0.copy(), 0
        for _ in range(k - kdone):
            rk = apply(Ks, rk)
        kdone = k
        rlk = ref[k * TAU]
        print(f"{n:4d} {tm:11.2f} {ts:11.2f} {ts:7.2f} {tdist(rc, rl):22.3f} {tdist(rk, rlk):30.2e}")
    print("   ⟹ n = 0, 1: ρ_S|n lệch Lindblad 0.1–0.3; n ≥ 2: lệch ≈ 0.02–0.03, ngang sai số O(τ) của kim lý tưởng — nhưng chỉ"
          " vì ρ_S đã gần trạng thái dừng (đồng hồ nào cũng 'khớp' khi hệ không còn đổi)")
    print("   ⟹ phép thử quyết định là độ phân giải: mỗi vạch n nhòe σ_t = 2–7.6 T0 ≫ 0.071 T0 (đòi hỏi 99%) và > chu trình 5 T0"
          " (hàng n = 30 bị cắt bởi cửa sổ 60 T0)")
    print("   ⟹ hai vai trò trong chứng minh tách rời: k (lý tưởng, đều, τ = T0/25) tái lập Lindblad; N_C (vật lý) đơn điệu"
          " nhưng thô — không đồng thời là cùng một kim\n")


def main():
    ls, H, Ls = channels()
    rho0 = np.outer(psi0(), psi0().conj())
    print(f"S: {H.shape[0]} chiều; τ = T0/25; L_c (Carnot/Davies), L_p (bơm)\n")
    step1(H, Ls)
    step3(H, Ls, rho0)
    step4(H, Ls, rho0)
    step5(H, Ls, rho0)


if __name__ == "__main__":
    main()
