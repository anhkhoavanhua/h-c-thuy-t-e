"""Bài toán người quan sát (measurement / observer problem) cho học thuyết E_dual (HocThuyetE (7).pdf).

Người quan sát = môi trường nhận năng lượng Carnot. Va chạm Carnot truyền xung Λ = m_eff u dọc pháp tuyến ∇f, nên môi
trường "đo" vận tốc tương đối û = (κ/I)·½{cosθ, p_θ} − p_x (giai_so_u_hermite_drude.basis).

(A) Einselection — "predictability sieve" (Zurek): mô hình lượng tử E_dual (rô-to ℓ ≤ 5 ⊗ hạt n ≤ 4, H Bảng I với
    −g x sinθ), đo liên tục û với cường độ Γ:  dρ/dt = −i[H, ρ] + Γ(ûρû − ½{û², ρ}).
    So độ tinh khiết Tr ρ² sau nửa chu kỳ T0/2 của ba họ trạng thái đầu thuần:
      E = 10 trạng thái riêng năng lượng thấp nhất (sóng dừng),
      U = 10 trạng thái riêng của û có |trị riêng| nhỏ nhất,
      C = 8 bó sóng "cổ điển" (rô-to định xứ quanh θ0 với ℓ0 ⊗ hạt kết hợp tại x = κ sinθ0, v = κℓ0 cosθ0) — hình ảnh hạt.
    Họ nào giữ thuần tốt nhất = trạng thái con trỏ (cái người quan sát "thấy").
(B) Quantum Darwinism: hệ = 3 mức rô-to của pha Quá khứ (|ℓ| = 1, 2, 3; trọng số Born p = 0.2, 0.3, 0.5), môi trường
    = N = 10 nucleon-qubit, mỗi qubit k ghi mức j bằng một phép quay không hoàn hảo góc j·φ_k (φ_k ∈ [0.25, 0.75] rad).
    Trạng thái toàn cục thuần (đơn nhất, không sụp đổ). Đo: hệ số kết hợp ngoài chéo của ρ_S theo số qubit; thông tin
    tương hỗ I(S:F) theo kích thước mảnh F; độ dư thừa R_δ (δ = 0.1); quy tắc Born khi một người quan sát đọc m qubit.
Chạy: python3 giai_so_nguoi_quan_sat.py [A|A2|D|B|B2]   (A ≈ 3 phút, A2 ≈ 1 phút, D ≈ 15 phút, B/B2 ≈ 20 s; B2: ghi mạnh hơn φ ∈ [0.6, 1.0])
"""
import itertools
import math
import random
import sys

import numpy as np

from giai_so_u_hermite_drude import basis

T0 = 2 * math.pi
KAP = 1 / 3


def purity_after(states, H, u, G, T):
    """RK4 theo lô cho nhiều ma trận mật độ cùng lúc."""
    u2 = u @ u
    rho = np.array([np.outer(s, s.conj()) for s in states])
    lam = np.max(np.abs(np.linalg.eigvalsh(u)))
    dt = min(0.005, 0.2 / max(G * lam**2, 1e-9))
    n = int(math.ceil(T / dt))
    dt = T / n

    def L(r):
        return -1j * (H @ r - r @ H) + G * (u @ r @ u - 0.5 * (u2 @ r + r @ u2))
    for _ in range(n):
        k1 = L(rho); k2 = L(rho + dt / 2 * k1); k3 = L(rho + dt / 2 * k2); k4 = L(rho + dt * k3)
        rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return np.real(np.einsum('kij,kji->k', rho, rho))


def coherent_vec(alpha, nmax):
    n = np.arange(nmax + 1)
    c = np.array([alpha**k / math.sqrt(math.factorial(k)) for k in n], complex) * math.exp(-abs(alpha)**2 / 2)
    return c / np.linalg.norm(c)


def families(lmax, nmax):
    """Ba họ trạng thái đầu: E (riêng năng lượng), U (riêng û), C (bó sóng cổ điển định xứ)."""
    ls, cos_t, sin_t, pth, H, u = basis(lmax=lmax, nmax=nmax)
    E, V = np.linalg.eigh(H)
    w, Vu = np.linalg.eigh(u)
    fam = {"E": [V[:, j] for j in range(10)], "U": [Vu[:, j] for j in np.argsort(np.abs(w))[:10]], "C": []}
    for th0 in (0, math.pi / 2, math.pi, 3 * math.pi / 2):
        for l0 in (2, 3):
            cr = np.exp(-(ls - l0)**2 / (4 * 0.7**2) - 1j * ls * th0)
            cr /= np.linalg.norm(cr)
            al = (KAP * math.sin(th0) + 1j * KAP * l0 * math.cos(th0)) / math.sqrt(2)
            fam["C"].append(np.kron(cr, coherent_vec(al, nmax)))
    return H, u, E, fam


def part_a():
    LMAX, NMAX = 5, 4
    H, u, E, fam = families(LMAX, NMAX)
    gaps = np.diff(E[:11])
    print(f"=== A. Predictability sieve: đo liên tục û với cường độ Γ (rô-to ℓ ≤ {LMAX} ⊗ hạt n ≤ {NMAX}), t = T0/2 ===")
    print(f"khoảng cách mức thấp: TB {gaps.mean():.3f}, min {gaps.min():.3f} ħω0")
    print(f"{'Γ (ω0)':>8} {'E: sóng dừng TB (min)':>22} {'U: trạng thái riêng û TB (min)':>31} "
          f"{'C: bó sóng cổ điển TB (min)':>28} {'thắng':>8}")
    for G in (0.02, 0.2, 0.5, 2.0, 20.0):
        res = {k: purity_after(f, H, u, G, T0 / 2) for k, f in fam.items()}
        win = max(res, key=lambda k: res[k].mean())
        print(f"{G:8.2f} {res['E'].mean():13.4f} ({res['E'].min():.4f}) {res['U'].mean():22.4f} ({res['U'].min():.4f}) "
              f"{res['C'].mean():19.4f} ({res['C'].min():.4f}) {win:>8}", flush=True)
    print("(Γ = 0.5 = tốc độ làm nguội Carnot Γc của mô hình V.C.1)\n")


def part_a2():
    print("=== A2. Sieve theo thời gian tại Γ = 0.5, kiểm cắt cụt cơ sở ===")
    for LM, NM in ((5, 4), (7, 6)):
        H, u, E, fam = families(LM, NM)
        var = {k: np.mean([np.real(v.conj() @ u @ u @ v) - np.real(v.conj() @ u @ v)**2 for v in f])
               for k, f in fam.items()}
        print(f"cơ sở ℓ ≤ {LM}, n ≤ {NM}: Var(û) TB  E {var['E']:.3f}  U {var['U']:.3f}  C {var['C']:.3f}")
        for tt in (0.02, 0.05, 0.1, 0.25):
            r = {k: purity_after(f, H, u, 0.5, tt * T0).mean() for k, f in fam.items()}
            print(f"  t = {tt:5.2f} T0: E {r['E']:.4f}  U {r['U']:.4f}  C {r['C']:.4f}  thắng {max(r, key=r.get)}")
    print()


def purity_after_ops(states, H, Ls, Ts):
    """Lindblad tổng quát (kênh Davies, bơm); một lần tích phân, trả độ tinh khiết TB tại các mốc Ts."""
    S = sum(L.conj().T @ L for L in Ls)
    rho = np.array([np.outer(s, s.conj()) for s in states])
    lam = max(np.max(np.linalg.eigvalsh(S)), np.max(np.abs(np.linalg.eigvalsh(H))))
    dt0 = min(0.01, 1.0 / lam)                      # RK4 ổn định tới ~2.8/λ

    def Lv(r):
        d = -1j * (H @ r - r @ H) - 0.5 * (S @ r + r @ S)
        for L in Ls:
            d = d + L @ r @ L.conj().T
        return d
    out, t = [], 0.0
    for T in Ts:
        n = max(1, int(math.ceil((T - t) / dt0)))
        dt = (T - t) / n
        for _ in range(n):
            k1 = Lv(rho); k2 = Lv(rho + dt / 2 * k1); k3 = Lv(rho + dt / 2 * k2); k4 = Lv(rho + dt * k3)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = T
        out.append(np.real(np.einsum('kij,kji->k', rho, rho)))
    return out


def davies_ops(lmax, nmax, H, u, Gc, Gp):
    """Kênh làm nguội Davies √Γc û₋ (như V.C.1) và (tùy chọn) bơm √Γp e^{iθ} f(ℓ̂), ℓ0 = 4."""
    E, V = np.linalg.eigh(H)
    Ls = [math.sqrt(Gc) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)]
    if Gp > 0:
        ls = np.arange(-lmax, lmax + 1)
        up = np.diag(np.ones(len(ls) - 1), -1).astype(complex)
        F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
        Ls.append(math.sqrt(Gp) * np.kron(up @ F, np.eye(nmax + 1)))
    return Ls


def part_d():
    print("=== D. Predictability sieve với kênh thật của V.C.1: làm nguội Davies √Γc û₋ (± bơm √0.2 e^{iθ} f(ℓ̂)) ===")
    for LM, NM in ((5, 4), (7, 6)):
        H, u, E, fam = families(LM, NM)
        print(f"-- cơ sở ℓ ≤ {LM}, n ≤ {NM}")
        Lc = davies_ops(LM, NM, H, u, 0.5, 0.0)[0]
        LL = Lc.conj().T @ Lc
        for k, f in fam.items():
            e = np.mean([np.real(v.conj() @ H @ v) for v in f])
            r = np.mean([np.real(v.conj() @ LL @ v) for v in f])
            ub = np.mean([abs(np.real(v.conj() @ u @ v)) for v in f])
            print(f"   họ {k}: ⟨H⟩ TB {e:.3f} ħω0, tốc độ nhảy Carnot ⟨L†L⟩ TB {r:.3f} (Γc = 0.5), |⟨û⟩| TB {ub:.3f}")
        for Gc, Gp in ((0.5, 0.0), (0.5, 0.2), (0.05, 0.0), (5.0, 0.0)):
            Ls = davies_ops(LM, NM, H, u, Gc, Gp)
            tts = (0.1, 0.25, 0.5, 1.0)
            res = {k: purity_after_ops(f, H, Ls, [tt * T0 for tt in tts]) for k, f in fam.items()}
            row = []
            for i, tt in enumerate(tts):
                r = {k: res[k][i].mean() for k in res}
                row.append(f"t={tt:4.2f}T0: E {r['E']:.3f} C {r['C']:.3f} U {r['U']:.3f} → {max(r, key=r.get)}")
            print(f"  Γc = {Gc:4.2f}, Γp = {Gp:3.1f}:  " + " | ".join(row), flush=True)
    print("(E = sóng dừng/trạng thái riêng năng lượng, C = bó sóng cổ điển định xứ, U = trạng thái riêng û)\n")


def entropy(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-13]
    return float(-np.sum(w * np.log2(w)))


def red_entropy(Psi, keep, N):
    """Entropy của tập trục 'keep' (0 = hệ, 1..N = qubit) của trạng thái thuần Psi có hình (3, 2, ..., 2)."""
    rest = [k for k in range(N + 1) if k not in keep]
    M = np.transpose(Psi, keep + rest).reshape(int(np.prod([Psi.shape[k] for k in keep])), -1)
    return entropy(M @ M.conj().T) if M.shape[0] <= M.shape[1] else entropy(M.conj().T @ M)


def part_b(N=10, seed=7, lo=0.25, hi=0.75):
    rng = np.random.default_rng(seed)
    p = np.array([0.2, 0.3, 0.5])
    phi = rng.uniform(lo, hi, N)
    rec = [[np.array([math.cos(j * f), math.sin(j * f)]) for f in phi] for j in range(3)]
    Psi = np.zeros((3,) + (2,) * N, complex)
    for j in range(3):
        v = np.array([1.0])
        for k in range(N):
            v = np.kron(v, rec[j][k])
        Psi[j] = math.sqrt(p[j]) * v.reshape((2,) * N)
    print(f"=== B. Quantum Darwinism: 3 mức rô-to (|ℓ| = 1, 2, 3; p = {p.tolist()}) + {N} nucleon-qubit, "
          f"φ_k ∈ [{lo}, {hi}] rad ===")
    print(f"trạng thái toàn cục thuần: S(toàn bộ) = {red_entropy(Psi, list(range(N + 1)), N):.2e} bit (đơn nhất, không sụp đổ)")
    print(f"{'số qubit đã ghi':>16} {'|ρ_S(1,2)|':>11} {'|ρ_S(2,3)|':>11} {'|ρ_S(1,3)|':>11}")
    for m in (0, 1, 2, 4, 6, 10):
        ov = lambda a, b: np.prod([rec[a][k] @ rec[b][k] for k in range(m)]) if m else 1.0
        print(f"{m:16d} {abs(math.sqrt(p[0]*p[1])*ov(0,1)):11.4f} {abs(math.sqrt(p[1]*p[2])*ov(1,2)):11.4f} "
              f"{abs(math.sqrt(p[0]*p[2])*ov(0,2)):11.4f}")
    HS = float(-np.sum(p * np.log2(p)))
    print(f"H(S) = entropy Born cổ điển = {HS:.4f} bit")
    print(f"{'|F|':>4} {'I(S:F) TB (bit)':>16} {'I/H(S)':>7}")
    Ivals = []
    for m in range(0, N + 1):
        frs = list(itertools.combinations(range(1, N + 1), m))
        random.Random(m).shuffle(frs)
        Is = []
        for F in frs[:6]:
            F = list(F)
            Is.append(red_entropy(Psi, [0], N) + (red_entropy(Psi, F, N) if F else 0.0) - red_entropy(Psi, [0] + F, N))
        Ivals.append(np.mean(Is))
        print(f"{m:4d} {Ivals[-1]:16.4f} {Ivals[-1]/HS:7.3f}")
    md = next((m for m, I in enumerate(Ivals) if I >= 0.9 * HS), None)
    print(f"m_δ (I ≥ 0.9 H) = {md} qubit ⟹ độ dư thừa R_δ = N/m_δ = {N/md:.1f} người quan sát độc lập có cùng thông tin\n")

    print("Quy tắc Born: một người quan sát đọc m qubit (cơ sở |0⟩,|1⟩), xác suất mỗi chuỗi bit tính TRỰC TIẾP từ |Ψ⟩⟩")
    print("(lấy vết hệ + các qubit còn lại; không giả định nhánh), đoán mức bằng hợp lý cực đại")
    for m in (1, 3, 5, 10):
        PB = np.sum(np.abs(Psi)**2, axis=tuple(range(m + 1, N + 1)))      # (3, 2^m): P(j, chuỗi bit)
        PB = PB.reshape(3, -1)
        Pbits = PB.sum(axis=0)
        guess = np.argmax(PB, axis=0)                                      # MAP theo P(j, bit) = p_j P(bit|j)
        freq = [float(Pbits[guess == k].sum()) for k in range(3)]
        acc = float(sum(PB[guess[i], i] for i in range(PB.shape[1])))
        print(f"  m = {m:2d}: tần suất kết luận = {[round(f,3) for f in freq]}  (Born |c_j|² = {p.tolist()}); "
              f"xác suất kết luận khớp nhánh = {acc:.1%}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "A2", "B", "B2"]:
        {"A": part_a, "A2": part_a2, "D": part_d, "B": part_b, "B2": lambda: part_b(lo=0.6, hi=1.0)}[name.upper()]()
