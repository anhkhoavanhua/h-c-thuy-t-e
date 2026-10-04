"""Số liệu cho Part III (dạng lượng tử vượt thời gian của E_dual, HocThuyetE (7).pdf, mô hình V.C.1, đồng hồ V.C.4).

Mô hình S: rô-to ℓ ≤ 8 ⊗ hạt n ≤ 6 (119 chiều, đúng V.C.1), H Bảng I với −g x sinθ; kênh L_c = √0.5 û₋ (Davies),
L_p = √0.2 e^{iθ} f(ℓ̂), ℓ0 = 4. Đồng hồ chuỗi xung τ = T0/25; môi trường: mỗi xung một ancilla 3 mức (mô hình va chạm
đơn nhất U, giai_so_wdw_pw_e7). Trạng thái lịch sử |Ψ⟩⟩ = N^{-1/2} Σ_k |k⟩ ⊗ U_k⋯U_1|ψ0⟩|0…0⟩.

(N1) Thoát bùng nổ 3^k một cách CHÍNH XÁC: mọi đại lượng người quan sát trong đo được chỉ phụ thuộc
       ρ_S(k) = Tr_E⟨k|Ψ⟩⟩ = Λ^k(ρ0),  Λ(ρ) = Σ_j K_j ρ K_j†,  K_j = (1⊗⟨j|) U (1⊗|0⟩)
     và thống kê bản ghi (nhánh). Kiểm Λ^k ≡ trạng thái lịch sử đầy đủ tới k = 7 (119·3^7 chiều), rồi kéo dài tới
     60 T0 (1500 xung); hội tụ τ → 0 về Lindblad (giới hạn liên tục của WdW); S(E) = S(ρ_S) theo k.
(N2) Bản ghi môi trường = nhánh của |Ψ⟩⟩, lấy mẫu chính xác theo Born (unraveling của cùng ánh xạ Kraus):
     trung bình nhánh tái lập Λ^k; đợt bản ghi Carnot, năng lượng mỗi bản ghi, thời gian chờ; nhánh luôn thuần.
(N3) Đồng hồ thực: rô-to tác dụng J = Iθ̇ (ħ) làm đồng hồ Page–Wootters cho hạt (g = 1, Ω = 0.3) — độ trùng với
     Schrödinger theo J (1 vòng và giới hạn WdW), xác định J tối thiểu để thời gian "cổ điển".
Chạy: python3 giai_so_part3.py [N1|N2|N3]   (N1 ≈ 3 phút, N2 ≈ 1 phút, N3 ≈ 2 phút)
"""
import math
import sys

import numpy as np

import giai_so_wdw_pw_e7 as W

T0 = 2 * math.pi
W.LMAX, W.NMAX = 8, 6


def setup(tau):
    W.TAU = tau
    ls, H, Ls = W.channels()
    U, d = W.collision_unitary(H, Ls)
    U4 = U.reshape(d, 3, d, 3)
    K = [U4[:, j, :, 0] for j in range(3)]
    return ls, H, Ls, U, d, K


def start_state(ls, d):
    """Bó sóng thỏa ràng buộc (pha Quá khứ, con trỏ của bể Carnot): ℓ0 = 3, θ0 = 0, hạt kết hợp α = iκℓ0/√2."""
    cr = np.exp(-(ls - 3)**2 / (4 * 0.7**2)); cr = cr / np.linalg.norm(cr)
    al = 1j * (1 / 3) * 3 / math.sqrt(2)
    n = np.arange(W.NMAX + 1)
    cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in n]) * math.exp(-abs(al)**2 / 2)
    return np.kron(cr, cx / np.linalg.norm(cx))


def lam(K, r):
    return sum(k @ r @ k.conj().T for k in K)


def td(a, b):
    return 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh(a - b))))


def obs(ls, d):
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(W.NMAX + 1))
    return Lop


def part_n1():
    tau = T0 / 25
    ls, H, Ls, U, d, K = setup(tau)
    psi0 = start_state(ls, d)
    rho0 = np.outer(psi0, psi0.conj())
    print(f"=== N1. Trạng thái lịch sử đầy đủ vs ánh xạ rút gọn (S {d} chiều, τ = T0/25) ===")
    print(f"Σ K_j†K_j − 1: {np.linalg.norm(sum(k.conj().T @ k for k in K) - np.eye(d)):.1e} (bảo toàn vết chính xác)")
    NT = 7
    psi = np.zeros((d,) + (3,) * NT, complex)
    psi[(slice(None),) + (0,) * NT] = psi0
    r = rho0.copy()
    worst = 0.0
    for k in range(NT):
        psi = W.apply_tick(psi, U, d, k)
        r = lam(K, r)
        M = psi.reshape(d, -1)
        worst = max(worst, td(M @ M.conj().T, r))
    print(f"k ≤ {NT}: max ‖Tr_E⟨k|Ψ⟩⟩ − Λ^k(ρ0)‖₁/2 = {worst:.1e}  (trạng thái lịch sử đủ {d}·3^{NT} = {d*3**NT} chiều)")

    print("\nHội tụ τ → 0 về Lindblad liên tục (giới hạn liên tục của ràng buộc WdW), khoảng cách vết tại t:")
    Lop = obs(ls, d)
    ref = {}
    rr = rho0.copy(); tprev = 0.0
    for tT in (1.0, 2.0, 4.0):
        seg = (tT - tprev) * T0
        rr = W.lindblad_ref(H, Ls, rr, seg, dt=seg / math.ceil(seg / 0.005)); tprev = tT   # bước chia đúng đoạn
        ref[tT] = rr.copy()
    print(f"{'τ':>9} {'t = 1 T0':>11} {'2 T0':>9} {'4 T0':>9}")
    for m in (25, 50, 100, 200):
        _, H2, Ls2, U2, d2, K2 = setup(T0 / m)
        r = rho0.copy(); out = []; k = 0
        for tT in (1.0, 2.0, 4.0):
            while k < round(tT * m):
                r = lam(K2, r); k += 1
            out.append(td(r, ref[tT]))
        print(f"{'T0/'+str(m):>9} {out[0]:11.2e} {out[1]:9.2e} {out[2]:9.2e}")

    print("\nKéo dài tới 60 T0 (1500 xung) bằng Λ^k — không cần lưu môi trường:")
    _, H, Ls, U, d, K = setup(tau)
    r = rho0.copy()
    print(f"{'t/T0':>6} {'⟨ℓ⟩':>7} {'⟨H⟩ (ħω0)':>10} {'Tr ρ_S²':>8} {'S(ρ_S) = S(E) (nat)':>20}")
    for k in range(0, 1501):
        if k in (0, 25, 50, 125, 250, 500, 1000, 1500):
            print(f"{k/25:6.0f} {np.real(np.trace(Lop @ r)):7.3f} {np.real(np.trace(H @ r)):10.3f} "
                  f"{np.real(np.trace(r @ r)):8.4f} {W.vn(r):20.4f}")
        r = lam(K, r)
    print()


def part_n2(M=200, ticks=1500):
    tau = T0 / 25
    ls, H, Ls, U, d, K = setup(tau)
    psi0 = start_state(ls, d)
    Lop = obs(ls, d)
    rng = np.random.default_rng(5)
    lfin, recs, dEc, purity_dev = [], [], [], 0.0
    for b in range(M):
        p = psi0.copy(); rc = []
        for k in range(ticks):
            outs = [kk @ p for kk in K]
            w = np.array([np.vdot(o, o).real for o in outs])
            j = rng.choice(3, p=w / w.sum())
            Eb = np.real(np.vdot(p, H @ p))
            p = outs[j] / math.sqrt(w[j])
            if j == 1:
                rc.append(k)
                if k > ticks // 5:
                    dEc.append(Eb - np.real(np.vdot(p, H @ p)))
        purity_dev = max(purity_dev, abs(np.vdot(p, p).real - 1))
        lfin.append(np.real(np.vdot(p, Lop @ p)))
        recs.append(rc)
    r = np.outer(psi0, psi0.conj())
    for k in range(ticks):
        r = lam(K, r)
    print(f"=== N2. Bản ghi môi trường = nhánh của |Ψ⟩⟩ ({M} nhánh lấy mẫu Born × {ticks} xung = {ticks/25:.0f} T0) ===")
    print(f"⟨ℓ⟩ cuối: trung bình nhánh {np.mean(lfin):.3f} ± {np.std(lfin)/math.sqrt(M):.3f}, Λ^k: "
          f"{np.real(np.trace(Lop @ r)):.3f}; mỗi nhánh thuần (lệch chuẩn hóa ≤ {purity_dev:.1e})")
    dE = np.array(dEc)
    print(f"năng lượng mỗi bản ghi Carnot: TB {dE.mean():.3f} ħω0, trung vị {np.median(dE):.3f}, 10–90% "
          f"{np.percentile(dE,10):.3f}–{np.percentile(dE,90):.3f} (V.C.1: ≈ 2.05)")
    sizes, waits = [], []
    for rc in recs:
        rc = [k for k in rc if k > ticks // 5]
        cur, last = 0, None
        for k in rc:
            if last is None or (k - last) * tau < T0:
                cur += 1
            else:
                sizes.append(cur); waits.append((k - last) * tau / T0); cur = 1
            last = k
    s, w = np.array(sizes), np.array(waits)
    print(f"đợt bản ghi Carnot (khoảng cách < T0): {len(s)} đợt, {s.mean():.2f} bản ghi/đợt (V.C.1: 1.57), "
          f"chờ giữa đợt TB {w.mean():.2f} T0, CV {w.std()/w.mean():.2f}")
    print(f"pha theo bản ghi: Quá khứ = khoảng không bản ghi ≥ T0, Hiện tại = đợt, Tương lai = sau đợt; thứ tự trên đồng hồ "
          f"chuỗi xung là đơn điệu theo định nghĩa (bộ đếm k)\n")


def clock_fidelity(J, g=1.0, Om=0.3, nmax=15):
    """Rô-to I = J/Ω, ℓ0 = J làm đồng hồ cho hạt; trả min F trên t ∈ {0.1..0.4} T_rot cho cửa sổ 1 vòng và WdW."""
    from giai_so_luong_tu import rotor_ops, ho_x
    I = J / Om
    Trot = 2 * math.pi / Om
    sl = max(1.0, math.sqrt(J / 3))
    ls = np.arange(J - 24, J + 25)
    cos_t, sin_t = rotor_ops(ls)
    x = ho_x(nmax)
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    H = np.kron(np.diag(ls**2 / (2 * I)), np.eye(nmax + 1)) + np.kron(np.eye(len(ls)), Hx) - g * np.kron(sin_t, x)
    E, Uv = np.linalg.eigh(H)
    c = np.exp(-(ls - J)**2 / (4 * sl * sl)).astype(complex); c /= np.linalg.norm(c)
    chi0 = np.zeros(nmax + 1, complex); chi0[0] = 1
    amp = Uv.conj().T @ np.kron(c, chi0)
    E0 = float(E[np.argmax(np.abs(amp))])
    res = {}
    for label, tw in (("1 vòng", 0.35 * Trot), ("WdW", 50 * Trot)):
        Psi = (Uv @ (np.exp(-((E - E0) * tw)**2 / 2) * amp)).reshape(len(ls), nmax + 1)
        Fs = []
        for f in (0.1, 0.2, 0.3, 0.4):
            t = f * Trot
            chi = np.exp(1j * ls * Om * t) @ Psi; chi /= np.linalg.norm(chi)
            psi = chi0.copy(); n_sub = 1500; dt = t / n_sub
            for i in range(n_sub):
                w_, v_ = np.linalg.eigh(Hx - g * math.sin(Om * (i + 0.5) * dt) * x)
                psi = v_ @ (np.exp(-1j * w_ * dt) * (v_.conj().T @ psi))
            Fs.append(abs(np.vdot(psi, chi))**2)
        res[label] = min(Fs)
    return res


def part_n3():
    print("=== N3. Đồng hồ rô-to có tác dụng J (ħ): độ trùng Page–Wootters với Schrödinger (g = 1, Ω = 0.3) ===")
    print(f"{'J (ħ)':>6} {'1 − F_min (1 vòng)':>19} {'1 − F_min (WdW)':>16}")
    rows = []
    for J in (3, 6, 12, 30, 60):
        r = clock_fidelity(J)
        rows.append((J, 1 - r["1 vòng"], 1 - r["WdW"]))
        print(f"{J:6d} {1-r['1 vòng']:19.4f} {1-r['WdW']:16.4f}", flush=True)
    a = np.array(rows)
    s1 = np.polyfit(np.log(a[:, 0]), np.log(np.maximum(a[:, 1], 1e-9)), 1)[0]
    s2 = np.polyfit(np.log(a[:, 0]), np.log(np.maximum(a[:, 2], 1e-9)), 1)[0]
    print(f"luật lũy thừa 1 − F ∝ J^p: p = {s1:.2f} (1 vòng), {s2:.2f} (WdW); mô hình ở J ≈ 3ħ\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["N1", "N2", "N3"]:
        {"N1": part_n1, "N2": part_n2, "N3": part_n3}[name.upper()]()
