"""Hai câu hỏi mở của đồng hồ photon trên học thuyết E_dual (bản thảo HocThuyetE (5).pdf).

Câu 1 — Page–Wootters lượng tử có môi trường (C + S + E), mô hình S như Outlook của bản thảo, cơ sở thu nhỏ
  (rô-to ℓ ∈ [−5, 5], hạt n ≤ 4): H với −g x sinθ; làm nguội Carnot L_c = √0.5 û₋ (Davies T=0);
  bơm pha L_p = √0.2 e^{iθ} f(ℓ̂), ℓ0 = 4.
  (1a) Kiểm tường minh: trạng thái lịch sử Σ_k |k⟩_C ⊗ |ψ_SE(k)⟩ (E = thanh ghi 3 mức mỗi nhịp: không nhảy /
       nhảy Carnot / nhảy bơm) cho ρ_S(k) = Tr_E trùng nghiệm Lindblad (3 nhịp).
  (1b) Đồng hồ lý tưởng đếm xung: độ tinh khiết của (rô-to+hạt) và của riêng hạt theo t.
  (1c) Điều kiện hóa cả bản ghi E (quỹ đạo lượng tử): độ tinh khiết của riêng hạt trước/sau bước nhảy Carnot.
Câu 2 — đồng hồ đếm vòng:
  (2a) OAM ℓ ∈ [−L, L]: trạng thái góc |φ⟩ tuần hoàn 2π ⟹ không đếm vòng; số thời điểm phân biệt được = 2L+1.
  (2b) Chuỗi xung chu kỳ T_rep: dùng mô phỏng cổ điển Bảng I (1 MeV) — tỉ lệ đợt mà lần gãy (Hiện tại) và
       lần tái khóa (Tương lai) rơi cùng một xung ⟹ không phân biệt được thứ tự.
Chạy: python3 giai_so_pw_luong_tu_mo.py [1|2]   (1 ≈ 20 s; 2 ≈ 40 s trên 8 lõi)
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

from giai_so_u_hermite_drude import basis

LMAX, NMAX = 5, 4


def model():
    ls, _, _, _, H, u = basis(lmax=LMAX, nmax=NMAX)
    dC, dS = len(ls), NMAX + 1
    E, V = np.linalg.eigh(H)
    Lc = math.sqrt(0.5) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    up = np.diag(np.ones(dC - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(dS))
    return ls, dC, dS, H, [Lc, Lp]


def purity(r):
    return float(np.real(np.trace(r @ r)))


def vn(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


def part1():
    ls, dC, dS, H, Ls = model()
    D = dC * dS
    psi0 = np.zeros(D, complex); psi0[int(np.where(ls == 0)[0][0]) * dS] = 1
    # ---- (1a) kiểm tường minh bằng thanh ghi môi trường, bước dt ----
    dt = 0.05
    K0 = np.eye(D) - 1j * H * dt - 0.5 * dt * sum(L.conj().T @ L for L in Ls)
    Ks = [K0] + [math.sqrt(dt) * L for L in Ls]           # toán tử Kraus mỗi nhịp; E ghi chỉ số Kraus
    branches = {(): psi0}
    rho_l = np.outer(psi0, psi0.conj())
    print("=== Câu 1a. Trạng thái lịch sử C⊗S⊗E tường minh (3 nhịp, E = thanh ghi 3 mức/nhịp) ===")
    for k in range(1, 4):
        branches = {key + (j,): Kj @ v for key, v in branches.items() for j, Kj in enumerate(Ks)}
        rho_hist = sum(np.outer(v, v.conj()) for v in branches.values())     # Tr_E ⟨k|Ψ⟩⟩
        rho_l = sum(Kj @ rho_l @ Kj.conj().T for Kj in Ks)
        print(f"nhịp k={k}: {len(branches)} nhánh E; ‖ρ_S(lịch sử) − ρ_S(Lindblad rời rạc)‖ = "
              f"{np.linalg.norm(rho_hist - rho_l):.2e}; Tr ρ_S² = {purity(rho_hist/np.trace(rho_hist).real):.6f}")
    print("(mỗi nhánh E là thuần; chỉ khi lấy vết trên E thì S mới hỗn hợp)\n")

    # ---- (1b) đồng hồ lý tưởng đếm xung: Lindblad liên tục ----
    def lind(r):
        d = -1j * (H @ r - r @ H)
        for L in Ls:
            LdL = L.conj().T @ L
            d += L @ r @ L.conj().T - 0.5 * (LdL @ r + r @ LdL)
        return d
    rho, t, h = np.outer(psi0, psi0.conj()), 0.0, 0.01
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(dS))
    print("=== Câu 1b. Đồng hồ đếm xung lý tưởng: ρ_S(t) = Tr_E ⟨t|Ψ⟩⟩ ===")
    print(f"{'t/T0':>6} {'⟨ℓ⟩':>6} {'Tr ρ²(rô-to+hạt)':>17} {'S_vN':>6} {'Tr ρ²(hạt)':>11} {'S_vN(hạt)':>10}")
    prev, maxrate, tmax = 1.0, 0.0, 0.0
    for mt in (0, 0.1, 0.25, 0.5, 1, 2, 4, 8, 16):
        T = mt * 2 * math.pi
        while t < T - 1e-9:
            k1 = lind(rho); k2 = lind(rho + h / 2 * k1); k3 = lind(rho + h / 2 * k2); k4 = lind(rho + h * k3)
            rho = rho + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += h
            p = purity(rho)
            if (prev - p) / h > maxrate:
                maxrate, tmax = (prev - p) / h, t
            prev = p
        rp = np.einsum('anam->nm', rho.reshape(dC, dS, dC, dS))
        print(f"{mt:6g} {float(np.real(np.trace(Lop @ rho))):6.2f} {purity(rho):17.4f} {vn(rho):6.3f} "
              f"{purity(rp):11.4f} {vn(rp):10.3f}")
    print(f"tốc độ giảm độ tinh khiết lớn nhất {maxrate:.3f}/đv t tại t = {tmax/(2*math.pi):.3f} T0 "
          f"(giảm TRƠN, không có bước nhảy)\n")

    # ---- (1c) điều kiện hóa cả bản ghi E: quỹ đạo lượng tử ----
    dtt = 0.01
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dtt)) @ np.linalg.inv(Rm)
    rng = random.Random(3)

    def ppur(v):
        M = (v / np.linalg.norm(v)).reshape(dC, dS)
        return purity(M.T @ M.conj())
    psi, t, r = psi0.copy(), 0.0, rng.random()
    bef, aft, betw, step = [], [], [], 0
    while t < 400:
        psi = U @ psi; t += dtt; step += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            j = 0 if rng.random() < ws[0] / sum(ws) else 1
            pb = ppur(psi)
            psi = Ls[j] @ psi; psi /= np.linalg.norm(psi)
            if j == 0 and t > 20:
                bef.append(pb); aft.append(ppur(psi))
            r = rng.random()
        elif t > 20 and step % 100 == 0:
            betw.append(ppur(psi))
    b, a, m = np.array(bef), np.array(aft), np.array(betw)
    print("=== Câu 1c. Điều kiện hóa cả bản ghi môi trường (1 quỹ đạo, 400 đv t ≈ 64 T0) ===")
    print("Tr ρ²(rô-to+hạt) = 1 mọi lúc (thuần)")
    print(f"Tr ρ²(riêng hạt): giữa các bước nhảy TB {m.mean():.4f} (min {m.min():.4f}); ngay trước bước Carnot TB "
          f"{b.mean():.4f}; ngay sau TB {a.mean():.4f}; {len(b)} bước Carnot\n")


def part2a():
    print("=== Câu 2a. Đồng hồ OAM ===")
    for L in (1, 5, 20):
        ls = np.arange(-L, L + 1)
        ket = lambda phi: np.exp(-1j * ls * phi) / math.sqrt(len(ls))
        same = abs(np.vdot(ket(0.3), ket(0.3 + 2 * math.pi)))**2
        orth = [abs(np.vdot(ket(0), ket(2 * math.pi * k / len(ls))))**2 for k in range(1, len(ls))]
        print(f"ℓ ∈ [−{L}, {L}]: |⟨φ|φ+2π⟩|² = {same:.6f}; số thời điểm trực giao trong một vòng = {len(ls)} "
              f"(max chồng lấp giữa chúng {max(orth):.1e}); độ phân giải 2π/{len(ls)}")
    print("⟹ OAM tăng độ phân giải TRONG một vòng nhưng vẫn tuần hoàn 2π: không đếm vòng\n")


def _traj(a):
    import giai_so_dong_co_langevin as Lm
    from giai_so_kramers_A import setup
    cfg, seed, T = a
    setup(cfg)
    return Lm.run(0.1, 1.05 * cfg[2], 0.1, T, 2e-3, seed, gx=0.7, w_start=0.82 * cfg[2])['events']


def part2b(N=8):
    from giai_so_burst_g import cfg_for
    T0 = 2 * math.pi
    cfg = cfg_for(0.1)
    with Pool() as pool:
        evs = pool.map(_traj, [(cfg, 61000 + i, 40 * T0) for i in range(N)])
    durs, gaps = [], []
    for ev in evs:
        for i, e in enumerate(ev):
            if e[0] != 'gãy' or e[1] < 5 * T0:
                continue
            prev = next((x for x in reversed(ev[:i]) if x[0] == 'khóa'), None)
            if prev is None or e[1] - prev[1] < T0:
                continue
            for j in range(i + 1, len(ev)):
                if ev[j][0] == 'khóa':
                    gaps.append(ev[j][1] - ev[j - 1][1])          # gãy → khóa ngay sau
                    nxt = next((y for y in ev[j + 1:] if y[0] == 'gãy'), None)
                    if nxt is not None and nxt[1] - ev[j][1] >= T0:
                        durs.append((ev[j][1] - e[1]) / T0)
                        break
    durs, gaps = np.array(durs), np.array(gaps) / T0
    print(f"=== Câu 2b. Đồng hồ chuỗi xung (Bảng I, 1 MeV; {len(durs)} đợt, {len(gaps)} cặp gãy→khóa) ===")
    print(f"thời lượng đợt: TB {durs.mean():.2f} T0, min {durs.min():.3f} T0; khoảng gãy→khóa kế tiếp: TB "
          f"{gaps.mean():.3f} T0, min {gaps.min():.4f} T0")
    print(f"{'T_rep (T0)':>11} {'xung/vòng rô-to':>16} {'đợt không xếp được thứ tự':>26} {'cặp gãy→khóa không phân biệt':>30}")
    for Tr in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
        # xác suất cùng rơi vào một xung với pha lệch ngẫu nhiên: max(0, 1 − Δ/T_rep)
        pd = float(np.mean(np.maximum(0, 1 - durs / Tr)))
        pg = float(np.mean(np.maximum(0, 1 - gaps / Tr)))
        print(f"{Tr:11.2f} {1/(Tr*2.575/(2*math.pi)):16.1f} {pd:26.1%} {pg:30.1%}")
    print("(xung/vòng tính với θ̇_d ≈ 2.575; xác suất 'cùng một xung' lấy trung bình trên pha ngẫu nhiên của chuỗi xung)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["1", "2"]:
        if name == "1":
            part1()
        else:
            part2a(); part2b()
