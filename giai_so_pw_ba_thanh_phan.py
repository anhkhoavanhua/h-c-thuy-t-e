"""Page–Wootters ba thành phần (Đồng hồ C + Hệ S + Môi trường E) cho mô hình lượng tử E_dual.

Mô hình S (như Outlook của bản thảo): rô-to ℓ ∈ [−8, 8] ⊗ dao động tử n ≤ 6; H với ghép −g x sinθ (dấu đúng);
làm nguội Carnot L_c = √0.5 û₋ (Davies T=0); bơm pha L_p = √0.2 e^{iθ} f(ℓ̂), ℓ0 = 4 (chế độ ⟨ℓ⟩ ≈ 2.4).
E = môi trường nhận các lượng tử (dilation của Lindblad / mô hình va chạm).

Góc nhìn 1 — đồng hồ lý tưởng đếm xung (C không tương tác): trạng thái điều kiện ρ_S(t) = Tr_E ⟨t|Ψ⟩⟩ = nghiệm Lindblad.
   Đo độ tinh khiết Tr ρ² của (rô-to+hạt) và của riêng hạt theo t; có "sụt đột ngột" không.
Góc nhìn 2 — điều kiện hóa thêm bản ghi của E (quỹ đạo lượng tử MCWF): S luôn thuần; đo độ tinh khiết của
   RIÊNG HẠT ngay trước/sau mỗi bước nhảy làm nguội (vướng víu rô-to–hạt).
Góc nhìn 3 — rô-to là đồng hồ (Page–Wootters như bản thảo): ở trạng thái dừng, độ tinh khiết của ρ_hạt(θ0) theo θ0.
Đơn vị ħ = m = ω0 = 1. Chạy: python3 giai_so_pw_ba_thanh_phan.py [1|2|3]
"""
import math
import random
import sys

import numpy as np

from giai_so_u_hermite_drude import basis

NMAX = 6


def model():
    ls, _, _, _, H, u = basis(lmax=8, nmax=NMAX)
    dC, dS = len(ls), NMAX + 1
    E, V = np.linalg.eigh(H)
    Lc = math.sqrt(0.5) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    up = np.diag(np.ones(dC - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(dS))
    return ls, dC, dS, H, [Lc, Lp]


def lind(rho, H, Ls):
    d = -1j * (H @ rho - rho @ H)
    for L in Ls:
        LdL = L.conj().T @ L
        d += L @ rho @ L.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return d


def ptrace_rotor(rho, dC, dS):
    return np.einsum('anam->nm', rho.reshape(dC, dS, dC, dS))


def purity(r):
    return float(np.real(np.trace(r @ r)))


def vn(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


def start(ls, dC, dS):
    rho = np.zeros((dC * dS, dC * dS), complex)
    i0 = int(np.where(ls == 0)[0][0]) * dS
    rho[i0, i0] = 1
    return rho


def view1():
    ls, dC, dS, H, Ls = model()
    rho, t, dt = start(ls, dC, dS), 0.0, 0.01
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(dS))
    print("=== Góc nhìn 1: đồng hồ lý tưởng (đếm xung), ρ_S(t) = Tr_E ⟨t|Ψ⟩⟩ ===")
    print(f"{'t':>6} {'t/T0':>6} {'⟨ℓ⟩':>6} {'Tr ρ²(rô-to+hạt)':>17} {'S_vN':>6} {'Tr ρ²(hạt)':>11} {'S_vN(hạt)':>10}")
    prev, maxdrop, tmax = 1.0, 0.0, 0.0
    marks = [0, 0.5, 1, 2, 4, 8, 16, 32, 64, 120]
    for mt in marks:
        while t < mt - 1e-9:
            k1 = lind(rho, H, Ls); k2 = lind(rho + dt / 2 * k1, H, Ls)
            k3 = lind(rho + dt / 2 * k2, H, Ls); k4 = lind(rho + dt * k3, H, Ls)
            rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            p = purity(rho)
            if (prev - p) / dt > maxdrop:
                maxdrop, tmax = (prev - p) / dt, t
            prev = p
        rp = ptrace_rotor(rho, dC, dS)
        print(f"{mt:6g} {mt/(2*math.pi):6.2f} {float(np.real(np.trace(Lop @ rho))):6.2f} {purity(rho):17.4f} "
              f"{vn(rho):6.3f} {purity(rp):11.4f} {vn(rp):10.3f}")
    print(f"tốc độ giảm độ tinh khiết lớn nhất: {maxdrop:.4f} mỗi đv t, tại t = {tmax:.2f} "
          f"(= {tmax/(2*math.pi):.2f} T0)\n")
    return rho


def view2(Tmax=400.0, dt=0.01, seed=11):
    ls, dC, dS, H, Ls = model()
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    rng = random.Random(seed)
    psi = np.zeros(dC * dS, complex)
    psi[int(np.where(ls == 0)[0][0]) * dS] = 1
    t, r = 0.0, rng.random()
    before, after, between = [], [], []
    k = 0

    def ppur(v):
        v = v / np.linalg.norm(v)
        M = v.reshape(dC, dS)
        return purity(M.T @ M.conj())

    while t < Tmax:
        psi = U @ psi
        t += dt
        k += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            j = 0 if rng.random() < ws[0] / sum(ws) else 1
            pb = ppur(psi)
            psi = Ls[j] @ psi
            psi /= np.linalg.norm(psi)
            if j == 0 and t > 20:
                before.append(pb); after.append(ppur(psi))
            r = rng.random()
        elif t > 20 and k % 200 == 0:
            between.append(ppur(psi))
    b, a, m = np.array(before), np.array(after), np.array(between)
    print(f"=== Góc nhìn 2: điều kiện hóa cả bản ghi môi trường (1 quỹ đạo MCWF, {Tmax:g} đv t) ===")
    print("độ tinh khiết của (rô-to+hạt) = 1 tại mọi thời điểm (trạng thái thuần theo định nghĩa)")
    print(f"độ tinh khiết của RIÊNG HẠT: giữa các bước nhảy TB {m.mean():.4f} (min {m.min():.4f}); "
          f"ngay trước bước làm nguội TB {b.mean():.4f}; ngay sau TB {a.mean():.4f} "
          f"({len(b)} bước làm nguội)\n")


def view3(rho):
    ls, dC, dS, H, Ls = model()
    a = np.diag(np.sqrt(np.arange(1, dS)), 1)
    x = (a + a.T) / math.sqrt(2)
    R = rho.reshape(dC, dS, dC, dS)
    ths = np.linspace(0, 2 * math.pi, 64, endpoint=False)
    pur, pw = [], []
    for th in ths:
        c = np.exp(1j * ls * th)
        rs = np.einsum('a,anbm,b->nm', c, R, c.conj())
        w = np.trace(rs).real
        rs /= w
        pur.append(purity(rs)); pw.append(w)
    pur, pw = np.array(pur), np.array(pw)
    print("=== Góc nhìn 3: rô-to là đồng hồ, trạng thái dừng (t = 120), ρ_hạt(θ0) ===")
    print(f"Tr ρ_hạt²(θ0): min {pur.min():.4f}, TB {np.average(pur, weights=pw):.4f}, max {pur.max():.4f}; "
          f"p(θ0) max/min = {pw.max()/pw.min():.2f}; θ0 và θ0+2π cho cùng ρ_hạt (tuần hoàn)\n")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3"]
    rho = None
    if "1" in which or "3" in which:
        rho = view1()
    if "2" in which:
        view2()
    if "3" in which:
        view3(rho)
