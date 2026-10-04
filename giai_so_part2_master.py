"""Số liệu cho Part II (động lực học hệ lượng tử mở của E_dual, HocThuyetE (7).pdf, mô hình V.C.1).

H0 = p_θ²/2I − ΔE cosθ + (a†a + ½) − g x sinθ (Bảng I, dấu đúng), kênh:
  L_p = √0.2 e^{iθ} f(ℓ̂), f = √max(0, 1 − ℓ²/16)   (bơm pha bão hòa, đúng V.C.1)
  L_c = √0.5 û₋                                    (Davies, làm nguội Carnot)
  L_x = √γx a, γx = 0.7                             (ma sát quang học, kênh MỚI đề xuất cho Part II; chạy có/không)

(M1) Liouvillian chính xác: L vec(ρ) = 0 ⟹ ρ_ss; quan sát dừng, tốc độ nhảy từng kênh, khe Liouvillian (tốc độ hồi phục
     chậm nhất); hội tụ cơ sở (ℓ ≤ 5, n ≤ 4: 3025²; ℓ ≤ 6, n ≤ 5: 5929²).
(M2) MCWF 100 quỹ đạo (cơ sở V.C.1 ℓ ≤ 8, n ≤ 6): phân rã bước nhảy (ℓ→ℓ−1, n→n−1), phổ năng lượng Carnot, số bước/đợt,
     thời gian chờ, độ tinh khiết riêng hạt trước/sau bước Carnot.
(M3) Phi Markov tại k_BT = 0.1 ħω0, γ = 0.7: ⟨x²⟩, ⟨p²⟩ của Lindblad (L = √γ a), Caldeira–Leggett (Markov nhiệt độ cao)
     và nghiệm chính xác chuyển động Brown lượng tử (Feynman–Vernon / Ford–Lewis–O'Connell, bể Drude ω_c).
(M4) Darwinism động lực trong mô hình va chạm thật (giai_so_wdw_pw_e7): qubit tham chiếu R vướng với hai nhánh của S;
     I(R:F) theo kích thước mảnh bản ghi F — so cặp bó sóng ràng buộc (C) với cặp trạng thái riêng năng lượng (E).
Chạy: python3 giai_so_part2_master.py [M1|M2|M3|M4]   (M1 ≈ 3 phút, M2 ≈ 4 phút, M3 vài giây, M4 ≈ 1 phút)
"""
import math
import random
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

from giai_so_u_hermite_drude import basis

T0 = 2 * math.pi
GX = 0.7


def model(lmax, nmax, with_lx):
    ls, cos_t, sin_t, pth, H, u = basis(lmax=lmax, nmax=nmax)
    d = len(ls)
    E, V = np.linalg.eigh(H)
    A = math.sqrt(0.5) * np.triu(V.conj().T @ u @ V, 1)
    Lc = V @ A @ V.conj().T
    up = np.diag(np.ones(d - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(nmax + 1))
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
    Ls = {"bơm": Lp, "Carnot": Lc}
    if with_lx:
        Ls["L_x"] = math.sqrt(GX) * np.kron(np.eye(d), a)
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(nmax + 1))
    Nn = np.kron(np.eye(d), a.conj().T @ a)
    return dict(ls=ls, d=d, H=H, E=E, V=V, A=A, Ls=Ls, Lop=Lop, Nn=Nn, nmax=nmax)


# ------------------------------------------------------------------ M1
def liouvillian(H, Ls):
    """Siêu toán tử theo vec hàng: vec(AρB) = (A ⊗ Bᵀ) vec(ρ)."""
    D = H.shape[0]
    I = np.eye(D)
    L = -1j * (np.kron(H, I) - np.kron(I, H.T))
    for Lj in Ls:
        LdL = Lj.conj().T @ Lj
        L += np.kron(Lj, Lj.conj()) - 0.5 * np.kron(LdL, I) - 0.5 * np.kron(I, LdL.T)
    return L


def steady(H, Ls):
    D = H.shape[0]
    L = liouvillian(H, Ls)
    b = np.zeros(D * D, complex)
    L[0, :] = np.eye(D).reshape(-1)               # thay một phương trình bằng Tr ρ = 1
    b[0] = 1
    rho = np.linalg.solve(L, b).reshape(D, D)
    return (rho + rho.conj().T) / 2


def part_m1():
    print("=== M1. Trạng thái dừng chính xác của phương trình chủ (Liouvillian) ===")
    print(f"{'cơ sở':>12} {'kênh':>18} {'⟨ℓ⟩':>7} {'⟨H⟩':>7} {'⟨n⟩':>7} {'Tr ρ²':>7} {'S (nat)':>8} "
          f"{'nhảy Carnot/T0':>15} {'nhảy bơm/T0':>12} {'nhảy L_x/T0':>12} {'λ_min ≠ 0':>10} {'min eig ρ':>10}")
    for lmax, nmax in ((5, 4), (6, 5)):
        for lx in (False, True):
            m = model(lmax, nmax, lx)
            rho = steady(m["H"], list(m["Ls"].values()))
            ev = lambda O: float(np.real(np.trace(O @ rho)))
            w = np.linalg.eigvalsh(rho); wp = w[w > 1e-14]
            rates = {k: ev(L.conj().T @ L) * T0 for k, L in m["Ls"].items()}
            if lmax == 5:
                lam = np.linalg.eigvals(liouvillian(m["H"], list(m["Ls"].values())))
                lam = lam[np.argsort(-lam.real)]
                gap = f"{-lam[1].real:10.4f}"
            else:
                gap = f"{'—':>10}"
            print(f"{f'ℓ≤{lmax},n≤{nmax}':>12} {'V.C.1 + L_x' if lx else 'V.C.1':>18} {ev(m['Lop']):7.3f} {ev(m['H']):7.3f} "
                  f"{ev(m['Nn']):7.4f} {np.real(np.trace(rho @ rho)):7.4f} "
                  f"{-np.sum(wp*np.log(wp)):8.4f} {rates['Carnot']:15.3f} {rates['bơm']:12.3f} "
                  f"{rates.get('L_x', 0.0):12.4f} {gap} {w.min():10.1e}", flush=True)
    print("(λ_min ≠ 0: phần thực của trị riêng Liouvillian thứ hai = tốc độ hồi phục chậm nhất, đơn vị ω0;"
          " so với mô hình va chạm τ = T0/25 của Part III: ⟨ℓ⟩ = 2.262, Tr ρ² = 0.241)\n")


# ------------------------------------------------------------------ M2
def _traj(args):
    seed, Tmax, dt, lx = args
    m = model(8, 6, lx)
    H, V, A = m["H"], m["V"], m["A"]
    names = list(m["Ls"]); Ls = [m["Ls"][k] for k in names]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    d, nm = m["d"], m["nmax"]
    rng = random.Random(seed)
    psi = np.zeros(d * (nm + 1), complex)
    psi[int(np.where(m["ls"] == 0)[0][0]) * (nm + 1)] = 1
    t, r = 0.0, rng.random()
    jumps, trans, purs = [], defaultdict(float), []

    def ppur(v):
        M = (v / np.linalg.norm(v)).reshape(d, nm + 1)
        rr = M.T @ M.conj()
        return float(np.real(np.trace(rr @ rr)))
    while t < Tmax:
        psi = U @ psi
        t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            k = int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random()))
            k = min(k, len(Ls) - 1)
            Eb = float(np.real(np.vdot(psi, H @ psi)))
            if names[k] == "Carnot" and t > Tmax / 5:
                c = V.conj().T @ psi
                Wt = np.abs(A * c[None, :])**2; Wt /= Wt.sum()
                for i, j in zip(*np.nonzero(Wt > 1e-6)):
                    trans[(int(i), int(j))] += float(Wt[i, j])
                pb = ppur(psi)
            psi = Ls[k] @ psi
            psi /= np.linalg.norm(psi)
            if t > Tmax / 5:
                jumps.append((t, names[k], Eb - float(np.real(np.vdot(psi, H @ psi)))))
                if names[k] == "Carnot":
                    purs.append((pb, ppur(psi)))
            r = rng.random()
    return jumps, dict(trans), purs


def part_m2(N=100, Tmax=1500.0, dt=0.01):
    m = model(8, 6, False)
    E, V = m["E"], m["V"]
    lab_l = np.array([float(np.real(V[:, j].conj() @ np.abs(m["Lop"]) @ V[:, j])) for j in range(len(E))])
    lab_n = np.array([float(np.real(V[:, j].conj() @ m["Nn"] @ V[:, j])) for j in range(len(E))])
    print(f"=== M2. MCWF {N} quỹ đạo × {Tmax:g} đv t (≈ {Tmax/T0:.0f} T0), cơ sở ℓ ≤ 8, n ≤ 6 ===")
    with Pool() as pool:
        for lx in (False, True):
            res = pool.map(_traj, [(2000 + i, Tmax, dt, lx) for i in range(N)])
            span = N * 0.8 * Tmax / T0
            allj = [j for rr in res for j in rr[0]]
            cnt = defaultdict(int)
            for j in allj:
                cnt[j[1]] += 1
            cE = np.array([j[2] for j in allj if j[1] == "Carnot"])
            trans = defaultdict(float)
            for rr in res:
                for kij, v in rr[1].items():
                    trans[kij] += v
            cls = defaultdict(lambda: [0.0, 0.0])
            for (i, j), wgt in trans.items():
                dl = round(lab_l[j]) - round(lab_l[i]); dn = round(lab_n[j]) - round(lab_n[i])
                key = ("ℓ→ℓ−1" if (dl == 1 and dn == 0) else "n→n−1" if (dl == 0 and dn == 1) else
                       "ℓ→ℓ−2" if (dl == 2 and dn == 0) else "khác")
                cls[key][0] += wgt; cls[key][1] += wgt * (E[j] - E[i])
            tot = sum(v[0] for v in cls.values()); etot = sum(v[1] for v in cls.values())
            sizes, waits = [], []
            for rr in res:
                cj = [j[0] for j in rr[0] if j[1] == "Carnot"]
                cur, last = 0, None
                for tj in cj:
                    if last is None or tj - last < T0:
                        cur += 1
                    else:
                        sizes.append(cur); waits.append((tj - last) / T0); cur = 1
                    last = tj
            s, wts = np.array(sizes), np.array(waits)
            P = np.array([p for rr in res for p in rr[2]])
            print(f"-- {'V.C.1 + L_x (γx = 0.7)' if lx else 'V.C.1'}: tốc độ nhảy/T0 " +
                  ", ".join(f"{k} {v/span:.3f}" for k, v in cnt.items()))
            print(f"   năng lượng mỗi bước Carnot: TB {cE.mean():.3f} ħω0, trung vị {np.median(cE):.3f}; phổ: " +
                  ", ".join(f"{a_}–{b_}: {np.mean((cE>=a_)&(cE<b_)):.1%}" for a_, b_ in ((0, .5), (.5, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 9))))
            print("   phân rã: " + ", ".join(f"{k} {cls[k][0]/tot:.1%} số bước / {cls[k][1]/etot:.1%} năng lượng (khe TB "
                                            f"{cls[k][1]/cls[k][0]:.2f})" for k in ("ℓ→ℓ−1", "ℓ→ℓ−2", "n→n−1", "khác") if k in cls))
            print(f"   đợt Carnot (khoảng cách < T0): {len(s)} đợt, TB {s.mean():.2f} bước/đợt, phân bố " +
                  ", ".join(f"{k}:{np.mean(s==k):.0%}" for k in range(1, 5)) + f"; chờ TB {wts.mean():.2f} T0, CV {wts.std()/wts.mean():.2f}")
            print(f"   độ tinh khiết riêng hạt quanh bước Carnot: trước TB {P[:,0].mean():.4f}, sau TB {P[:,1].mean():.4f}; "
                  f"tăng ở {np.mean(P[:,1] > P[:,0]):.0%} số bước", flush=True)
    print()


# ------------------------------------------------------------------ M3
def part_m3(kT=0.1, gam=0.7):
    from giai_so_bom_pha import qbm_moments
    print(f"=== M3. Phi Markov: dao động tử ħω0, γ = {gam}, k_BT = {kT} ħω0 (1 MeV) ===")
    nth = 1 / math.expm1(1 / kT)
    print(f"{'phương trình':>40} {'⟨x²⟩':>8} {'⟨p²⟩':>8} {'⟨x²⟩⟨p²⟩':>10} {'ν = √tích':>9} {'trị riêng ρ nhỏ nhất':>21}")

    def row(name, x2, p2):
        nu = math.sqrt(x2 * p2)
        lam1 = (1 / (nu + 0.5)) * ((nu - 0.5) / (nu + 0.5))      # trị riêng n = 1 của trạng thái Gauss nhiệt
        lmin = min(1 / (nu + 0.5), lam1)
        print(f"{name:>40} {x2:8.4f} {p2:8.4f} {x2*p2:10.4f} {nu:9.4f} {lmin:21.4f}")
    row("Lindblad L = √γ a (+ n_th)", nth + 0.5, nth + 0.5)
    row("Caldeira–Leggett (Markov, nhiệt độ cao)", kT, kT)
    for E in (4.5, 50, 100):
        x2, p2 = qbm_moments(E / 10, kT, gam)
        row(f"chính xác QBM, Drude ħω_c = {E:g} MeV", x2, p2)
    print("(Robertson–Schrödinger: ⟨x²⟩⟨p²⟩ ≥ 1/4; CL vi phạm ⟹ ρ có trị riêng âm; Lindblad thỏa nhưng bỏ qua ép"
          " do liên kết mạnh; QBM chính xác là chuẩn phi Markov cho phần tuyến tính — tương đương HEOM hội tụ)\n")


# ------------------------------------------------------------------ M4
def part_m4():
    import giai_so_wdw_pw_e7 as W
    W.LMAX, W.NMAX = 4, 3
    NT = 7
    ls, H, Ls = W.channels()
    E, V = np.linalg.eigh(H)
    n = np.arange(W.NMAX + 1)

    def packet(th0, l0):
        cr = np.exp(-(ls - l0)**2 / (4 * 0.7**2) - 1j * ls * th0); cr /= np.linalg.norm(cr)
        al = ((1 / 3) * math.sin(th0) + 1j * (1 / 3) * l0 * math.cos(th0)) / math.sqrt(2)
        cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in n]) * math.exp(-abs(al)**2 / 2)
        return np.kron(cr, cx / np.linalg.norm(cx))
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(W.NMAX + 1))
    il = [j for j in range(len(E)) if abs(np.real(V[:, j].conj() @ Lop @ Lop @ V[:, j]) - 9) < 1.0]
    pairs = {"C: bó sóng ràng buộc θ0 = 0 / π (ℓ0 = 3)": (packet(0, 3), packet(math.pi, 3)),
             "E: trạng thái riêng cơ bản / mức |ℓ| = 3": (V[:, 0], V[:, il[0]])}

    def ent(M):
        w = np.linalg.eigvalsh(M); w = w[w > 1e-13]
        return float(-np.sum(w * np.log2(w)))
    print(f"=== M4. Darwinism động lực trong mô hình va chạm thật (S 36 chiều, {NT} bản ghi 3 mức) ===")
    for tau_div in (25, 5):
        W.TAU = T0 / tau_div
        U, d = W.collision_unitary(H, Ls)
        print(f"-- τ = T0/{tau_div} ({NT} bản ghi = {NT/tau_div:.2f} T0)")
        for name, (A_, B_) in pairs.items():
            psi = np.zeros((d,) + (3,) * NT + (2,), complex)
            psi[(slice(None),) + (0,) * NT + (0,)] = A_ / math.sqrt(2)
            psi[(slice(None),) + (0,) * NT + (1,)] = B_ / math.sqrt(2)
            for k in range(NT):
                psi = W.apply_tick(psi, U, d, k)
            Iv = []
            for m_ in range(0, NT + 1):
                F = list(range(1, m_ + 1))                       # mảnh = m bản ghi đầu
                rest = [ax for ax in range(1, NT + 1) if ax not in F]
                def red(keep):
                    other = [ax for ax in range(NT + 2) if ax not in keep]
                    Mx = np.transpose(psi, keep + other).reshape(int(np.prod([psi.shape[k] for k in keep])), -1)
                    return ent(Mx @ Mx.conj().T) if Mx.shape[0] <= Mx.shape[1] else ent(Mx.conj().T @ Mx)
                SR = red([NT + 1]); SF = red(F) if F else 0.0; SRF = red([NT + 1] + F)
                Iv.append(SR + SF - SRF)
            ov = abs(np.vdot(A_, B_))
            print(f"   {name}: ⟨A|B⟩ = {ov:.3f}; I(R:F) theo |F| = 0…{NT} (bit): " + " ".join(f"{x:.3f}" for x in Iv) +
                  f";  I(R:F)/H(R) cả môi trường = {Iv[-1]:.3f}", flush=True)
    print("(H(R) = 1 bit; Darwinism đòi I(R:F) ≈ 1 bit đã ở mảnh nhỏ (bình nguyên); I ≈ 0 = môi trường không biết nhánh nào)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["M1", "M2", "M3", "M4"]:
        {"M1": part_m1, "M2": part_m2, "M3": part_m3, "M4": part_m4}[name.upper()]()
