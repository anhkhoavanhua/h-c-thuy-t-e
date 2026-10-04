"""Tách các bước nhảy làm nguội của mô hình MCWF V.C.1 (HocThuyetE (7).pdf) thành phần rô-to (ℓ → ℓ−1) và phần
hạt (n → n−1).

Mô hình y hệt V.C.1 / giai_so_bom_pha.py cấu hình 'sat': rô-to ℓ ∈ [−8, 8] ⊗ dao động tử n ≤ 6, H với −g x sinθ,
bơm L_p = √0.2 e^{iθ} f(ℓ̂) (ℓ0 = 4), làm nguội Davies L_c = √0.5 û₋; 16 quỹ đạo × 1500 đv t (≈ 240 T0).
Ở mỗi bước làm nguội ghi:
  (1) ΔE toàn phần và phần thay đổi của từng khối: rô-to ⟨ℓ²/2I − ΔE cosθ⟩, hạt ⟨n + ½⟩, ghép ⟨−g x sinθ⟩;
      Δ⟨|ℓ|⟩, Δ⟨n⟩.
  (2) phân rã theo chuyển dời riêng j → i: trọng số |A_ij c_j|² (A = L_c trong cơ sở riêng của H), mỗi trạng
      thái riêng gắn nhãn (|ℓ|, n) trội; xếp chuyển dời vào lớp theo (Δ|ℓ|, Δn).
So với thang rô-to trên đa tạp ràng buộc (schrodinger-rang-buoc-kiem-chung.md: ℓ 3→2 = 2.32 ħω0) và với ħω0.
Chạy: python3 giai_so_mcwf_tach_buoc.py   (≈ 2 phút, 8 lõi)
"""
import math
import random
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

from giai_so_u_hermite_drude import basis, DE, G

NMAX = 6
T0 = 2 * math.pi


def model():
    ls, cos_t, sin_t, pth, H, u = basis()
    d = len(ls)
    E, V = np.linalg.eigh(H)
    A = math.sqrt(0.5) * np.triu(V.conj().T @ u @ V, 1)            # L_c trong cơ sở riêng: j → i < j
    Lc = V @ A @ V.conj().T
    up = np.diag(np.ones(d - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Lp = math.sqrt(0.2) * np.kron(up @ F, np.eye(NMAX + 1))
    Ip = np.eye(NMAX + 1)
    a = np.diag(np.sqrt(np.arange(1, NMAX + 1)), 1)
    x = (a + a.T) / math.sqrt(2)
    Hrot = np.kron(pth @ pth / 2 - DE * cos_t, Ip)
    Hpar = np.kron(np.eye(d), a.T @ a + 0.5 * Ip)
    Hcpl = -G * np.kron(sin_t, x)
    Labs = np.kron(np.diag(np.abs(ls).astype(float)), Ip)
    Nn = np.kron(np.eye(d), np.diag(np.arange(NMAX + 1, dtype=float)))
    return dict(ls=ls, d=d, H=H, E=E, V=V, A=A, Ls=[Lp, Lc], Hrot=Hrot, Hpar=Hpar, Hcpl=Hcpl, Labs=Labs, Nn=Nn)


def _traj(args):
    seed, Tmax, dt = args
    m = model()
    H, Ls, V, A = m["H"], m["Ls"], m["V"], m["A"]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    ops = [m[k] for k in ("H", "Hrot", "Hpar", "Hcpl", "Labs", "Nn")]
    ev = lambda p, O: float(np.real(np.vdot(p, O @ p)))
    rng = random.Random(seed)
    psi = np.zeros(m["d"] * (NMAX + 1), complex)
    psi[int(np.where(m["ls"] == 0)[0][0]) * (NMAX + 1)] = 1
    t, r = 0.0, rng.random()
    recs, trans = [], defaultdict(float)
    while t < Tmax:
        psi = U @ psi
        t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            if k == 1 and t > Tmax / 5:
                before = [ev(psi, O) for O in ops]
                c = V.conj().T @ psi
                W = np.abs(A * c[None, :])**2                     # trọng số j → i
                W /= W.sum()
                for i, j in zip(*np.nonzero(W > 1e-6)):
                    trans[(int(i), int(j))] += float(W[i, j])
            psi = Ls[k] @ psi
            psi /= np.linalg.norm(psi)
            if k == 1 and t > Tmax / 5:
                after = [ev(psi, O) for O in ops]
                recs.append([b - a for b, a in zip(before, after)])
            r = rng.random()
    return recs, dict(trans)


def main(N=16, Tmax=1500.0, dt=0.01):
    m = model()
    E, V = m["E"], m["V"]
    lab_l = np.array([float(np.real(V[:, j].conj() @ m["Labs"] @ V[:, j])) for j in range(len(E))])
    lab_n = np.array([float(np.real(V[:, j].conj() @ m["Nn"] @ V[:, j])) for j in range(len(E))])
    with Pool() as pool:
        res = pool.map(_traj, [(900 + i, Tmax, dt) for i in range(N)])
    recs = np.array([r for rr, _ in res for r in rr])
    trans = defaultdict(float)
    for _, tr in res:
        for kij, v in tr.items():
            trans[kij] += v
    nj = len(recs)
    print(f"=== MCWF V.C.1 ({N} quỹ đạo × {Tmax:g} đv t), {nj} bước làm nguội ===")
    names = ["ΔE toàn phần", "Δ rô-to ⟨ℓ²/2I − ΔE cosθ⟩", "Δ hạt ⟨n + ½⟩", "Δ ghép ⟨−g x sinθ⟩", "Δ⟨|ℓ|⟩", "Δ⟨n⟩"]
    print("(1) thay đổi kỳ vọng qua mỗi bước (dương = mất đi):")
    for k, nm in enumerate(names):
        col = recs[:, k]
        print(f"   {nm:>28}: TB {col.mean():+.3f}, trung vị {np.median(col):+.3f}, 10–90% "
              f"{np.percentile(col,10):+.3f}…{np.percentile(col,90):+.3f}")
    fr = recs[:, 1].sum() / recs[:, 0].sum()
    fp = recs[:, 2].sum() / recs[:, 0].sum()
    print(f"   tỉ phần năng lượng mất: rô-to {fr:.1%}, hạt {fp:.1%}, ghép {1-fr-fp:.1%}")

    print("(2) phân rã theo chuyển dời riêng j → i (trọng số |A_ij c_j|²):")
    cls = defaultdict(lambda: [0.0, 0.0])
    gaps = []
    for (i, j), wgt in trans.items():
        dl = round(lab_l[j]) - round(lab_l[i])
        dn = round(lab_n[j]) - round(lab_n[i])
        key = ("rô-to ℓ→ℓ−1" if (dl == 1 and dn == 0) else "hạt n→n−1" if (dl == 0 and dn == 1) else
               "rô-to ℓ→ℓ−2" if (dl == 2 and dn == 0) else "trộn / khác")
        cls[key][0] += wgt
        cls[key][1] += wgt * (E[j] - E[i])
        gaps.append((E[j] - E[i], wgt))
    tot = sum(v[0] for v in cls.values())
    print(f"   {'lớp':>14} {'tỉ phần bước':>13} {'khe TB (ħω0)':>13} {'tỉ phần năng lượng':>19}")
    etot = sum(v[1] for v in cls.values())
    for key in ("rô-to ℓ→ℓ−1", "rô-to ℓ→ℓ−2", "hạt n→n−1", "trộn / khác"):
        if key in cls:
            wsum, esum = cls[key]
            print(f"   {key:>14} {wsum/tot:13.1%} {esum/wsum:13.3f} {esum/etot:19.1%}")
    g = np.array(gaps)
    hist, edges = np.histogram(g[:, 0], bins=[0, 0.5, 1.5, 2.5, 3.5, 5, 20], weights=g[:, 1])
    print("   phổ khe năng lượng (trọng số): " + ", ".join(f"{edges[k]:.1f}–{edges[k+1]:.1f}: {hist[k]/hist.sum():.1%}"
                                               for k in range(len(hist))))
    top = sorted(trans.items(), key=lambda kv: -kv[1])[:8]
    print("   8 chuyển dời mạnh nhất (|ℓ|, n trội; khe):")
    for (i, j), wgt in top:
        print(f"     ({lab_l[j]:.2f}, {lab_n[j]:.2f}) → ({lab_l[i]:.2f}, {lab_n[i]:.2f}): khe {E[j]-E[i]:.3f}, "
              f"trọng số {wgt/tot:.1%}")


if __name__ == "__main__":
    main()
