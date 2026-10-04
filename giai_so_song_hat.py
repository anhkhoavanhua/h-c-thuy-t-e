"""Kiểm khẳng định "Quá khứ = sóng, Hiện tại = hạt, ngưỡng = đồng thời" cho E_dual (HocThuyetE (7).pdf, Bảng I).

(A) Cổ điển (Langevin Bảng I, k_BT = 1 MeV): năng lượng mỗi va chạm Carnot ΔE = ½ m_eff u² ⟹ bước sóng de Broglie
    của chuyển động tương đối λ = h/(m_eff u) = 2πℓ0/√(2 m_eff ΔE/ħω0); so với vùng chuyển động κ = ℓ0/3 và
    bán kính hạt nhân. "Hạt định xứ tại điểm va chạm" đòi λ ≪ κ.
(B) Lượng tử (mô hình MCWF V.C.1, giai_so_mcwf_tach_buoc.model): chia thời gian thành
      Quá khứ  = đã > T0 kể từ bước làm nguội gần nhất,
      Hiện tại = trong T0 sau một bước làm nguội;
    đo trên trạng thái điều kiện: độ xác định pha rô-to |⟨e^{iθ}⟩|, Δℓ, Δx của hạt (so với κ), và tỉ số tham gia
    trong cơ sở riêng PR = 1/Σ|c_j|⁴ (PR ≈ 1–2: sóng dừng/trạng thái dừng; PR lớn: bó sóng chuyển động).
(C) Mật độ mức: số mức của mô hình 2 bậc tự do dưới E so với số trạng thái exciton Ericson của hạt nhân A = 70.
Chạy: python3 giai_so_song_hat.py [A|B|C]   (A ≈ 40 s, B ≈ 30 s, 8 lõi)
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
L0_FM, HW = 2.0364, 10.0


def _lang(a):
    from giai_so_dong_ho_photon import _traj
    return _traj(a)[1]


def part_a(N=8):
    from giai_so_burst_g import cfg_for
    cfg = cfg_for(0.1)
    with Pool() as pool:
        evs = pool.map(_lang, [(cfg, 0.1, 91000 + i, 60 * T0, 0.05) for i in range(N)])
    dE = np.array([e[3] for ev in evs for e in ev if e[0] == 'khóa' and e[1] > 5 * T0 and e[3] > 1e-6])
    print(f"=== A. Bước sóng de Broglie của chuyển động tương đối ở va chạm Carnot ({len(dE)} va chạm, 1 MeV) ===")
    print(f"ΔE_Carnot: TB {dE.mean():.3f} ħω0 ({dE.mean()*HW:.2f} MeV), trung vị {np.median(dE):.3f}, "
          f"10–90% {np.percentile(dE,10):.3f}–{np.percentile(dE,90):.3f}")
    kap, Rn = 1 / 3, 1.2 * 70**(1 / 3) / L0_FM
    print(f"{'m_eff':>6} {'λ trung vị (fm)':>16} {'λ 10% (fm)':>11} {'λ/κ trung vị':>13} {'λ/R_hạt nhân':>13} "
          f"{'va chạm có λ < κ':>17}")
    for me in (0.9, 1.0):                                 # m_eff ∈ [mI/(I+mκ²), m] = [0.9, 1.0]
        lam = 2 * math.pi / np.sqrt(2 * me * dE)          # theo ℓ0
        print(f"{me:6.2f} {np.median(lam)*L0_FM:16.2f} {np.percentile(lam,10)*L0_FM:11.2f} "
              f"{np.median(lam)/kap:13.1f} {np.median(lam)/Rn:13.2f} {np.mean(lam < kap):17.1%}")
    print(f"(κ = {kap*L0_FM:.3f} fm, R_hạt nhân = {Rn*L0_FM:.2f} fm; λ < κ cần ΔE > {(2*math.pi/kap)**2/2/0.9:.0f} ħω0)\n")


def _mcwf(args):
    from giai_so_mcwf_tach_buoc import model, NMAX
    seed, Tmax, dt = args
    m = model()
    H, Ls, V, E = m["H"], m["Ls"], m["V"], m["E"]
    d = m["d"]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    Ip = np.eye(NMAX + 1)
    up = np.diag(np.ones(d - 1), -1).astype(complex)
    Eith = np.kron(up, Ip)                               # e^{iθ}
    Lop = np.kron(np.diag(m["ls"].astype(float)), Ip)
    a = np.diag(np.sqrt(np.arange(1, NMAX + 1)), 1)
    X = np.kron(np.eye(d), (a + a.T) / math.sqrt(2))
    rng = random.Random(seed)
    psi = np.zeros(d * (NMAX + 1), complex)
    psi[int(np.where(m["ls"] == 0)[0][0]) * (NMAX + 1)] = 1
    t, r, last = 0.0, rng.random(), -1e9
    rows = {"past": [], "present": []}
    step = 0
    while t < Tmax:
        psi = U @ psi
        t += dt
        step += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            psi = Ls[k] @ psi
            psi /= np.linalg.norm(psi)
            if k == 1:
                last = t
            r = rng.random()
        if t > Tmax / 5 and step % 20 == 0:
            p = psi / math.sqrt(np.vdot(psi, psi).real)
            ev = lambda O: np.vdot(p, O @ p)
            c = V.conj().T @ p
            lm = ev(Lop).real
            xm = ev(X).real
            rows["past" if t - last > T0 else "present"].append(
                (abs(ev(Eith)), math.sqrt(max(0.0, ev(Lop @ Lop).real - lm**2)),
                 math.sqrt(max(0.0, ev(X @ X).real - xm**2)), 1 / np.sum(np.abs(c)**4), abs(lm)))
    return rows


def part_b(N=16, Tmax=1500.0, dt=0.01):
    with Pool() as pool:
        res = pool.map(_mcwf, [(900 + i, Tmax, dt) for i in range(N)])
    print(f"=== B. Trạng thái điều kiện lượng tử (MCWF V.C.1, {N} quỹ đạo × {Tmax:g}) ===")
    print(f"{'pha':>9} {'tỉ lệ t/g':>9} {'|⟨e^iθ⟩| TB':>12} {'(90%)':>7} {'Δℓ TB':>7} {'Δx TB (ℓ0)':>11} {'Δx/κ':>6} "
          f"{'PR TB':>7} {'(trung vị)':>10} {'|⟨ℓ⟩| TB':>9}")
    tot = sum(len(rr[k]) for rr in res for k in rr)
    for key, nm in (("past", "Quá khứ"), ("present", "Hiện tại")):
        A = np.array([row for rr in res for row in rr[key]])
        print(f"{nm:>9} {len(A)/tot:9.1%} {A[:,0].mean():12.3f} {np.percentile(A[:,0],90):7.3f} {A[:,1].mean():7.3f} "
              f"{A[:,2].mean():11.3f} {A[:,2].mean()*3:6.2f} {A[:,3].mean():7.2f} {np.median(A[:,3]):10.2f} "
              f"{A[:,4].mean():9.3f}")
    print("(|⟨e^iθ⟩| = 1: pha rô-to xác định, = 0: trải đều vòng; Δx dao động tử cơ bản = 0.707 ℓ0; PR = 1: một trạng thái dừng)\n")


def part_c():
    from giai_so_u_hermite_drude import basis
    E = np.linalg.eigvalsh(basis(lmax=12, nmax=12)[4])
    E = E - E[0]
    g = 70 / 13.0
    print("=== C. Mật độ mức: mô hình 2 bậc tự do vs exciton Ericson (A = 70, g = A/13 MeV⁻¹) ===")
    print(f"{'E (MeV)':>8} {'số mức mô hình ≤ E':>19} {'ρ Ericson tổng (/MeV)':>22} {'khoảng cách mức mô hình (MeV)':>30}")
    for Em in (5, 10, 20, 30, 40):
        nmod = int(np.sum(E * HW <= Em))
        rho = sum(g**n * Em**(n - 1) / (math.factorial(p) * math.factorial(n - p) * math.factorial(n - 1))
                  for n in range(1, 16) for p in range(0, n + 1))
        sel = np.sort(E[E * HW <= Em + 5]) * HW
        sp = np.mean(np.diff(sel[sel > Em - 5])) if np.sum(sel > Em - 5) > 1 else float('nan')
        print(f"{Em:8d} {nmod:19d} {rho:22.3g} {sp:30.2f}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
