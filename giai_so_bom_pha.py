"""Giải số hai đề xuất cho phần lượng tử:
(A) Tần số cắt bể Drude ħω_c ≳ 50–100 MeV: tỉ lệ ma sát giữ lại tại ω0 và Ω = 2.63ω0 theo ω_c; ω_c tối thiểu để
    giữ ≥ 90%; cân bằng chính xác ⟨x²⟩, ⟨p²⟩; số học "m_σ/4" và ngưỡng 2m_π.
(B) Bơm pha không kết hợp + làm nguội Davies (quỹ đạo lượng tử MCWF):
    - dạng nguyên văn: L_p = √Γ0 e^{iθ} (tốc độ vô hướng — không bão hòa)
    - dạng sửa: L_p = √Γ0 e^{iθ} f(ℓ̂), f(ℓ) = √max(0, 1 − ℓ²/ℓ0²), ℓ0 = 4
    - làm nguội: L_c = √Γc û₋ (Davies T=0 của vận tốc tương đối)
    Đo: ⟨ℓ⟩ dừng, năng lượng mỗi bước làm nguội, số bước làm nguội mỗi đợt (khoảng cách < T0),
    thời gian chờ giữa các đợt (CV), so với cổ điển (hình học TB 5.1; mũ; 0.40 ħω0).
Đơn vị ħ = m = ω0 = 1 (ħω0 = 10 MeV), I = 1, κ = 1/3, ΔE = 1, g = 0.1.
Chạy: python3 giai_so_bom_pha.py [A|B]   (B ≈ 3–5 phút, 8 lõi)
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

from giai_so_u_hermite_drude import basis

T0 = 2 * math.pi


def qbm_moments(wc, kT=0.1, gam=0.7):
    w = np.linspace(1e-6, 200.0, 2000001)
    gh = gam * wc / (wc - 1j * w)
    chi = 1 / (1 - w**2 - 1j * w * gh)
    cth = 1 / np.tanh(w / (2 * kT))
    return (np.trapezoid(cth * chi.imag, w) / math.pi, np.trapezoid(w**2 * cth * chi.imag, w) / math.pi)


def part_a():
    print("=== A. Tần số cắt bể Drude ===")
    print(f"{'ħω_c (MeV)':>11} {'giữ tại ω0':>11} {'giữ tại Ω=2.63':>15}")
    for E in (4.5, 25, 50, 79, 100, 280):
        wc = E / 10
        print(f"{E:11g} {wc*wc/(wc*wc+1):11.1%} {wc*wc/(wc*wc+2.63**2):15.1%}")
    wc90 = 2.63 * 3
    print(f"Giữ ≥ 90% tại Ω cần ω_c ≥ 3Ω = {wc90:.2f} ⟹ ħω_c ≥ {wc90*10:.0f} MeV")
    for E in (50, 100, 280):
        x2, p2 = qbm_moments(E / 10)
        print(f"cân bằng chính xác (k_BT = 1 MeV), ħω_c = {E} MeV: ⟨x²⟩ = {x2:.4f}, ⟨p²⟩ = {p2:.4f}, tích = {x2*p2:.4f}")
    print(f"m_σ c² = 400–500 MeV ⟹ m_σ/4 = {400/4:.0f}–{500/4:.0f} MeV (đề xuất ghi 50–100); ngưỡng phổ liên tục 2π: "
          f"2m_πc² = {2*139.57:.0f} MeV\n")


def _traj(args):
    Gp, Gc, l0, mode, seed, Tmax, dt = args
    ls, cos_t, sin_t, pth, H, u = basis()
    nmax = 6
    E, V = np.linalg.eigh(H)
    Lc = math.sqrt(Gc) * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    d = len(ls)
    up = np.diag(np.ones(d - 1), -1).astype(complex)              # e^{iθ}: ℓ → ℓ+1
    if mode == 'literal':
        F = np.eye(d)
    else:
        F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / l0**2)))
    Lp = math.sqrt(Gp) * np.kron(up @ F, np.eye(nmax + 1))
    Ls = [Lp, Lc]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    Lsq = np.kron(np.diag(ls.astype(float)), np.eye(nmax + 1))
    rng = random.Random(seed)
    psi = np.zeros(d * (nmax + 1), complex)
    psi[int(np.where(ls == 0)[0][0]) * (nmax + 1)] = 1
    t, r = 0.0, rng.random()
    jumps, lsum, lcount = [], 0.0, 0
    while t < Tmax:
        psi = U @ psi
        t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            Eb = float(np.real(np.vdot(psi, H @ psi)))
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            psi = Ls[k] @ psi
            psi /= np.linalg.norm(psi)
            Ea = float(np.real(np.vdot(psi, H @ psi)))
            jumps.append((t, k, Eb - Ea))
            r = rng.random()
        if t > Tmax / 5 and int(t / dt) % 50 == 0:
            p = psi / math.sqrt(np.vdot(psi, psi).real)
            lsum += float(np.real(np.vdot(p, Lsq @ p))); lcount += 1
    return jumps, lsum / max(lcount, 1)


def part_b(N=16, Tmax=1500.0, dt=0.01):
    print(f"=== B. Bơm pha + làm nguội Davies (MCWF, {N} quỹ đạo × {Tmax:g} đv t ≈ {Tmax/T0:.0f} T0) ===")
    configs = [(0.05, 0.5, 4, 'literal'), (0.05, 0.5, 4, 'sat'), (0.2, 0.5, 4, 'sat')]
    with Pool() as pool:
        for Gp, Gc, l0, mode in configs:
            res = pool.map(_traj, [(Gp, Gc, l0, mode, 900 + i, Tmax, dt) for i in range(N)])
            cool_E, bursts, gaps, npump, ncool, lm = [], [], [], 0, 0, []
            for jumps, lbar in res:
                lm.append(lbar)
                jj = [j for j in jumps if j[0] > Tmax / 5]
                npump += sum(1 for j in jj if j[1] == 0)
                cj = [j for j in jj if j[1] == 1]
                ncool += len(cj)
                cool_E += [j[2] for j in cj]
                cur, last = 0, None
                for j in cj:
                    if last is None or j[0] - last < T0:
                        cur += 1
                    else:
                        bursts.append(cur); gaps.append(j[0] - last); cur = 1
                    last = j[0]
            span = N * (Tmax * 0.8) / T0
            name = {"literal": "nguyên văn (không bão hòa)", "sat": "bão hòa f(ℓ̂)"}[mode]
            print(f"-- {name}, Γ0 = {Gp}, Γc = {Gc}, ℓ0 = {l0}: ⟨ℓ⟩ dừng TB {np.mean(lm):.2f} (min {min(lm):.2f}, max {max(lm):.2f}); "
                  f"nhảy bơm {npump/span:.2f}/T0, nhảy làm nguội {ncool/span:.2f}/T0")
            if cool_E:
                cE = np.array(cool_E)
                print(f"   năng lượng mỗi bước làm nguội: TB {cE.mean():.3f} ħω0, trung vị {np.median(cE):.3f}, 10–90% "
                      f"{np.percentile(cE,10):.3f}–{np.percentile(cE,90):.3f}")
            if bursts:
                b = np.array(bursts); g = np.array(gaps) / T0
                cnt = {k: int(np.sum(b == k)) for k in range(1, min(int(b.max()), 8) + 1)}
                print(f"   đợt làm nguội (khoảng cách < T0): {len(b)} đợt, bước/đợt TB {b.mean():.2f}, CV {b.std()/b.mean():.2f} "
                      f"(hình học cùng TB: CV {math.sqrt(1-1/b.mean()):.2f}); phân bố " + ", ".join(f"{k}:{v}" for k, v in cnt.items()))
                print(f"   thời gian chờ giữa đợt: TB {g.mean():.2f} T0, CV {g.std()/g.mean():.2f} (mũ: 1)", flush=True)
    print("tham chiếu cổ điển: bước/đợt hình học TB 5.12 (CV 0.90); chờ mũ TB ≈ 4.6 T0 (CV ≈ 1); xung TB 0.40 ħω0\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
