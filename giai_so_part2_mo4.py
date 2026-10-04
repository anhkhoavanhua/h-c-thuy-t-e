"""Ba câu hỏi còn lại của Part II (E_dual, HocThuyetE (7).pdf).

R1  Nhiễu nhiệt cho rô-to lượng tử: mô hình chặn gãy (giếng R = 3, Carnot chỉ trên trạng thái không liên kết p* = 0.3) với
    trần bơm ℓ0 = 3.2 (bơm một mình dừng ở |ℓ| = 3, dưới ngưỡng gãy) + bể nhiệt rô-to Davies ghép qua cosθ:
    xuống Γt|X_ij|², lên Γt|X_ij|² e^{−ω/k_BT}. Quét k_BT → thời gian lưu (TB, CV), va chạm/đợt; so Part I
    (4.5–6 T0, CV 1.0–1.1, 5.12). Kiểm định Arrhenius ln(lưu) ∝ ΔE/k_BT.
R1b Đối chứng: bể T ≈ 0 và quét Γt (tách nhiễu lượng tử khỏi nhiễu nhiệt).
Chạy: python3 giai_so_part2_mo4.py [R1|R1B]
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

import giai_so_part2_mo as P
from giai_so_part2_mo3 import burst_stats

T0 = 2 * math.pi


def r1_model(kT, Gt=0.5, G0=2.0, l0=3.2, Gc=0.5, pstar=0.3):
    m, H, Pout = P.a1_model()
    E, V = np.linalg.eigh(H)
    po = np.array([float(np.real(V[:, j].conj() @ Pout @ V[:, j])) for j in range(len(E))])
    ue = V.conj().T @ m["u"] @ V
    Lc = V @ (math.sqrt(Gc) * np.triu(ue, 1) * (po > pstar)[None, :]) @ V.conj().T
    Ls = [Lc, P.pump(m, G0, l0)]
    if kT > 0:
        X = V.conj().T @ np.kron(m["cos"], np.eye(m["nmax"] + 1)) @ V
        w = E[:, None] - E[None, :]
        Ls.append(V @ (math.sqrt(Gt) * np.triu(X, 1)) @ V.conj().T)                       # xuống
        Ls.append(V @ (math.sqrt(Gt) * np.tril(X, -1) * np.exp(-np.clip(w, 0, None) / (2 * kT))) @ V.conj().T)  # lên
    return m, H, Ls, V


def _traj(args):
    kT, seed, Tmax, dt = args[:4]
    Gt = args[4] if len(args) > 4 else 0.5
    m, H, Ls, V = r1_model(kT, Gt=Gt)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = V[:, 0].astype(complex)
    t, r, jumps = 0.0, rng.random(), []
    while t < Tmax:
        psi = U @ psi; t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            k = min(int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random())), len(Ls) - 1)
            if k == 0:
                Eb = np.real(np.vdot(psi, H @ psi)); lb = np.real(np.vdot(psi, Labs @ psi))
            psi = Ls[k] @ psi; psi /= np.linalg.norm(psi)
            if k == 0:
                jumps.append((t, Eb - np.real(np.vdot(psi, H @ psi)), lb))
            r = rng.random()
    return jumps


def part_r1(N=16, Tmax=2400.0, dt=0.01):
    m, H, Ls, V = r1_model(0.0)
    E, V = np.linalg.eigh(H)
    L2 = m["Lop"] @ m["Lop"]
    l2 = np.array([np.real(V[:, j].conj() @ L2 @ V[:, j]) for j in range(len(E))])
    e3 = E[np.argmin(np.abs(np.sqrt(l2) - 3) + 10 * (E > 10))]
    e4 = E[np.argmin(np.abs(np.sqrt(l2) - 4) + 10 * (E > 15))]
    print("=== R1. Nhiễu nhiệt cho rô-to (trần bơm ℓ0 = 3.2, Γ0 = 2, bể nhiệt rô-to Γt = 0.5 qua cosθ) ===")
    print(f"khe rô-to |ℓ| 3 → 4 trong H + V_b: ΔE ≈ {e4 - e3:.2f} ħω0 = {10*(e4-e3):.0f} MeV ⟹ Arrhenius ở 1 MeV: "
          f"e^(−ΔE/kT) = {math.exp(-(e4-e3)/0.1):.1e}")
    print(f"{'k_BT':>8} {'đợt':>5} {'va chạm/đợt':>12} {'lưu TB (T0)':>12} {'CV':>5} {'ΔE/va chạm':>11} {'⟨|ℓ|⟩ khởi phát':>16}")
    rows = []
    with Pool() as pool:
        for kT in (0.0, 0.1, 0.3, 0.5, 1.0):
            st = burst_stats(pool.map(_traj, [(kT, 11000 + i, Tmax, dt) for i in range(N)]), Tmax / 6)
            rows.append((kT, st))
            print(f"{kT*10:6.0f}MeV {st['n']:5d} {st['size']:12.2f} {st['dwell']:12.2f} {st['cv']:5.2f} {st['eJ']:11.3f} "
                  f"{st['lon']:16.2f}", flush=True)
    ok = [(kT, st) for kT, st in rows if kT > 0 and st["n"] > 10 and np.isfinite(st["dwell"])]
    if len(ok) >= 2:
        x = np.array([1 / kT for kT, _ in ok]); y = np.log([st["dwell"] for _, st in ok])
        sl = np.polyfit(x, y, 1)[0]
        print(f"Arrhenius: ln(lưu) theo 1/k_BT có hệ số góc {sl:.2f} ħω0 (rào hiệu dụng {10*sl:.1f} MeV)")
    print("(Part I: lưu 4.5–6 T0, CV 1.0–1.1, 5.12 va chạm/đợt, rào hiệu dụng ≈ 0.4 MeV)\n")


def part_r1b(N=16, Tmax=2400.0, dt=0.01):
    """Đối chứng: bể T ≈ 0 (chỉ nhảy xuống) và quét cường độ ghép Γt — tách nhiễu lượng tử khỏi nhiễu nhiệt."""
    print("=== R1b. Bể rô-to ở T ≈ 0 và theo cường độ ghép Γt ===")
    print(f"{'k_BT':>9} {'Γt':>5} {'đợt':>5} {'va chạm/đợt':>12} {'lưu TB (T0)':>12} {'CV':>5} {'ΔE/va chạm':>11}")
    with Pool() as pool:
        for kT, Gt in ((1e-4, 0.5), (0.1, 0.1), (0.1, 0.5), (0.1, 2.0)):
            st = burst_stats(pool.map(_traj, [(kT, 12000 + i, Tmax, dt, Gt) for i in range(N)]), Tmax / 6)
            print(f"{kT*10:7.3f}MeV {Gt:5.1f} {st['n']:5d} {st['size']:12.2f} {st['dwell']:12.2f} {st['cv']:5.2f} {st['eJ']:11.3f}",
                  flush=True)
    print()


if __name__ == "__main__":
    for a in sys.argv[1:] or ["R1", "R1B"]:
        {"R1": part_r1, "R1B": part_r1b}[a.upper()]()
