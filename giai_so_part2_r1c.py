"""Hiệu chỉnh kết hợp (Part II, R1c): trần bơm ℓ0 + bể rô-to Γt — so thống kê đợt với Part I.
Chạy: python3 giai_so_part2_r1c.py   (≈ 10 phút, 8 lõi)
"""
import sys, math, random, numpy as np
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from multiprocessing import Pool
import giai_so_part2_mo4 as R
from giai_so_part2_mo3 import burst_stats
T0 = 2*math.pi
def traj(args):
    l0, G0, Gt, seed = args
    m, H, Ls, V = R.r1_model(0.1, Gt=Gt, G0=G0, l0=l0)
    dt, Tmax = 0.01, 2400.0
    Heff = H - 0.5j*sum(L.conj().T@L for L in Ls); w, Rm = np.linalg.eig(Heff)
    U = Rm@np.diag(np.exp(-1j*w*dt))@np.linalg.inv(Rm); Labs = np.abs(m["Lop"])
    rng = random.Random(seed); psi = V[:,0].astype(complex); t, r, jumps = 0.0, rng.random(), []
    while t < Tmax:
        psi = U@psi; t += dt; nrm = np.vdot(psi,psi).real
        if nrm < r:
            psi /= math.sqrt(nrm); ws = np.array([np.linalg.norm(L@psi)**2 for L in Ls])
            k = min(int(np.searchsorted(np.cumsum(ws/ws.sum()), rng.random())), len(Ls)-1)
            if k == 0: Eb = np.real(np.vdot(psi,H@psi)); lb = np.real(np.vdot(psi,Labs@psi))
            psi = Ls[k]@psi; psi /= np.linalg.norm(psi)
            if k == 0: jumps.append((t, Eb-np.real(np.vdot(psi,H@psi)), lb))
            r = rng.random()
    return jumps
if __name__ == "__main__":
    print(f"{'ℓ0':>4} {'Γ0':>4} {'Γt':>4} {'đợt':>5} {'va chạm/đợt':>12} {'lưu (T0)':>9} {'CV':>5} {'ΔE/va chạm':>11} {'ℓ khởi phát':>12}")
    with Pool() as pool:
        for l0, G0, Gt in ((4.0, 2.0, 0.5), (4.0, 2.0, 1.0), (4.0, 1.0, 1.0), (3.6, 2.0, 1.0)):
            st = burst_stats(pool.map(traj, [(l0, G0, Gt, 13000+i) for i in range(16)]), 400.0)
            print(f"{l0:4.1f} {G0:4.1f} {Gt:4.1f} {st['n']:5d} {st['size']:12.2f} {st['dwell']:9.2f} {st['cv']:5.2f} {st['eJ']:11.3f} {st['lon']:12.2f}", flush=True)
    print("(Part I: 5.12 va chạm/đợt; lưu 4.5–6 T0; CV 1.0–1.1; 0.40–0.46 ħω0; ngưỡng 3.56)")
