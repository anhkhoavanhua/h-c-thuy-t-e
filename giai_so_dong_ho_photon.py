"""Đồng hồ photon Page–Wootters (Moreva et al., Phys. Rev. A 89, 052122, 2014) trên học thuyết E_dual
(bản thảo HocThuyetE (5).pdf: τ = 10, γx = 0.70, g = 0.1, f = 1.05, k_BT = 1 MeV).

(A) Tái lập thí nghiệm lý tưởng: photon đồng hồ C và photon hệ S (phân cực H/V), trạng thái singlet,
    Ĥ = ω(σ_y^C + σ_y^S) (bản làm trễ quay phân cực). Kiểm:
      – "người quan sát ngoài": Ĥ|Ψ⟩⟩ = 0, trạng thái toàn cục không đổi khi quay;
      – "người quan sát trong": đọc C ở góc φ ⟹ trạng thái điều kiện của S = e^{−iωσ_y t}|ψ0⟩ với t = φ/ω;
      – độ phân giải của đồng hồ qubit: ⟨φ|φ'⟩ = cos(φ − φ'), chu kỳ π.
(B) Động cơ E_dual làm đồng hồ: số đọc = pha tích lũy Φ(t) của rô-to (như góc phân cực của photon đồng hồ).
    Mô phỏng Langevin với tham số Bảng I; với mỗi đợt (Hiện tại) và lần tái khóa kết thúc đợt (Tương lai), đọc Φ tại
    các biến cố. Đo: tỉ lệ đợt mà đồng hồ đọc "Tương lai trước Hiện tại" (Φ_khóa ≤ Φ_gãy); tỉ lệ đợt kéo dài ≥ nửa vòng
    (π, chu kỳ của đồng hồ phân cực) và ≥ một vòng (2π); sai số khi suy thời gian từ số đọc t̂ = Φ/⟨Φ̇⟩.
Chạy: python3 giai_so_dong_ho_photon.py [A|B]   (B ≈ 3 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

import giai_so_dong_co_langevin as L
from giai_so_burst_g import cfg_for
from giai_so_kramers_A import setup

T0 = 2 * math.pi


def part_a():
    sy = np.array([[0, -1j], [1j, 0]])
    I2 = np.eye(2)
    H = np.kron(sy, I2) + np.kron(I2, sy)                    # ω = 1
    H_, V_ = np.array([1, 0]), np.array([0, 1])
    Psi = (np.kron(H_, V_) - np.kron(V_, H_)) / math.sqrt(2)  # singlet
    print("=== A. Thí nghiệm Moreva lý tưởng (qubit đồng hồ + qubit hệ) ===")
    print(f"‖Ĥ|Ψ⟩⟩‖ = {np.linalg.norm(H @ Psi):.2e}  ⟹ trạng thái toàn cục là nghiệm Wheeler–DeWitt")
    w, U = np.linalg.eigh(H)
    maxdev = 0.0
    for t in np.linspace(0, 3, 13):
        Ut = U @ np.diag(np.exp(-1j * w * t)) @ U.conj().T
        maxdev = max(maxdev, 1 - abs(np.vdot(Psi, Ut @ Psi))**2)
    print(f"người quan sát ngoài: max(1 − |⟨Ψ|e^(−iĤt)|Ψ⟩|²) trên t ∈ [0,3] = {maxdev:.2e}  (tĩnh)")
    psi0 = None
    fids = []
    for t in np.linspace(0, math.pi, 25):
        phi_state = np.array([math.cos(t), math.sin(t)])      # C đo ở phân cực góc φ = ωt
        cond = np.array([np.vdot(np.kron(phi_state, e), Psi) for e in (H_, V_)])
        cond /= np.linalg.norm(cond)
        if psi0 is None:
            psi0 = cond.copy()
        ref = (np.cos(t) * I2 - 1j * np.sin(t) * sy) @ psi0    # e^{−iσ_y t}|ψ0⟩
        fids.append(abs(np.vdot(ref, cond))**2)
    print(f"người quan sát trong: độ trùng khớp với Schrödinger trên φ ∈ [0, π]: min {min(fids):.6f}")
    print("độ phân giải: ⟨φ|φ'⟩ = cos(φ−φ'); chỉ φ và φ+π/2 trực giao ⟹ đồng hồ qubit phân biệt tối đa 2 thời điểm;"
          " φ và φ+π cho cùng trạng thái (chu kỳ π)\n")


def _traj(a):
    cfg, kT, seed, T, every = a
    setup(cfg)
    wc = cfg[2]
    r = L.run(0.1, 1.05 * wc, kT, T, 2e-3, seed, gx=0.7, w_start=0.82 * wc, sample_every=every)
    return np.array(r['samples']), r['events']


def phi_at(samples, every, t):
    i = min(int(t / every), len(samples) - 2)
    f = t / every - i
    return samples[i] * (1 - f) + samples[i + 1] * f


def part_b(N=8):
    cfg = cfg_for(0.1)
    T, warm, every = 70 * T0, 10 * T0, 0.05
    print("=== B. Động cơ E_dual làm đồng hồ photon (tham số Bảng I của bản thảo) ===")
    print(f"{'k_BT':>7} {'đợt':>5} {'Tương lai TRƯỚC Hiện tại':>25} {'đợt ≥ π (nửa vòng)':>19} {'đợt ≥ 2π':>9} "
          f"{'lệch t̂ − t RMS (T0)':>20} {'lệch max (T0)':>14}")
    with Pool() as pool:
        for kT in (0.03, 0.1, 0.5):
            res = pool.map(_traj, [(cfg, kT, 51000 + i, T, every) for i in range(N)])
            nb = rev = half = full = 0
            errs, errmax = [], 0.0
            for samp, ev in res:
                # đợt: từ lần gãy mở đợt tới lần khóa đầu tiên kéo dài ≥ T0
                for i, e in enumerate(ev):
                    if e[0] != 'gãy' or e[1] < warm:
                        continue
                    prev = next((x for x in reversed(ev[:i]) if x[0] == 'khóa'), None)
                    if prev is None or e[1] - prev[1] < T0:
                        continue
                    end = None
                    for j in range(i + 1, len(ev)):
                        if ev[j][0] == 'khóa':
                            nxt = next((y for y in ev[j + 1:] if y[0] == 'gãy'), None)
                            if nxt is None or nxt[1] - ev[j][1] >= T0:
                                end = ev[j][1] if nxt is not None else None
                                break
                    if end is None or end > T - every:
                        continue
                    nb += 1
                    d = phi_at(samp, every, end) - phi_at(samp, every, e[1])
                    rev += d <= 0
                    half += abs(d) >= math.pi
                    full += abs(d) >= 2 * math.pi
                # suy thời gian từ số đọc đồng hồ
                k0 = int(warm / every)
                ph = samp[k0:] - samp[k0]
                tt = np.arange(len(ph)) * every
                rate = ph[-1] / tt[-1]
                if abs(rate) > 1e-9:
                    e_ = (ph / rate - tt) / T0
                    errs.append(float(np.sqrt(np.mean(e_**2))))
                    errmax = max(errmax, float(np.max(np.abs(e_))))
            print(f"{kT*10:5.1f}MeV {nb:5d} {rev:>12d}/{nb:<4d} ({rev/max(nb,1):5.1%}) {half/max(nb,1):19.1%} "
                  f"{full/max(nb,1):9.1%} {np.mean(errs):20.2f} {errmax:14.2f}", flush=True)
    print("(Φ đo liên tục, không lấy modulo; 'Tương lai TRƯỚC Hiện tại' = Φ tại lần tái khóa ≤ Φ tại lần gãy)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
