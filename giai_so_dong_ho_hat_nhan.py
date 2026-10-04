"""Đồng hồ hạt nhân ²²⁹Th làm đồng hồ Page–Wootters trong vũ trụ của học thuyết E_dual (HocThuyetE (7).pdf).

Đồng hồ hạt nhân thật: chuyển dời đồng phân ²²⁹ᵐTh, E_c = 8.355733 eV (Tiedau et al., PRL 132, 182501, 2024;
ν = 2 020 407 384 335 kHz, Zhang et al., Nature 633, 63, 2024). Động cơ: T0 = 2π/ω0 = 4.136·10⁻²² s (ħω0 = 10 MeV).

(N1) Đòi hỏi: từ mô phỏng Langevin Part I (Bảng I, 1 MeV), đọc đồng hồ có sai số Gauss σ_t ⟹ xác suất xếp đúng
     Quá khứ (tái khóa trước) < Hiện tại (gãy) < Tương lai (tái khóa sau) cho từng đợt; tìm σ_t cho 90%, 99%.
(N2) Giới hạn lượng tử Mandelstam–Tamm σ_t ≥ ħ/(2ΔE): trải năng lượng cần có của đồng hồ; số hạt nhân Th cần
     (N độc lập: ΔE = √N E_c/2; N vướng víu GHZ: ΔE = N E_c/2); so với một hạt nhân Th.
(N3) Page–Wootters với đồng hồ thang M mức (nhân Fejér, độ phân giải Δ): trạng thái điều kiện của động cơ lượng tử
     (MCWF mô hình V.C.1, một bản ghi môi trường) — độ thấy được của đợt V(Δ) = std(tín hiệu làm trơn)/std(gốc),
     và khoảng cách tới trạng thái dừng; ở độ phân giải của Th: vũ trụ trông tĩnh.
Chạy: python3 giai_so_dong_ho_hat_nhan.py   (≈ 3 phút, 8 lõi)
"""
import math
import random
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
T0_S = 4.1357e-22
HBAR_EVS = 6.582119569e-16
E_TH = 8.355733                       # eV


def _lang(seed):
    from giai_so_dong_ho_photon import _traj
    from giai_so_burst_g import cfg_for
    return _traj((cfg_for(0.1), 0.1, seed, 80 * T0, 0.05))[1]


def bursts_from(evs):
    out = []
    for ev in evs:
        for i, e in enumerate(ev):
            if e[0] != "gãy" or e[1] < 5 * T0:
                continue
            prev = next((x for x in reversed(ev[:i]) if x[0] == "khóa"), None)
            if prev is None or e[1] - prev[1] < T0:
                continue
            for j in range(i + 1, len(ev)):
                if ev[j][0] == "khóa":
                    nxt = next((y for y in ev[j + 1:] if y[0] == "gãy"), None)
                    if nxt is None or nxt[1] - ev[j][1] >= T0:
                        out.append((prev[1], e[1], ev[j][1]))
                        break
    return np.array(out) / T0                                  # đơn vị T0


def part_n1_n2():
    with Pool() as pool:
        evs = pool.map(_lang, range(21000, 21016))
    B = bursts_from(evs)
    dwell, dur = B[:, 1] - B[:, 0], B[:, 2] - B[:, 1]
    rng = np.random.default_rng(3)
    print(f"=== N1. Đòi hỏi độ phân giải đồng hồ (Langevin Bảng I, 1 MeV; {len(B)} đợt) ===")
    print(f"lưu Quá khứ: trung vị {np.median(dwell):.2f} T0, 10% {np.percentile(dwell,10):.2f}; đợt Hiện tại: trung vị "
          f"{np.median(dur):.2f} T0, 10% {np.percentile(dur,10):.3f}")
    print(f"{'σ_t (T0)':>9} {'σ_t (s)':>10} {'xếp đúng 3 pha':>15} {'Hiện tại sau Quá khứ':>21} {'Tương lai sau Hiện tại':>23}")
    sig = np.logspace(-3, 3, 25)
    good = []
    for s in sig:
        noise = rng.normal(0, s, size=(len(B), 400, 3))
        t = B[:, None, :] + noise
        ok3 = np.mean((t[:, :, 0] < t[:, :, 1]) & (t[:, :, 1] < t[:, :, 2]))
        ok1 = np.mean(t[:, :, 0] < t[:, :, 1]); ok2 = np.mean(t[:, :, 1] < t[:, :, 2])
        good.append(ok3)
        if s in sig[::3]:
            print(f"{s:9.3f} {s*T0_S:10.2e} {ok3:15.1%} {ok1:21.1%} {ok2:23.1%}")
    good = np.array(good)
    s90 = float(np.exp(np.interp(-0.90, -good, np.log(sig))))
    s99 = float(np.exp(np.interp(-0.99, -good, np.log(sig))))
    print(f"(đoán ngẫu nhiên thứ tự 3 sự kiện = 1/6 = 16.7%)")
    print(f"σ_t cho 90% xếp đúng: {s90:.3f} T0 = {s90*T0_S:.2e} s; cho 99%: {s99:.4f} T0 = {s99*T0_S:.2e} s\n")

    print("=== N2. Giới hạn Mandelstam–Tamm σ_t ≥ ħ/(2ΔE) và đồng hồ hạt nhân ²²⁹Th ===")
    tick = HBAR_EVS * 2 * math.pi / E_TH
    s_th = HBAR_EVS / E_TH                                      # một hạt nhân, ΔE = E_c/2
    print(f"²²⁹ᵐTh: E_c = {E_TH} eV, chu kỳ dao động 1/ν = {tick:.3e} s = {tick/T0_S:.3e} T0; "
          f"σ_t tối thiểu của MỘT hạt nhân = ħ/E_c = {s_th:.2e} s = {s_th/T0_S:.2e} T0")
    for lab, s in (("90%", s90), ("99%", s99)):
        dE = HBAR_EVS / (2 * s * T0_S)                          # eV
        n_sql = (2 * dE / E_TH)**2; n_ghz = 2 * dE / E_TH
        print(f"{lab}: cần ΔE ≥ {dE/1e6:.1f} MeV ⟹ N hạt nhân Th độc lập ≥ {n_sql:.2e}; vướng víu GHZ ≥ {n_ghz:.2e}")
    print(f"thang năng lượng cần ≈ {HBAR_EVS/(2*s90*T0_S)/1e6:.0f}–{HBAR_EVS/(2*s99*T0_S)/1e6:.0f} MeV ≈ "
          f"{HBAR_EVS/(2*s90*T0_S)/1e7:.1f}–{HBAR_EVS/(2*s99*T0_S)/1e7:.1f} ħω0 — thang của chính dao động hạt nhân "
          "(chuyển dời gamma MeV), lớn hơn chuyển dời đồng phân Th (eV) khoảng 10⁶ lần\n")
    return B, s90


def mcwf_series(seed, Tmax=600.0, dt=0.01):
    from giai_so_part2_master import model
    m = model(8, 6, False)
    H, Ls = m["H"], list(m["Ls"].values())
    names = list(m["Ls"])
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(R)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = np.zeros(H.shape[0], complex); psi[int(np.where(m["ls"] == 0)[0][0]) * 7] = 1
    t, r, k = 0.0, rng.random(), 0
    ts, lv, jc = [], [], []
    while t < Tmax:
        psi = U @ psi; t += dt; k += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            j = int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random())); j = min(j, len(Ls) - 1)
            psi = Ls[j] @ psi; psi /= np.linalg.norm(psi)
            if names[j] == "Carnot":
                jc.append(t)
            r = rng.random()
        if k % 5 == 0:
            p = psi / np.linalg.norm(psi)
            ts.append(t); lv.append(np.real(np.vdot(p, Labs @ p)))
    return np.array(ts), np.array(lv), np.array(jc)


def fejer(tau, Delta, M=40):
    """|K(τ)|² của đồng hồ thang M mức, độ phân giải Δ (không điểm đầu tiên): đọc 'thời gian' bằng phép chiếu lên |t⟩."""
    eps = 2 * math.pi / (M * Delta)
    x = eps * tau / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        k = np.where(np.abs(np.sin(x)) < 1e-12, 1.0, (np.sin(M * x) / (M * np.sin(x)))**2)
    return k


def part_n3():
    print("=== N3. Page–Wootters với đồng hồ thang M mức: động cơ lượng tử V.C.1 nhìn qua đồng hồ độ phân giải Δ ===")
    with Pool() as pool:
        res = pool.map(mcwf_series, range(4))
    print(f"{'Δ (T0)':>8} {'Δ (s)':>10} {'thấy được đợt V(Δ)':>19} {'tương quan với gốc':>19} {'đồng hồ cần (MeV)':>18}")
    for Dl in (0.02, 0.1, 0.5, 2.0, 10.0, 50.0):
        Vs, Cs = [], []
        for ts, lv, jc in res:
            m = ts > 100
            t, x = ts[m], lv[m]
            dtt = t[1] - t[0]
            half = min(int(3 * Dl * T0 / dtt), len(t) // 2 - 1)
            tau = np.arange(-half, half + 1) * dtt
            ker = fejer(tau, Dl * T0); ker /= ker.sum()
            sm = np.convolve(x - x.mean(), ker, mode="same")
            core = slice(half, len(x) - half) if half > 0 else slice(None)
            Vs.append(np.std(sm[core]) / np.std(x[core] - x.mean()))
            Cs.append(np.corrcoef(sm[core], x[core])[0, 1])
        dE = HBAR_EVS / (2 * Dl * T0_S) / 1e6
        print(f"{Dl:8.2f} {Dl*T0_S:10.2e} {np.mean(Vs):19.3f} {np.mean(Cs):19.3f} {dE:18.2f}")
    s_th = HBAR_EVS / E_TH / T0_S
    print(f"(đồng hồ ²²⁹Th một hạt nhân: Δ ≈ {s_th:.1e} T0 ≫ mọi dòng trên ⟹ V → 0: người quan sát dùng đồng hồ Th chỉ thấy"
          " trạng thái dừng — vũ trụ tĩnh)\n")


if __name__ == "__main__":
    part_n1_n2()
    part_n3()
