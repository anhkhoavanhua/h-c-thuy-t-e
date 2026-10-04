"""Giải số đề xuất:
(A) Biên M_det(f) = 1 (đề xuất f = 1.14) — lưới mịn f = 1.10…1.20 không nhiễu.
(B) "f ≥ 1.14: liên kết không tồn tại ổn định dù một chu kỳ" — cho phép gãy, không nhiễu: đếm gãy,
    độ dài các pha khóa.
(C) Tích lũy exciton: E_exc(k) = Σ xung Carnot trong một đợt; P(xung đơn ≥ S_n), P(E_exc(2) ≥ S_n),
    P(tổng đợt ≥ S_n), phổ ε_max = E_exc − S_n; tỉ lệ đợt N ≥ 3 (đề xuất ≈ 50%). S_n = 8 MeV.

Đơn vị m=ω0=ħ=ℓ0=1, ħω0 = 10 MeV, τ=10, γx=0.70, g=0.1.
Chạy: python3 giai_so_exciton.py [A|B|C]   (≈ 8 phút tổng, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import analyse, cfg_for
from giai_so_kramers_A import setup
from giai_so_saturation import _det_margin

TAU, GX, HW0, SN = 10.0, 0.70, 10.0, 8.0
T0 = L.T0


def part_a():
    cfg = cfg_for(0.1)
    fs = [1.10 + 0.01 * i for i in range(11)]
    with Pool() as pool:
        res = pool.map(_det_margin, [(cfg, f) for f in fs])
    print("=== A. Lề ứng suất tất định M_det(f) ===")
    prev = None
    for f, M, *_ in res:
        print(f"f={f:.2f}: M_det = {M:.4f}")
        if prev and prev[1] < 1 <= M:
            fc = prev[0] + (1 - prev[1]) * (f - prev[0]) / (M - prev[1])
            print(f"   ⟹ M_det = 1 tại f ≈ {fc:.4f}")
        prev = (f, M)
    print()


def _det_break(a):
    cfg, f, T, warm = a
    setup(cfg)
    wc = cfg[2]
    ev = L.run(1 / TAU, f * wc, 0.0, T, 2e-3, 0, gx=GX, w_start=0.82 * wc)['events']
    locks, last = [], None
    for e in ev:
        if e[0] == 'khóa':
            last = e[1]
        elif last is not None:
            if e[1] > warm:
                locks.append(e[1] - last)
            last = None
    nb = sum(1 for e in ev if e[0] == 'gãy' and e[1] > warm)
    return f, nb, locks


def part_b():
    cfg = cfg_for(0.1)
    T, warm = 50 * T0, 10 * T0
    print("=== B. Cho phép gãy, KHÔNG nhiễu, 40 T0 sau khởi động ===")
    with Pool() as pool:
        res = pool.map(_det_break, [(cfg, f, T, warm) for f in (1.05, 1.10, 1.15, 1.20, 1.30)])
    for f, nb, locks in res:
        if locks:
            print(f"f={f:.2f}: {nb} lần gãy; pha khóa: {len(locks)}, dài nhất {max(locks)/T0:.3f} T0, "
                  f"TB {sum(locks)/len(locks)/T0:.3f} T0, số pha ≥ T0: {sum(1 for w in locks if w >= T0)}")
        else:
            print(f"f={f:.2f}: {nb} lần gãy (khóa suốt, không gãy)")
    print()


def _run(a):
    cfg, kT, seed, T = a
    setup(cfg)
    wc = cfg[2]
    return L.run(1 / TAU, 1.05 * wc, kT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc)['events']


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))]


def part_c(N=16):
    cfg = cfg_for(0.1)
    T, warm = 70 * T0, 10 * T0
    print(f"=== C. Tích lũy exciton trong đợt (f=1.05, S_n = {SN} MeV, {N} quỹ đạo × 60 T0) ===")
    with Pool() as pool:
        for kT in (0.1, 0.5):
            bs = [b for ev in pool.map(_run, [(cfg, kT, 27000 + i, T) for i in range(N)]) for b in analyse(ev, warm)[1]]
            pulses = [p * HW0 for b in bs for p in b['pulses']]
            n = len(bs)
            ge3 = sum(1 for b in bs if b['n'] >= 3) / n
            two = [sum(b['pulses'][:2]) * HW0 for b in bs if b['n'] >= 2]
            tot = [b['E'] * HW0 for b in bs]
            open_ = [e - SN for e in tot if e >= SN]
            k_open = []
            for b in bs:
                acc = 0.0
                for k, p in enumerate(b['pulses'], 1):
                    acc += p * HW0
                    if acc >= SN:
                        k_open.append(k)
                        break
            print(f"k_BT={kT*HW0:g} MeV: {n} đợt; tỉ lệ N≥3 = {ge3:.1%} (hình học (1−p)² với p=1/⟨N⟩: "
                  f"{(1-1/(sum(b['n'] for b in bs)/n))**2:.1%})")
            print(f"   xung đơn: TB {sum(pulses)/len(pulses):.2f} MeV, trung vị {pct(pulses,.5):.2f}, 10–90%: "
                  f"{pct(pulses,.1):.2f}–{pct(pulses,.9):.2f}; P(xung ≥ S_n) = {sum(1 for p in pulses if p >= SN)/len(pulses):.1%}")
            print(f"   E_exc(2) (đợt N≥2): TB {sum(two)/len(two):.2f} MeV, trung vị {pct(two,.5):.2f}; "
                  f"P(E_exc(2) ≥ S_n) = {sum(1 for e in two if e >= SN)/len(two):.1%}")
            print(f"   tổng đợt ≥ S_n: {len(open_)}/{n} = {len(open_)/n:.1%}; bước mở kênh k: TB "
                  f"{sum(k_open)/len(k_open):.2f}, trung vị {pct(k_open,.5)}; "
                  f"ε_max = E_exc − S_n: TB {sum(open_)/len(open_):.2f}, 10–90%: {pct(open_,.1):.2f}–{pct(open_,.9):.2f} MeV",
                  flush=True)
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
