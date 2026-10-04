"""Giải số đề xuất: (1) thời gian lưu dạng mũ do ghép cặp g, τ_lưu ∝ 1/g, g → 0 thì τ_lưu → ∞;
(2) "biến cố E gãy" = đợt vi chấn N_sub ≈ 5–10 xung, ΔE_biến cố ≈ N_sub × ΔE_xung ≈ 20–25 MeV ≈ thang GDR.

(A) Tất định (không nhiễu) theo g: có gãy không? (kiểm "g là nguồn kích hoạt tất định").
(B) Langevin 1 MeV theo g ∈ {0, 0.05, 0.1, 0.2}: thời gian lưu pha đồng bộ dài (chờ ≥ T0), CV phần dư.
(C) Cấu trúc đợt vi chấn (g = 0.1; 1 và 5 MeV): số lần khóa trong đợt, thời lượng, tổng Carnot mỗi đợt,
    Carnot mỗi xung trong đợt; so với GDR thực nghiệm.
Một đợt = chuỗi gãy/khóa giữa hai lần khóa kéo dài ≥ T0.

Đơn vị: m=ω0=ħ=ℓ0=1, 1 đơn vị năng lượng = ħω0 = 10 MeV; τ=10, γx=0.70, Ω0 = 1.05·θ̇_gãy(g).
Chạy: python3 giai_so_burst_g.py [A|B|C]   (≈ 10 phút tổng, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params
from giai_so_kramers_A import setup

TAU, GX, HW0 = 10.0, 0.70, 10.0
T0 = L.T0


def cfg_for(g):
    base, R = params(1.0, 3.0, I=1.0, g=g)
    return base, R, omega_crit_general(base, R)


def _run(a):
    cfg, kT, seed, T = a
    setup(cfg)
    wc = cfg[2]
    return L.run(1 / TAU, 1.05 * wc, kT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc)['events']


def analyse(events, warm):
    """Tách pha đồng bộ dài và đợt vi chấn từ chuỗi sự kiện."""
    dwell, last_lock = [], None
    for e in events:
        if e[0] == 'gãy':
            if last_lock is not None and e[1] - last_lock[1] >= T0 and e[1] > warm:
                dwell.append(e[1] - last_lock[1])      # một pha đồng bộ dài vừa kết thúc
            last_lock = None
        else:
            last_lock = e
    # đợt vi chấn: mở ở lần gãy kết thúc pha dài, khép ở lần khóa đầu tiên kéo dài ≥ T0
    bursts = []
    for i, e in enumerate(events):
        if e[0] != 'gãy' or e[1] <= warm:
            continue
        prev_lock = next((x for x in reversed(events[:i]) if x[0] == 'khóa'), None)
        if prev_lock is None or e[1] - prev_lock[1] < T0:
            continue
        # e là lần gãy mở đợt; gom các lần khóa tới lần khóa đầu tiên kéo dài ≥ T0 (hoặc tới cuối chuỗi)
        n, E, pulses, t_end, closed = 0, 0.0, [], None, False
        for j in range(i + 1, len(events)):
            x = events[j]
            if x[0] == 'khóa':
                n += 1
                E += x[3]
                pulses.append(x[3])
                t_end = x[1]
                nxt = next((y for y in events[j + 1:] if y[0] == 'gãy'), None)
                if nxt is None or nxt[1] - x[1] >= T0:
                    closed = nxt is not None
                    break
        if closed:
            bursts.append(dict(n=n, E=E, dur=t_end - e[1], pulses=pulses))
    return dwell, bursts


def part_a():
    print("=== A. Tất định (không nhiễu), 50 T0, theo g ===")
    for g in (0.0, 0.1, 0.2, 0.3):
        cfg = cfg_for(g)
        setup(cfg)
        ev = L.run(1 / TAU, 1.05 * cfg[2], 0.0, 50 * T0, 2e-3, 0, gx=GX, w_start=0.82 * cfg[2])['events']
        print(f"g={g}: θ̇_gãy = {cfg[2]:.4f}; số lần gãy = {sum(e[0]=='gãy' for e in ev)}")
    print("V_eff đúng (rút gọn từ L): −ΔE cosθ + (½kκ² − gκ) sin²θ;  đề xuất ghi ½(kκ² − gκ) sin²θ (hệ số g chỉ bằng ½)\n")


def part_b(N=8):
    print(f"=== B. Thời gian lưu pha đồng bộ dài theo g (1 MeV, {N} quỹ đạo × 60 T0 sau khởi động) ===")
    T, warm = 70 * T0, 10 * T0
    with Pool() as pool:
        for g in (0.0, 0.05, 0.1, 0.2):
            cfg = cfg_for(g)
            evs = pool.map(_run, [(cfg, 0.1, 17000 + i, T) for i in range(N)])
            dw = [w for ev in evs for w in analyse(ev, warm)[0]]
            nb = sum(sum(e[0] == 'gãy' and e[1] > warm for e in ev) for ev in evs)
            if len(dw) > 1:
                res = [w - T0 for w in dw]
                mr = sum(res) / len(res)
                cv = math.sqrt(sum((x - mr)**2 for x in res) / (len(res) - 1)) / mr
                print(f"g={g:<4}: θ̇_gãy={cfg[2]:.4f}; tổng gãy {nb}; pha dài {len(dw)}; τ_lưu TB = {sum(dw)/len(dw)/T0:.3f} T0 "
                      f"(phần dư TB {mr/T0:.3f} T0, CV {cv:.3f});  g·τ_lưu = {g*sum(dw)/len(dw)/T0:.3f}", flush=True)
            else:
                print(f"g={g:<4}: θ̇_gãy={cfg[2]:.4f}; tổng gãy {nb}; pha dài {len(dw)}", flush=True)
    print()


def part_c(N=16):
    print(f"=== C. Cấu trúc đợt vi chấn (g=0.1, {N} quỹ đạo × 60 T0), ħω0 = {HW0} MeV ===")
    T, warm = 70 * T0, 10 * T0
    cfg = cfg_for(0.1)
    with Pool() as pool:
        for kT in (0.1, 0.5):
            evs = pool.map(_run, [(cfg, kT, 19000 + i, T) for i in range(N)])
            bs = [b for ev in evs for b in analyse(ev, warm)[1]]
            ns = sorted(b['n'] for b in bs)
            Es = sorted(b['E'] * HW0 for b in bs)
            durs = [b['dur'] / T0 for b in bs]
            pul = [p * HW0 for b in bs for p in b['pulses']]
            q = lambda v, f: v[min(len(v) - 1, int(f * len(v)))]
            print(f"k_BT={kT*10:g} MeV: {len(bs)} đợt; số khóa/đợt TB {sum(ns)/len(ns):.2f} (trung vị {q(ns,.5)}, "
                  f"10–90%: {q(ns,.1)}–{q(ns,.9)}); thời lượng TB {sum(durs)/len(durs):.3f} T0 "
                  f"= {sum(durs)/len(durs)*2*math.pi*6.582e-22/HW0:.2e} s")
            print(f"   tổng Carnot/đợt: TB {sum(Es)/len(Es):.2f} MeV (trung vị {q(Es,.5):.2f}, 10–90%: {q(Es,.1):.2f}–{q(Es,.9):.2f});"
                  f"  Carnot/xung trong đợt: TB {sum(pul)/len(pul):.3f} MeV, xung lớn nhất TB/đợt "
                  f"{sum(max(b['pulses']) for b in bs)/len(bs)*HW0:.2f} MeV", flush=True)
    print("GDR thực nghiệm (thang): E ≈ 31.2A^-1/3 + 20.6A^-1/6 MeV; A=56: "
          f"{31.2*56**(-1/3)+20.6*56**(-1/6):.1f}, A=208: {31.2*208**(-1/3)+20.6*208**(-1/6):.1f} MeV; "
          f"80A^-1/3: A=16 {80*16**(-1/3):.1f}, A=56 {80*56**(-1/3):.1f}, A=208 {80*208**(-1/3):.1f} MeV")
    print(f"Độ rộng GDR Γ ≈ 4–5 MeV ⟹ ħ/Γ ≈ {6.582e-22/4.5:.2e} s   (đề xuất: ~1e-21 s)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
