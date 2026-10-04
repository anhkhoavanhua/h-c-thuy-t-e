"""Giải số hai đề xuất và một ghi chú phương pháp:
(A) Phổ Ericson P_n(ε) ∝ ε(1 − ε/E_max)^(n−2): trung bình đúng (tích phân số) so với ⟨ε⟩ = E_max/(n+1) của đề xuất.
(B) "Cắt nguồn lái trong pha tự do ⟹ ngân sách E_budget ≈ ½I(θ̇_gãy² − θ̇_cuối²) ≈ 15–25 MeV, loại đuôi trăm MeV":
    mô phỏng lái bật/tắt ở pha tự do; tổng Carnot mỗi đợt, ngân sách quay đo được, E_exc(3), bước mở kênh,
    và ⟨ε⟩ Ericson (đúng) trên từng đợt mở kênh với n = 2k ± 1.
(C) Mức hợp lệ bán cổ điển: tác dụng quay Iθ̇/ħ, năng lượng hạt so với ħω0/2, xung Carnot so với ħω0.

Đơn vị m=ω0=ħ=ℓ0=1, ħω0 = 10 MeV, f = 1.05, τ=10, γx=0.70, g=0.1, S_n = 8 MeV.
Chạy: python3 giai_so_ngan_sach_ericson.py [A|B|C]   (B ≈ 6 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import cfg_for
from giai_so_kramers_A import setup

TAU, GX, HW0, SN = 10.0, 0.70, 10.0, 8.0
T0 = L.T0


def ericson_mean(Emax, n, M=20000):
    """⟨ε⟩ của P(ε) ∝ ε(Emax − ε)^(n−2) trên [0, Emax], tích phân điểm giữa."""
    h = Emax / M
    num = den = 0.0
    for i in range(M):
        e = (i + 0.5) * h
        w = e * (Emax - e) ** (n - 2)
        num += e * w
        den += w
    return num / den


def part_a():
    print("=== A. Trung bình phổ Ericson ===")
    for Emax in (4.0,):
        for n in (3, 4, 5):
            print(f"E_max={Emax}, n={n}: ⟨ε⟩ số = {ericson_mean(Emax, n):.4f}; Beta(2,n−1): 2E/(n+1) = "
                  f"{2*Emax/(n+1):.4f}; đề xuất E/(n+1) = {Emax/(n+1):.4f}")
    print()


def _run(a):
    cfg, kT, seed, T, drive_free = a
    setup(cfg)
    wc = cfg[2]
    return L.run(1 / TAU, 1.05 * wc, kT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc, drive_free=drive_free)['events']


def bursts_with_budget(events, warm):
    """Đợt = chuỗi từ lần gãy kết thúc pha khóa ≥ T0 tới lần khóa đầu tiên kéo dài ≥ T0."""
    out = []
    for i, e in enumerate(events):
        if e[0] != 'gãy' or e[1] <= warm:
            continue
        prev = next((x for x in reversed(events[:i]) if x[0] == 'khóa'), None)
        if prev is None or e[1] - prev[1] < T0:
            continue
        pulses, w_end, closed = [], None, False
        for j in range(i + 1, len(events)):
            x = events[j]
            if x[0] == 'khóa':
                pulses.append(x[3] * HW0)
                nxt = next((y for y in events[j + 1:] if y[0] == 'gãy'), None)
                if nxt is None or nxt[1] - x[1] >= T0:
                    closed, w_end = nxt is not None, x[2]
                    break
        if closed:
            budget = 0.5 * (e[2]**2 - w_end**2) * HW0          # ½I(θ̇_gãy² − θ̇_cuối²), I = 1
            out.append(dict(pulses=pulses, total=sum(pulses), budget=budget))
    return out


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))]


def part_b(N=16):
    cfg = cfg_for(0.1)
    T, warm = 70 * T0, 10 * T0
    print(f"=== B. Ngân sách năng lượng: lái bật / tắt trong pha tự do (f=1.05, {N} quỹ đạo × 60 T0) ===")
    with Pool() as pool:
        for drive_free in (True, False):
            for kT in (0.1, 0.5):
                bs = [b for ev in pool.map(_run, [(cfg, kT, 29000 + i, T, drive_free) for i in range(N)])
                      for b in bursts_with_budget(ev, warm)]
                tot = [b['total'] for b in bs]
                bud = [b['budget'] for b in bs]
                e3 = [sum(b['pulses'][:3]) for b in bs if len(b['pulses']) >= 3]
                opened, eps_lo, eps_hi = [], [], []
                for b in bs:
                    acc = 0.0
                    for k, p in enumerate(b['pulses'], 1):
                        acc += p
                        if acc >= SN:
                            Emax = acc - SN
                            opened.append((k, Emax))
                            eps_lo.append(2 * Emax / (2 * k))           # n = 2k − 1 ⟹ n + 1 = 2k
                            eps_hi.append(2 * Emax / (2 * k + 2))       # n = 2k + 1 ⟹ n + 1 = 2k + 2
                            break
                tag = "BẬT" if drive_free else "TẮT"
                print(f"lái {tag} ở pha tự do, k_BT={kT*HW0:g} MeV: {len(bs)} đợt")
                print(f"   tổng Carnot/đợt: TB {sum(tot)/len(tot):.2f}, trung vị {pct(tot,.5):.2f}, 10–90%: "
                      f"{pct(tot,.1):.2f}–{pct(tot,.9):.2f}, max {max(tot):.1f} MeV")
                print(f"   ngân sách quay ½I(θ̇_gãy² − θ̇_cuối²): TB {sum(bud)/len(bud):.2f}, 10–90%: "
                      f"{pct(bud,.1):.2f}–{pct(bud,.9):.2f} MeV;  tổng/ngân sách TB {sum(tot)/max(sum(bud),1e-9):.2f}")
                if e3:
                    print(f"   E_exc(3): TB {sum(e3)/len(e3):.2f}, trung vị {pct(e3,.5):.2f} MeV")
                if opened:
                    em = [o[1] for o in opened]
                    print(f"   mở kênh: {len(opened)}/{len(bs)}; E_exc − S_n tại bước mở: TB {sum(em)/len(em):.2f}, "
                          f"trung vị {pct(em,.5):.2f} MeV;  ⟨ε⟩ Ericson (đúng) TB {sum(eps_hi)/len(eps_hi):.2f}–"
                          f"{sum(eps_lo)/len(eps_lo):.2f} MeV (n = 2k+1 … 2k−1)", flush=True)
    print()


def part_c():
    wc = cfg_for(0.1)[2]
    print("=== C. Mức hợp lệ bán cổ điển (ħ = ω0 = 1) ===")
    print(f"Mô-men động lượng pha tại ngưỡng: Iθ̇_gãy/ħ = {wc:.3f}  (số lượng tử ~3)")
    v = (1 / 3) * wc
    print(f"Năng lượng hạt đồng bộ cực đại ½m(κθ̇_gãy)² = {0.5*v*v:.3f} ħω0  so với điểm không ħω0/2 = 0.5")
    print(f"Xung Carnot TB ≈ 0.40 ħω0 < ħω0: truyền năng lượng nhỏ hơn một lượng tử dao động của hạt\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "C", "B"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
