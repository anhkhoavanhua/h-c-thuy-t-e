"""Giải số đề xuất "nhịp gãy = quá trình Poisson băng ngưỡng (Rice) quanh θ̇_d".

(A) Dự đoán Rice/OU theo đúng công thức đề xuất: σ_u² = γθk_BT/(I_eff|F'|) (+ phần nhiễu hạt chiếu qua κcosθ),
    Δu = θ̇_gãy − θ̇_d, thừa số exp(−Δu²/2σ_u²) và tỉ số r(1 MeV)/r(5 MeV) nó dự đoán.
(B) Pha đồng bộ KHÔNG cho gãy, đo theo k_BT: θ̇ tại mặt cắt θ ≡ 0 (mod 2π) — trung bình, độ lệch chuẩn;
    lề ứng suất mỗi vòng M = max|λ|/R; tần suất vòng có M ≥ 1 (tốc độ vượt ngưỡng "nếu được phép gãy").
    Kiểm "điểm dừng dời lên gần ngưỡng khi tăng k_BT" và so tốc độ vượt ngưỡng với r đo được.
(C) Phân bố thời gian chờ khóa→gãy (cho phép gãy): hệ số biến thiên CV (=1 nếu mũ/Poisson),
    hàm sống P(t>k·t̄) so với e^{−k}, tỉ lệ gãy lại tức thì (< 0.2 T0).

Chạy: python3 giai_so_rice_poisson.py [A|B|C]   (B, C ≈ 3–4 phút mỗi phần, 8 lõi)
"""
import math
import random
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_kramers_A import setup, std_cfg
from giai_so_lien_ket_e_gay import rk4_step

TAU, GX = 10.0, 0.70
KAP = 1 / 3
I_EFF = 1 + 0.5 * KAP**2


def part_a(cfg):
    wc = cfg[2]
    gth = 1 / TAU
    Om = 1.05 * wc
    wd = Om * math.sqrt(1 - GX * KAP**2 / (2 * gth))
    fp = 2 * gth * wd**2 / Om**2
    du = wc - wd
    print(f"=== A. Rice/OU theo công thức đề xuất (τ={TAU}, γx={GX}, θ̇_d={wd:.4f}, Δu={du:.4f}, |F'|={fp:.4f}) ===")
    out = {}
    for kT in (0.1, 0.5):
        s_th = gth * kT / (I_EFF * fp)
        s_x = 0.5 * GX * KAP**2 * kT / (I_EFF * fp)
        for name, s2 in (("chỉ ξθ (đề xuất)", s_th), ("ξθ + ξx", s_th + s_x)):
            ex = du * du / (2 * s2)
            out[(kT, name)] = ex
            print(f"k_BT={kT*10:g} MeV, {name}: σ_u² = {s2:.4f}, σ_u = {math.sqrt(s2):.4f}, "
                  f"Δu²/2σ² = {ex:.3f}, exp(−·) = {math.exp(-ex):.4f}")
    for name in ("chỉ ξθ (đề xuất)", "ξθ + ξx"):
        ratio = math.exp(-out[(0.1, name)] + out[(0.5, name)])
        print(f"Rice dự đoán r(1 MeV)/r(5 MeV) [{name}] = {ratio:.3f}   (đo được: 0.1576/0.1403 = {0.1576/0.1403:.3f})")
    print()
    return wd


def _locked(a):
    cfg, kT, seed, T, warm, dt = a
    setup(cfg)
    wc = cfg[2]
    gam, Om = 1 / TAU, 1.05 * wc
    rng = random.Random(seed)
    sig, sigx = math.sqrt(2 * gam * kT), math.sqrt(2 * GX * kT)
    base = L.BASE
    s = [0.0, 0.82 * wc, 0.0, base['kappa'] * 0.82 * wc]
    t, rev_start, M, sec_w, margins = 0.0, 0.0, 0.0, [], []
    n_rev = 0
    while t < T:
        s1 = rk4_step(lambda z: L.f_con(z, gam, Om, GX), s, dt)
        c = math.cos(s1[0])
        ieff = base['I'] + base['m'] * base['kappa']**2 * c * c
        s1[1] += (sig * math.sqrt(dt) * rng.gauss(0, 1) + base['kappa'] * c * sigx * math.sqrt(dt) * rng.gauss(0, 1)) / ieff
        s1[3] = base['kappa'] * c * s1[1]
        t += dt
        lam = abs(L.acc_con(s1, gam, Om, GX)[2]) / L.R
        M = max(M, lam)
        if math.floor(s1[0] / (2 * math.pi)) != math.floor(s[0] / (2 * math.pi)):   # qua mặt cắt θ ≡ 0
            if t > warm:
                sec_w.append(s1[1])
                margins.append(M)
                n_rev += 1
            M = 0.0
        s = s1
    return sec_w, margins, T - warm


def part_b(cfg, N=8):
    print(f"=== B. Pha đồng bộ không cho gãy (τ={TAU}, γx={GX}); {N} quỹ đạo × 60 T0 sau khởi động ===")
    print(f"{'k_BT':>6} {'θ̇(θ≡0) TB':>11} {'σ':>7} {'(θ̇_gãy−TB)/σ':>14} {'M TB':>7} {'σ_M':>7} "
          f"{'P(M≥1)/vòng':>12} {'vòng/T0':>8} {'vượt ngưỡng /T0':>16}")
    T, warm = 70 * L.T0, 10 * L.T0
    with Pool() as pool:
        for kT in (0.1, 0.3, 0.5, 0.8):
            res = pool.map(_locked, [(cfg, kT, 13000 + i, T, warm, 2e-3) for i in range(N)])
            ws = [w for r in res for w in r[0]]
            ms = [m for r in res for m in r[1]]
            span = sum(r[2] for r in res) / L.T0
            mw = sum(ws) / len(ws)
            sw = math.sqrt(sum((w - mw)**2 for w in ws) / (len(ws) - 1))
            mm = sum(ms) / len(ms)
            sm = math.sqrt(sum((m - mm)**2 for m in ms) / (len(ms) - 1))
            p = sum(1 for m in ms if m >= 1) / len(ms)
            print(f"{kT*10:4g}MeV {mw:11.4f} {sw:7.4f} {(cfg[2]-mw)/sw:14.2f} {mm:7.4f} {sm:7.4f} {p:12.4f} "
                  f"{len(ms)/span:8.3f} {p*len(ms)/span:16.3f}", flush=True)
    print("(so với r đo được ở chế độ cho phép gãy: 0.158 / 0.143 / 0.140 / 0.146 mỗi T0)\n")


def _waits(a):
    cfg, kT, seed, T, warm = a
    setup(cfg)
    wc = cfg[2]
    r = L.run(1 / TAU, 1.05 * wc, kT, T, 2e-3, seed, gx=GX, w_start=0.82 * wc)
    out, last = [], None
    for e in r['events']:
        if e[0] == 'khóa':
            last = e[1]
        elif e[0] == 'gãy' and last is not None and e[1] > warm:
            out.append(e[1] - last)
            last = None
    return out


def part_c(cfg, N=16):
    print(f"=== C. Phân bố thời gian chờ khóa→gãy (τ={TAU}, γx={GX}); {N} quỹ đạo × 60 T0 ===")
    T, warm = 70 * L.T0, 10 * L.T0
    with Pool() as pool:
        for kT in (0.1, 0.5):
            ws = [w / L.T0 for r in pool.map(_waits, [(cfg, kT, 15000 + i, T, warm) for i in range(N)]) for w in r]
            n = len(ws)
            m = sum(ws) / n
            sd = math.sqrt(sum((w - m)**2 for w in ws) / (n - 1))
            imm = sum(1 for w in ws if w < 0.2) / n
            longw = [w for w in ws if w >= 1.0]
            ml = sum(longw) / len(longw)
            sdl = math.sqrt(sum((w - ml)**2 for w in longw) / (len(longw) - 1))
            # kiểm dạng mũ cho phần dư (w − 1) của các lần chờ ≥ T0
            res = [w - 1.0 for w in longw]
            mr = sum(res) / len(res)
            surv = [(k, sum(1 for x in res if x > k * mr) / len(res), math.exp(-k)) for k in (0.5, 1, 2, 3)]
            print(f"k_BT={kT*10:g} MeV: {n} lần chờ, TB {m:.3f} T0, CV = {sd/m:.3f}; gãy lại tức thì (<0.2 T0): {imm:.1%}")
            print(f"   chờ ≥ T0: {len(longw)} lần, TB {ml:.3f} T0, CV {sdl/ml:.3f}; phần dư (w−1): TB {mr:.3f} T0, CV "
                  f"{math.sqrt(sum((x-mr)**2 for x in res)/(len(res)-1))/mr:.3f}; "
                  f"P(>k·TB) đo/mũ: " + ", ".join(f"k={k}: {a:.3f}/{b:.3f}" for k, a, b in surv), flush=True)
    print()


if __name__ == "__main__":
    cfg = std_cfg()
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": lambda: part_a(cfg), "B": lambda: part_b(cfg), "C": lambda: part_c(cfg)}[name.upper()]()
