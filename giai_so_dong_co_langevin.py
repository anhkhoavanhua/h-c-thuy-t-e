"""Giải số hai đề xuất cho động cơ E_kép:
(1) Mô-men lái tự kích Rayleigh–Van der Pol: Q_θ = γθ(1 − θ̇²/Ω0²)θ̇, Ω0 = f·θ̇_gãy ⟹ chu trình gãy/khóa tuần hoàn?
(2) Lực Langevin ξ(t), ⟨ξξ⟩ = 2γθ k_BT δ, k_BT = 1 MeV = 0.1 đơn vị ⟹ đo khuếch tán pha và tỉ lệ tái khóa (Monte Carlo).

Bộ thông số hạt nhân: I=1, m=k=1, ΔE=1, κ=1/3, R=3, g=0.1, e=0 (1 đơn vị năng lượng = ħω0 = 10 MeV, T0 = 2π).
Tích phân: RK4 cho phần tất định + cú hích ngẫu nhiên Euler–Maruyama sau mỗi bước (tách toán tử).
Trong pha liên kết, cú hích mô-men σ√dt·N chia cho quán tính hiệu dụng I + mκ²cos²θ, và x' = κcosθ·θ' được giữ đúng.

Chạy: python3 giai_so_dong_co_langevin.py [A|B|C|Q]   (Q: như B với nhiệt độ hiệu dụng ħω0/2 = 5 MeV)
         (A ≈ 1–2 phút; B, C chạy song song nhiều lõi, ≈ 5–10 phút)
"""
import math
import random
import sys
from multiprocessing import Pool

from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params
from giai_so_he_e_kep import rhs as rhs_free
from giai_so_lien_ket_e_gay import constrained_acc, energy_p, rk4_step
from giai_so_va_cham_tai_lien_ket import impulse

BASE, R = params(1.0, 3.0, I=1.0, g=0.1)
WC = omega_crit_general(BASE, R)
T0 = 2 * math.pi
KT = 0.1  # 1 MeV / (ħω0 = 10 MeV)


def drive(w, gam, Om):
    """Mô-men lái ròng Rayleigh–Van der Pol (đã gồm tiêu tán)."""
    return gam * (1.0 - w * w / (Om * Om)) * w


def acc_con(s, gam, Om, gx=0.0):
    return constrained_acc(s, **BASE, gth=-gam * (1.0 - s[1]**2 / Om**2), gx=gx)  # −gth·w = Q_drive


def f_con(s, gam, Om, gx=0.0):
    a_th, a_x, _ = acc_con(s, gam, Om, gx)
    return (s[1], a_th, s[3], a_x)


def f_free(s, gam, Om, gx=0.0):
    d = rhs_free(s, BASE['I'], BASE['m'], BASE['k'], BASE['g'], BASE['dE'])
    return (d[0], d[1] + drive(s[1], gam, Om) / BASE['I'], d[2], d[3] - gx * s[3] / BASE['m'])


def fval(s):
    return s[2] - BASE['kappa'] * math.sin(s[0])


def locate(step, s, dt, hit, n=40):
    lo, hi = 0.0, dt
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if hit(step(s, mid)):
            hi = mid
        else:
            lo = mid
    return hi, step(s, hi)


def run(gam, Om, kT, T, dt, seed, breakable=True, sample_every=None, w_start=None, gx=0.0, drive_free=True):
    """Một quỹ đạo. Trả về sự kiện, công lái, tổn hao Carnot, mẫu θ(t) (nếu yêu cầu).
    gx > 0: thêm cản γx·x' và lực Langevin ξx (⟨ξxξx⟩ = 2γx k_BT δ) lên hạt.
    drive_free=False: tắt mô-men lái khi hạt ở pha tự do (nhiễu ξθ vẫn giữ)."""
    rng = random.Random(seed)
    sig = math.sqrt(2.0 * gam * kT)
    sigx = math.sqrt(2.0 * gx * kT)
    w0 = Om if w_start is None else w_start
    s = [0.0, w0, 0.0, BASE['kappa'] * w0]
    t, locked, ref, events = 0.0, True, 0.0, []
    W_drive, carnot = 0.0, 0.0
    samples, next_sample = [], 0.0
    sc = lambda q, h: rk4_step(lambda z: f_con(z, gam, Om, gx), q, h)
    g_free = gam if drive_free else 0.0
    sf = lambda q, h: rk4_step(lambda z: f_free(z, g_free, Om, gx), q, h)
    lam_hit = lambda q: abs(acc_con(q, gam, Om, gx)[2]) >= R
    while t < T:
        if sample_every is not None and t >= next_sample:
            samples.append(s[0])
            next_sample += sample_every
        h = dt
        if locked:
            if breakable and lam_hit(s):
                events.append(('gãy', t, s[1]))
                locked, ref = False, 0.0
                continue
            s1 = sc(s, dt)
            if breakable and lam_hit(s1):
                h, s1 = locate(sc, s, dt, lam_hit)
        else:
            s1 = sf(s, dt)
            f1 = fval(s1)
            if ref == 0.0:
                ref = math.copysign(1.0, f1)
            elif f1 * ref <= 0:
                h, s1 = locate(sf, s, dt, lambda q: fval(q) * ref <= 0)
                s1, Lam, u = impulse(s1, BASE)       # e = 0: khóa dẻo
                c = math.cos(s1[0])
                meff = BASE['m'] * BASE['I'] / (BASE['I'] + BASE['m'] * BASE['kappa']**2 * c * c)
                carnot += 0.5 * meff * u * u
                events.append(('khóa', t + h, s1[1], 0.5 * meff * u * u))   # (loại, t, θ̇ sau khóa, Carnot)
                locked = True
        W_drive += drive(0.5 * (s[1] + s1[1]), gam if locked else g_free, Om) * (s1[0] - s[0])
        s, t = list(s1), t + h
        if sig > 0.0 or sigx > 0.0:                   # cú hích Langevin
            kick = sig * math.sqrt(h) * rng.gauss(0.0, 1.0)
            kx = sigx * math.sqrt(h) * rng.gauss(0.0, 1.0) if sigx > 0.0 else 0.0
            if locked:                                # ξx chiếu lên θ qua ∂x/∂θ = κcosθ
                c = math.cos(s[0])
                s[1] += (kick + BASE['kappa'] * c * kx) / (BASE['I'] + BASE['m'] * BASE['kappa']**2 * c * c)
                s[3] = BASE['kappa'] * c * s[1]
            else:
                s[1] += kick / BASE['I']
                s[3] += kx / BASE['m']
    return dict(events=events, W=W_drive, carnot=carnot, s=s, samples=samples, E=energy_p(s, BASE))


def part_a():
    print(f"=== A. Mô-men lái Van der Pol, không nhiễu (θ̇_gãy = {WC:.5f}) ===")
    for tau in (10.0, 50.265):
        gam = BASE['I'] / tau
        r0 = run(gam, 0.9 * WC, 0.0, 20 * T0, 1e-3, 0, breakable=False, w_start=0.5)
        print(f"τ={tau}: kiểm chu trình giới hạn (không gãy, Ω0=0.9θ̇_gãy={0.9*WC:.4f}, xuất phát θ̇=0.5): "
              f"θ̇(20T0) = {r0['s'][1]:.4f}")
        for f in (1.05, 1.10):
            Om = f * WC
            r = run(gam, Om, 0.0, 50 * T0, 1e-3, 0)
            ev = r['events']
            br = [e[1] for e in ev if e[0] == 'gãy']
            lk = [e for e in ev if e[0] == 'khóa']
            gaps = [b - a for a, b in zip(br, br[1:])]
            tail = gaps[-5:] if len(gaps) >= 5 else gaps
            wl = [e[2] for e in lk]
            print(f"  f={f}: Ω0={Om:.4f}; gãy={len(br)}, khóa={len(lk)} trong 50T0; "
                  f"khoảng cách giữa hai lần gãy: đầu {gaps[0] if gaps else float('nan'):.3f}, "
                  f"5 cuối {', '.join(f'{g:.3f}' for g in tail)} (TB {sum(tail)/max(len(tail),1):.3f} = "
                  f"{sum(tail)/max(len(tail),1)/T0:.3f} T0)")
            if wl:
                print(f"         θ̇ ngay sau khóa: min {min(wl):.4f}, max {max(wl):.4f} (so θ̇_gãy {WC:.4f});"
                      f"  công lái = {r['W']:.4f}, tổn hao Carnot = {r['carnot']:.4f}, "
                      f"chênh = {r['W']-r['carnot']:+.4f} (≈ ΔE hệ = {r['E'] - (0.5*(1+BASE['kappa']**2)*Om*Om - 1):+.4f})")
    print()


def _phase_worker(a):
    gam, Om, seed, T, dt, every, kT = a
    return run(gam, Om, kT, T, dt, seed, breakable=False, sample_every=every)['samples']


def part_b(N=64, kT=KT):
    print(f"=== B. Khuếch tán pha dưới nhiễu Langevin (k_BT = {kT} = {kT*10:g} MeV, không gãy, N = {N}) ===")
    T, dt, every = 10 * T0, 2e-3, 0.5
    with Pool() as pool:
        for tau in (1.0, 10.0):
            gam = BASE['I'] / tau
            Om = 0.9 * WC
            ref = run(gam, Om, 0.0, T, dt, 0, breakable=False, sample_every=every)['samples']
            sims = pool.map(_phase_worker, [(gam, Om, 1000 + i, T, dt, every, kT) for i in range(N)])
            L = min(len(ref), *(len(x) for x in sims))
            var = [sum((x[j] - ref[j])**2 for x in sims) / N for j in range(L)]
            t1 = next((j * every for j in range(L) if var[j] >= 1.0), float('nan'))
            slope = var[L - 1] / ((L - 1) * every)
            prop = BASE['I'] / (2 * kT * tau)                 # công thức đề xuất τ_deph = I/(2 k_BT τ)
            lc = 2 * BASE['I'] / (kT * tau)                   # pha chu trình giới hạn: 2I/(k_BT τ)
            print(f"τ={tau} (γθ={gam}): ⟨Δθ²⟩(10T0) = {var[L-1]:.4f}; ⟨Δθ²⟩/t = {slope:.5f} ⟹ τ_deph (⟨Δθ²⟩=1) đo = "
                  f"{t1 if not math.isnan(t1) else 1/slope:.2f} ({(t1 if not math.isnan(t1) else 1/slope)/T0:.2f} T0"
                  f"{'' if not math.isnan(t1) else ', ngoại suy'});"
                  f"  đề xuất I/(2k_BTτ) = {prop:.2f} ({prop/T0:.3f} T0);  chu trình giới hạn 2I/(k_BTτ) = {lc:.2f} ({lc/T0:.3f} T0)")
    print()


def _relock_worker(a):
    gam, Om, seed, T, dt, kT = a
    r = run(gam, Om, kT, T, dt, seed)
    ev = r['events']
    b = next((e[1] for e in ev if e[0] == 'gãy'), None)
    k = next((e[1] for e in ev if e[0] == 'khóa' and b is not None and e[1] > b), None)
    return (b, k, sum(e[0] == 'gãy' for e in ev), r['s'][1])


def part_c(N=64, kT=KT, tau_deph_measured=None):
    print(f"=== C. Monte Carlo tái khóa: lái Van der Pol + Langevin k_BT={kT:g} ({kT*10:g} MeV), Ω0 = 1.05θ̇_gãy, N = {N}, cửa sổ 10T0 ===")
    T, dt = 10 * T0, 2e-3
    with Pool() as pool:
        for tau in (1.0, 10.0):
            gam = BASE['I'] / tau
            Om = 1.05 * WC
            res = pool.map(_relock_worker, [(gam, Om, 5000 + i, T, dt, kT) for i in range(N)])
            broke = [r for r in res if r[0] is not None]
            rel = [r[1] - r[0] for r in broke if r[1] is not None]
            prop = BASE['I'] / (2 * kT * tau)
            if tau_deph_measured:
                meas = tau_deph_measured[tau]
                print(f"   τ={tau}: tái khóa trong τ_deph ĐO = {meas}: {sum(1 for r in res if r[0] is not None and r[1] is not None and r[1]-r[0] <= meas)}/{sum(1 for r in res if r[0] is not None)}")
            within = sum(1 for d in rel if d <= prop)
            nb = [r[2] for r in res]
            print(f"τ={tau}: gãy ≥1 lần: {len(broke)}/{N}; tái khóa sau lần gãy đầu: {len(rel)}/{len(broke)}; "
                  f"thời gian tái khóa: min {min(rel):.3f}, TB {sum(rel)/len(rel):.3f}, max {max(rel):.3f} "
                  f"({sum(rel)/len(rel)/T0:.3f} T0 TB);  trong τ_deph đề xuất {prop:.2f}: {within}/{len(broke)};"
                  f"  số lần gãy/quỹ đạo TB {sum(nb)/N:.1f}" if rel else f"τ={tau}: gãy {len(broke)}/{N}, không tái khóa")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c, "Q": lambda: (part_b(kT=0.5), part_c(kT=0.5, tau_deph_measured={1.0: 2.0, 10.0: 4.0}))}[name.upper()]()
