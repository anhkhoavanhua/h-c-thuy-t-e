"""Giải số hệ E_kép – hạt với liên kết holonom f = x − κ sinθ = 0 và mô hình "E gãy".

Pha 1 (có liên kết, n=1): tích phân hệ 4 biến (θ, θ', x, x'), mỗi bước giải hệ tuyến tính 3×3
    I θ'' + κ cosθ · λ = −ΔE sinθ + g x cosθ          (4)
    m x''        − λ    = −k x + g sinθ                 (5)
    x'' − κ cosθ · θ'' = −κ θ'² sinθ                    (f'' = 0)
để lấy θ'', x'', λ — độc lập với công thức (6), dùng để kiểm chứng (6).
Khi |λ| chạm ngưỡng R: định vị thời điểm gãy bằng chia đôi bước, rồi bỏ liên kết.
Pha 2 (tự do, n=2): tích phân (1)–(2) với điều kiện đầu = trạng thái lúc gãy.

Chạy: python3 giai_so_lien_ket_e_gay.py   (chỉ thư viện chuẩn)
"""
import math

from giai_so_he_e_kep import energy, rhs as rhs_free

P = dict(I=1.0, m=1.0, k=4.0, g=0.1, dE=0.5, kappa=0.2)


def solve3(a, b):
    """Khử Gauss 3×3 có chọn trụ."""
    a = [row[:] + [bi] for row, bi in zip(a, b)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(a[r][c]))
        a[c], a[piv] = a[piv], a[c]
        for r in range(c + 1, 3):
            f = a[r][c] / a[c][c]
            a[r] = [u - f * v for u, v in zip(a[r], a[c])]
    xs = [0.0] * 3
    for r in (2, 1, 0):
        xs[r] = (a[r][3] - sum(a[r][j] * xs[j] for j in range(r + 1, 3))) / a[r][r]
    return xs


def constrained_acc(s, I, m, k, g, dE, kappa, gth=0.0, gx=0.0):
    """gth, gx: hệ số cản Rayleigh Φ = ½gth·θ'² + ½gx·x'² (mặc định 0 = không cản)."""
    th, w, x, v = s
    c, sn = math.cos(th), math.sin(th)
    return solve3([[I, 0.0, kappa * c], [0.0, m, -1.0], [-kappa * c, 1.0, 0.0]],
                  [-dE * sn + g * x * c - gth * w, -k * x + g * sn - gx * v, -kappa * w * w * sn])  # θ'', x'', λ


def rhs_con(s, p):
    a_th, a_x, _ = constrained_acc(s, **p)
    return (s[1], a_th, s[3], a_x)


def lam_closed(th, w, I, m, k, g, dE, kappa):
    """λ(θ, θ') sau khi khử θ'' bằng phương trình rút gọn (dạng đóng, kiểm bằng sympy)."""
    c, sn = math.cos(th), math.sin(th)
    return (-dE * kappa * m * c - I * w * w * kappa * m - I * g + I * k * kappa
            + g * kappa**2 * m * c * c) * sn / (I + kappa**2 * m * c * c)


def lam_eq6(th, w, a_th, m, k, g, kappa, **_):
    """Công thức (6) của bản gốc, dùng θ'' lấy từ nghiệm số."""
    return m * kappa * (a_th * math.cos(th) - w * w * math.sin(th)) + (k * kappa - g) * math.sin(th)


def rk4_step(f, s, h):
    k1 = f(s)
    k2 = f([a + h / 2 * b for a, b in zip(s, k1)])
    k3 = f([a + h / 2 * b for a, b in zip(s, k2)])
    k4 = f([a + h * b for a, b in zip(s, k3)])
    return [a + h / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4)]


def energy_p(s, p):
    return energy((0.0, *s), p['I'], p['m'], p['k'], p['g'], p['dE'])


def run(w_start, R, p=P, t_end=20.0, dt=1e-4, t_after=10.0):
    """Pha 1 có liên kết tới khi |λ| = R (hoặc t_end), rồi pha 2 tự do trong t_after."""
    fc = lambda s: rhs_con(s, p)
    s = [0.0, w_start, 0.0, p['kappa'] * w_start]  # thỏa f=0 và f'=0 tại t=0
    t, E0 = 0.0, energy_p(s, p)
    stats = dict(err6=0.0, errc=0.0, drift=0.0, lam_max=0.0, work_x=0.0, work_th=0.0)
    lam = constrained_acc(s, **p)[2]
    while t < t_end:
        s_new = rk4_step(fc, s, dt)
        a_th, a_x, lam_new = constrained_acc(s_new, **p)
        if abs(lam_new) >= R:
            lo, hi = 0.0, dt  # chia đôi bước để tìm đúng t_gãy với |λ| = R
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                if abs(constrained_acc(rk4_step(fc, s, mid), **p)[2]) >= R:
                    hi = mid
                else:
                    lo = mid
            s_b = rk4_step(fc, s, hi)
            _accumulate_work(stats, s, s_b, lam, constrained_acc(s_b, **p)[2], p)
            t += hi
            s = s_b
            break
        _accumulate_work(stats, s, s_new, lam, lam_new, p)
        s, t, lam = s_new, t + dt, lam_new
        th, w, x, _ = s
        stats['err6'] = max(stats['err6'], abs(lam - lam_eq6(th, w, a_th, **p)))
        stats['errc'] = max(stats['errc'], abs(lam - lam_closed(th, w, **p)))
        stats['drift'] = max(stats['drift'], abs(x - p['kappa'] * math.sin(th)))
        stats['lam_max'] = max(stats['lam_max'], abs(lam))
        stats['Econ'] = max(stats.get('Econ', 0.0), abs(energy_p(s, p) - E0) / abs(E0))
    else:
        return dict(broke=False, t=t, stats=stats, E0=E0)

    # ----- Tại thời điểm gãy -----
    th, w, x, v = s
    a_th_c, a_x_c, lam_b = constrained_acc(s, **p)
    _, a_th_f, _, a_x_f = rhs_free(s, p['I'], p['m'], p['k'], p['g'], p['dE'])
    res = dict(broke=True, t=t, th=th, w=w, x=x, v=v, lam=lam_b, stats=stats, E0=E0,
               E_b=energy_p(s, p), a_x_con=a_x_c, a_x_free=a_x_f, a_th_con=a_th_c, a_th_free=a_th_f)

    # ----- Pha 2: tự do (1)–(2) -----
    ff = lambda q: rhs_free(q, p['I'], p['m'], p['k'], p['g'], p['dE'])
    q, tt, dev_max, Ex_max, Ef = list(s), 0.0, 0.0, 0.0, 0.0
    marks = {1e-1: None, 1e-2: None, 1e-3: None}
    h = 1e-5
    while tt < 1e-1 - 1e-12:  # đoạn ngắn đầu tiên bước nhỏ để đo Δp(Δt)
        q = rk4_step(ff, q, h)
        tt += h
        for d in marks:
            if marks[d] is None and tt >= d - 1e-12:
                marks[d] = p['m'] * (q[3] - v)
    while tt < t_after:
        q = rk4_step(ff, q, dt)
        tt += dt
        dev_max = max(dev_max, abs(q[2] - p['kappa'] * math.sin(q[0])))
        Ex_max = max(Ex_max, 0.5 * p['m'] * q[3]**2 + 0.5 * p['k'] * q[2]**2)
        Ef = max(Ef, abs(energy_p(q, p) - res['E_b']) / abs(res['E_b']))
    res.update(dp=marks, dev_max=dev_max, Ex_b=0.5 * p['m'] * v**2 + 0.5 * p['k'] * x**2,
               Ex_max=Ex_max, Efree=Ef)
    return res


def _accumulate_work(st, s0, s1, l0, l1, p):
    """Công của lực liên kết lên hạt ∫λ dx và lên pha ∫(−λκcosθ) dθ (quy tắc hình thang)."""
    lm = 0.5 * (l0 + l1)
    st['work_x'] += lm * (s1[2] - s0[2])
    cm = 0.5 * (math.cos(s0[0]) + math.cos(s1[0]))
    st['work_th'] += -lm * p['kappa'] * cm * (s1[0] - s0[0])


def lam_rot_max(w, p=P, n=20000):
    """Ước lượng max|λ| trên một vòng quay (θ' lấy từ bảo toàn năng lượng rút gọn)."""
    I, m, k, g, dE, ka = p['I'], p['m'], p['k'], p['g'], p['dE'], p['kappa']
    E = 0.5 * (I + m * ka**2) * w * w - dE  # tại θ = 0, hạt bị kéo theo với x' = κθ'
    best = 0.0
    for i in range(n):
        th = 2 * math.pi * i / n
        c, s = math.cos(th), math.sin(th)
        Ueff = -dE * c + (0.5 * k * ka**2 - g * ka) * s * s
        w2 = 2 * (E - Ueff) / (I + m * ka**2 * c * c)
        if w2 > 0:
            best = max(best, abs(lam_closed(th, math.sqrt(w2), **p)))
    return best


if __name__ == "__main__":
    print("Thông số:", ", ".join(f"{k}={v}" for k, v in P.items()))
    print("\n=== 1. Kiểm chứng pha có liên kết (không gãy, R = ∞), θ'(0)=5 ===")
    r = run(5.0, float('inf'), t_end=20.0)
    st = r['stats']
    print(f"max|λ_DAE − λ_(6)|        = {st['err6']:.2e}")
    print(f"max|λ_DAE − λ_dạng-đóng|  = {st['errc']:.2e}")
    print(f"max|x − κ sinθ| (trôi liên kết) = {st['drift']:.2e}")
    print(f"trôi năng lượng tương đối = {st['Econ']:.2e}")
    print(f"max|λ| trên [0,20] = {st['lam_max']:.5f}")
    print(f"∫λ dx (công lên hạt) = {st['work_x']:+.6f};  ∫(−λκcosθ)dθ (công lên pha) = {st['work_th']:+.6f};"
          f"  tổng = {st['work_x']+st['work_th']:+.2e}")

    print("\n=== 2. max|λ| theo tốc độ quay (so với hệ số θ'² giải tích mκ) ===")
    for w in (2.0, 3.0, 5.0, 8.0, 12.0):
        lm = lam_rot_max(w)
        print(f"θ'(0)={w:5.1f}: max|λ| = {lm:9.4f}   mκθ'(0)² = {P['m']*P['kappa']*w*w:9.4f}   tỉ số = {lm/(P['m']*P['kappa']*w*w):.4f}")

    print("\n=== 3. Mô phỏng E gãy: R = 3.0 ===")
    for w in (3.0, 5.0, 8.0):
        r = run(w, 3.0)
        if not r['broke']:
            print(f"θ'(0)={w}: KHÔNG gãy trên [0,20] (max|λ| = {r['stats']['lam_max']:.4f} < R)")
            continue
        st = r['stats']
        print(f"θ'(0)={w}: GÃY tại t = {r['t']:.6f}, θ = {r['th']:.6f} rad ({math.degrees(r['th']) % 360:.2f}°), "
              f"θ' = {r['w']:.5f}, λ = {r['lam']:+.6f}")
        print(f"   x = {r['x']:+.6f} (κsinθ = {P['kappa']*math.sin(r['th']):+.6f}), p = m x' = {P['m']*r['v']:+.6f}")
        print(f"   bước nhảy gia tốc hạt: x''_liên-kết = {r['a_x_con']:+.6f} → x''_tự-do = {r['a_x_free']:+.6f}"
              f"  (Δx'' = {r['a_x_free']-r['a_x_con']:+.6f}, −λ/m = {-r['lam']/P['m']:+.6f})")
        print(f"   bước nhảy gia tốc pha:  θ''_liên-kết = {r['a_th_con']:+.6f} → θ''_tự-do = {r['a_th_free']:+.6f}")
        print(f"   Δp = p(t_gãy+Δt) − p(t_gãy):  Δt=0.1 → {r['dp'][1e-1]:+.3e},  0.01 → {r['dp'][1e-2]:+.3e},"
              f"  0.001 → {r['dp'][1e-3]:+.3e}   (chặn trên R·Δt = {3.0*0.1:.1e}, {3.0*0.01:.1e}, {3.0*0.001:.1e})")
        print(f"   E toàn phần: đầu = {r['E0']:.10f}, lúc gãy = {r['E_b']:.10f} "
              f"(Δ = {r['E_b']-r['E0']:+.1e}); trôi sau gãy = {r['Efree']:.1e}")
        print(f"   ∫λ dx đến lúc gãy = {st['work_x']:+.6f};  ∫(−λκcosθ)dθ = {st['work_th']:+.6f};  tổng = {st['work_x']+st['work_th']:+.1e}")
        print(f"   năng lượng riêng của hạt ½mx'²+½kx²: lúc gãy = {r['Ex_b']:.6f}, max 10 đơn vị t sau = {r['Ex_max']:.6f}")
        print(f"   độ lệch khỏi quỹ đạo đồng bộ max|x − κsinθ| sau gãy = {r['dev_max']:.6f} (κ = {P['kappa']})")
