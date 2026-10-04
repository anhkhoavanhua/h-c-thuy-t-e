"""Giải số: (A) tốc độ gãy tới hạn khi κ = ℓ·E_kép/E_vũtrụ, R = E_vũtrụ/ℓ;
(B) va chạm tái liên kết (khóa pha đột ngột) và chu kỳ gãy → tự do → tái khóa.

Phương trình va chạm đúng cho hệ 2 bậc tự do (xung lực Λ tác dụng lên CẢ hạt và pha):
    m (v⁺ − v⁻) = Λ,   I (θ'⁺ − θ'⁻) = −Λ κ cosθ,   v⁺ = κ cosθ θ'⁺
 ⟹ Λ = m_eff (κ cosθ θ'⁻ − v⁻),  m_eff = m I / (I + m κ² cos²θ)
    Năng lượng mất (định lý Carnot) = ½ m_eff (κ cosθ θ'⁻ − v⁻)² ≥ 0.

Chạy: python3 giai_so_va_cham_tai_lien_ket.py   (chỉ thư viện chuẩn)
"""
import math

from giai_so_he_e_kep import rhs as rhs_free
from giai_so_lien_ket_e_gay import constrained_acc, energy_p, lam_rot_max, rk4_step


def phys(p):
    """Bỏ khóa R (ngưỡng) để truyền vào các hàm chỉ nhận thông số vật lý."""
    return {k: v for k, v in p.items() if k != 'R'}


def omega_crit(p, lo=0.0, hi=100.0):
    """θ'(0) nhỏ nhất để max|λ| trên quỹ đạo đồng bộ chạm R (chia đôi)."""
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if lam_rot_max(mid, phys(p), n=4000) >= p['R']:
            hi = mid
        else:
            lo = mid
    return hi


def part_a():
    print("=== A. Tốc độ gãy tới hạn: κ = ℓ E_kép/E_vt, R = E_vt/ℓ ⟹ R·κ = E_kép ===")
    base = dict(I=1.0, m=1.0, k=4.0, g=0.1, dE=0.5)
    ell, Ek = 1.0, 0.5
    print(f"ℓ={ell}, m={base['m']}, E_kép={Ek} (= ΔE), I={base['I']}, k={base['k']}, g={base['g']}")
    print(f"{'E_vt':>6} {'κ':>8} {'R':>6} {'R·κ':>6} {'θ̇ công thức':>13} {'θ̇ chính xác':>13} {'sai lệch':>9}")
    for Ev in (1.0, 2.0, 5.0, 10.0, 20.0):
        p = dict(base, kappa=ell * Ek / Ev, R=Ev / ell)
        wf = Ev / (ell * math.sqrt(base['m'] * Ek))
        wc = omega_crit(p)
        print(f"{Ev:6.1f} {p['kappa']:8.4f} {p['R']:6.1f} {p['R']*p['kappa']:6.2f} {wf:13.5f} {wc:13.5f} {(wf-wc)/wc:+9.2%}")
    # kiểm tra lại bằng mô phỏng: ngay dưới / trên ngưỡng
    p = dict(base, kappa=ell * Ek / 5.0, R=5.0 / ell)
    wc = omega_crit(p)
    for w in (0.99 * wc, 1.01 * wc):
        br = simulate(p, w, t_end=20.0, max_events=1)
        print(f"  kiểm chứng E_vt=5: θ'(0)={w:.4f} ({w/wc:.2f}×θ̇_tới hạn) → "
              f"{'GÃY tại t=%.4f' % br['events'][0]['t'] if br['events'] else 'không gãy trên [0,20]'}")
    print()


def impulse(s, p, exact=True, e=0.0):
    """Trạng thái sau va chạm với hệ số hồi phục e (e=0: khóa dẻo, ḟ⁺ = −e ḟ⁻).
    exact=False dùng công thức bản gốc (θ' giữ nguyên, e=0)."""
    th, w, x, v = s
    c = math.cos(th)
    u = p['kappa'] * c * w - v                           # vận tốc tương đối trước va chạm
    if exact:
        meff = p['m'] * p['I'] / (p['I'] + p['m'] * p['kappa']**2 * c * c)
        Lam = (1.0 + e) * meff * u
        w_new = w - Lam * p['kappa'] * c / p['I']
        return [th, w_new, x, v + Lam / p['m']], Lam, u
    Lam, w_new = p['m'] * u, w
    return [th, w_new, x, p['kappa'] * c * w_new], Lam, u


def locate(step, s, dt, hit):
    """Chia đôi trong một bước để tìm thời điểm sự kiện; trả (h, trạng thái)."""
    lo, hi = 0.0, dt
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if hit(step(s, mid)):
            hi = mid
        else:
            lo = mid
    return hi, step(s, hi)


def simulate(p, w_start, t_end=30.0, dt=1e-4, max_events=10**9, e=0.0, lock_tol=1e-3):
    """Chu kỳ: đồng bộ (n=1) → gãy khi |λ|=R → tự do (n=2) → va chạm khi f=0 → ...
    e=0: khóa ngay. e>0: nảy (ḟ⁺ = −e ḟ⁻), chỉ khóa khi |u| < lock_tol (cắt chuỗi Zeno)."""
    fc = lambda s: rhs_con_p(s, p)
    gth, gx = p.get('gth', 0.0), p.get('gx', 0.0)   # cản Rayleigh (tùy chọn)

    def ff(s):
        d = rhs_free(s, p['I'], p['m'], p['k'], p['g'], p['dE'])
        return (d[0], d[1] - gth * s[1] / p['I'], d[2], d[3] - gx * s[3] / p['m'])
    step_c = lambda s, h: rk4_step(fc, s, h)
    step_f = lambda s, h: rk4_step(ff, s, h)
    fval = lambda s: s[2] - p['kappa'] * math.sin(s[0])
    s, t, locked, events = [0.0, w_start, 0.0, p['kappa'] * w_start], 0.0, True, []
    ref_sign = 0.0
    while t < t_end and len(events) < max_events:
        if locked:
            if abs(constrained_acc(s, **phys(p))[2]) >= p['R']:          # gãy ngay khi vừa khóa
                events.append(dict(kind='gãy', t=t, th=s[0], w=s[1], lam=constrained_acc(s, **phys(p))[2], E=energy_p(s, p)))
                locked, ref_sign = False, 0.0
                continue
            s1 = step_c(s, dt)
            if abs(constrained_acc(s1, **phys(p))[2]) >= p['R']:
                h, s1 = locate(step_c, s, dt, lambda q: abs(constrained_acc(q, **phys(p))[2]) >= p['R'])
                t += h
                s = s1
                events.append(dict(kind='gãy', t=t, th=s[0], w=s[1], lam=constrained_acc(s, **phys(p))[2], E=energy_p(s, p)))
                locked, ref_sign = False, 0.0
                continue
            s, t = s1, t + dt
        else:
            s1 = step_f(s, dt)
            f1 = fval(s1)
            if ref_sign == 0.0:
                ref_sign = math.copysign(1.0, f1)
            elif f1 * ref_sign <= 0:                                 # quỹ đạo tự do cắt lại mặt f = 0
                h, s1 = locate(step_f, s, dt, lambda q: fval(q) * ref_sign <= 0)
                t += h
                E_before = energy_p(s1, p)
                u0 = p['kappa'] * math.cos(s1[0]) * s1[1] - s1[3]
                lock = e == 0.0 or abs(u0) < lock_tol
                s_new, Lam, u = impulse(s1, p, e=0.0 if lock else e)
                s_paste, Lam_paste, _ = impulse(s1, p, exact=False)
                th, c = s1[0], math.cos(s1[0])
                meff = p['m'] * p['I'] / (p['I'] + p['m'] * p['kappa']**2 * c * c)
                events.append(dict(
                    kind='khóa' if lock else 'nảy', t=t, th=th, f=fval(s1), v_minus=s1[3], w_minus=s1[1], u=u,
                    v_plus=s_new[3], w_plus=s_new[1], Lam=Lam, Lam_paste=Lam_paste, E_before=E_before,
                    E_after=energy_p(s_new, p), carnot=0.5 * (1.0 - (0.0 if lock else e)**2) * meff * u * u,
                    dT_paste=0.5 * p['m'] * (s1[3]**2 - s_paste[3]**2),
                    E_after_paste=energy_p(s_paste, p)))
                s, locked, ref_sign = s_new, lock, 0.0
                continue
            s, t = s1, t + dt
    return dict(events=events, t=t, s=s, locked=locked, E=energy_p(s, p))


def rhs_con_p(s, p):
    a_th, a_x, _ = constrained_acc(s, **phys(p))
    return (s[1], a_th, s[3], a_x)


def part_b(I, t_end=30.0):
    p = dict(I=I, m=1.0, k=4.0, g=0.1, dE=0.5, kappa=0.2, R=3.0)
    print(f"=== B. Chu kỳ gãy / tái khóa: {', '.join(f'{k}={v}' for k, v in p.items())}, θ'(0)=5, t∈[0,{t_end}] ===")
    r = simulate(p, 5.0, t_end=t_end)
    E0 = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * 25.0 - p['dE']
    print(f"E ban đầu = {E0:.6f};  θ̇ tới hạn (chính xác) = {omega_crit(p):.5f}")
    for i, e in enumerate(r['events'], 1):
        if e['kind'] == 'gãy':
            print(f"{i:2d}. GÃY  t={e['t']:9.5f}  θ={math.degrees(e['th']) % 360:7.2f}°  θ'={e['w']:8.5f}  λ={e['lam']:+.5f}  E={e['E']:.6f}")
        else:
            print(f"{i:2d}. KHÓA t={e['t']:9.5f}  θ={math.degrees(e['th']) % 360:7.2f}°  |f|={abs(e['f']):.1e}  "
                  f"v⁻={e['v_minus']:+.5f} → v⁺={e['v_plus']:+.5f}   θ'⁻={e['w_minus']:.5f} → θ'⁺={e['w_plus']:.5f}")
            print(f"      Λ đúng={e['Lam']:+.6f}  Λ bản gốc m(κcθ̇−v⁻)={e['Lam_paste']:+.6f}  (tỉ số {e['Lam_paste']/e['Lam']:.5f})")
            print(f"      E mất thực = {e['E_before']-e['E_after']:.6f}   Carnot ½m_eff u² = {e['carnot']:.6f}   "
                  f"ΔT bản gốc ½m(v⁻²−v⁺²) = {e['dT_paste']:+.6f}   (E mất nếu giữ θ' như bản gốc: {e['E_before']-e['E_after_paste']:+.6f})")
    n_b = sum(e['kind'] == 'gãy' for e in r['events'])
    n_k = sum(e['kind'] == 'khóa' for e in r['events'])
    print(f"Kết thúc t={r['t']:.3f}: {n_b} lần gãy, {n_k} lần khóa; trạng thái cuối "
          f"{'ĐỒNG BỘ (n=1)' if r['locked'] else 'TỰ DO (n=2)'}, θ'={r['s'][1]:.5f}, E={r['E']:.6f} "
          f"(mất tổng {E0 - r['E']:.6f} = {(E0 - r['E'])/E0:.2%})\n")


if __name__ == "__main__":
    part_a()
    part_b(1.0)
    part_b(0.05, t_end=10.0)
