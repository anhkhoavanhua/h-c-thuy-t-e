"""Giải số hệ E_kép – hạt vật chất (phương trình Lagrăng (1)–(2)) và so với nghiệm giải tích.

(1) I*θ'' + ΔE*sinθ - g*x*cosθ = 0
(2) m*x'' + k*x = g*sinθ

Chỉ dùng thư viện chuẩn Python (RK4 bước cố định), chạy: python3 giai_so_he_e_kep.py
"""
import math


def rhs(s, I, m, k, g, dE):
    th, w, x, v = s
    return (w, (-dE * math.sin(th) + g * x * math.cos(th)) / I,
            v, (-k * x + g * math.sin(th)) / m)


def rk4(s0, p, t_end, dt):
    """Tích phân RK4, trả về danh sách (t, θ, θ', x, x')."""
    out, s, t = [(0.0, *s0)], list(s0), 0.0
    for _ in range(int(round(t_end / dt))):
        k1 = rhs(s, *p)
        k2 = rhs([a + dt / 2 * b for a, b in zip(s, k1)], *p)
        k3 = rhs([a + dt / 2 * b for a, b in zip(s, k2)], *p)
        k4 = rhs([a + dt * b for a, b in zip(s, k3)], *p)
        s = [a + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4)]
        t += dt
        out.append((t, *s))
    return out


def energy(r, I, m, k, g, dE):
    _, th, w, x, v = r
    return 0.5 * I * w**2 + 0.5 * m * v**2 - dE * math.cos(th) + 0.5 * k * x**2 - g * x * math.sin(th)


def energy_drift(tr, p):
    e0 = energy(tr[0], *p)
    return max(abs(energy(r, *p) - e0) for r in tr) / abs(e0)


def mean_rotation_speed(E, I, dE, n=200000):
    """Ω̄ = 2π / ∮ dθ/θ', với θ' = √(2(E + ΔE cosθ)/I) (tích phân năng lượng của (1) khi g = 0)."""
    h = 2 * math.pi / n
    T = sum(h / math.sqrt(2 * (E + dE * math.cos((i + 0.5) * h)) / I) for i in range(n))
    return 2 * math.pi / T


def case_a():
    I, m, k, g, dE = 1.0, 1.0, 4.0, 0.1, 0.5
    w_start = 5.0                            # θ'(0), tại θ=0 (đỉnh tốc độ)
    w0 = math.sqrt(k / m)
    E_rot = 0.5 * I * w_start**2 - dE        # năng lượng của E_kép (bỏ qua g·x, rất nhỏ)
    Om = mean_rotation_speed(E_rot, I, dE)   # tốc độ quay trung bình thực
    C = g / (m * (w0**2 - Om**2))            # biên độ nghiệm riêng
    P = g * Om / (w0**2 - Om**2)             # biên độ xung lượng cưỡng bức
    # Chọn A = 0: x(0) = 0, x'(0) = C*Ω̄  → chỉ còn nghiệm cưỡng bức
    p = (I, m, k, g, dE)
    tr = rk4((0.0, w_start, 0.0, C * Om), p, 20.0, 1e-3)
    om_num = tr[-1][1] / tr[-1][0]
    err_x = max(abs(x - C * math.sin(Om * t)) for t, _, _, x, _ in tr)
    err_th = max(abs(th - Om * t) for t, th, _, _, _ in tr)
    wmin, wmax = min(r[2] for r in tr), max(r[2] for r in tr)
    zeros = [tr[i][0] for i in range(1, len(tr)) if tr[i - 1][4] * tr[i][4] < 0][:4]
    zeros_th = [(math.pi / 2 + n * math.pi) / Om for n in range(4)]
    print("=== Trường hợp A: quay gần đều, không cộng hưởng ===")
    print(f"I={I}, m={m}, k={k}, g={g}, ΔE={dE}, θ(0)=0, θ'(0)={w_start}, x(0)=0, x'(0)={C*Om:.6f}")
    print(f"ω0 = √(k/m) = {w0:.6f}")
    print(f"E_rot = ½Iθ'(0)² − ΔE = {E_rot:.4f};  E_rot/ΔE = {E_rot/dE:.1f} (≫ 1; vượt phân cách 2ΔE = {2*dE})")
    print(f"√(2E/I) = {math.sqrt(2*E_rot/I):.6f};  Ω̄ chính xác (tích phân chu kỳ) = {Om:.6f};  Ω̄ đo = θ(20)/20 = {om_num:.6f}")
    print(f"Ước lượng δΩ/Ω ~ ΔE/(IΩ²) = {dE/(I*Om**2):.4f}; đo được θ' ∈ [{wmin:.5f}, {wmax:.5f}]"
          f" → δΩ/Ω = {(wmax-wmin)/2/Om:.4f}")
    print(f"C = g/(m(ω0²−Ω̄²)) = {C:.7f}   |  biên độ p cưỡng bức = gΩ̄/(ω0²−Ω̄²) = {P:.7f}")
    print(f"max|θ_số − Ω̄t| trên [0,20] = {err_th:.5f} rad (chỉ là gợn tuần hoàn, không trôi)")
    print(f"max|x_số − C·sinΩ̄t| = {err_x:.2e}  (so với |C| = {abs(C):.2e}, sai số tương đối {err_x/abs(C):.2%})")
    print("t tại p=0 (số):       ", ", ".join(f"{z:.4f}" for z in zeros))
    print("t tại p=0 (giải tích):", ", ".join(f"{z:.4f}" for z in zeros_th))
    print(f"Trôi năng lượng tương đối max = {energy_drift(tr, p):.2e}")
    print(f"⟨p⟩ trung bình trên [0,20] = {sum(r[4] for r in tr)/len(tr)*m:.2e}  (biên độ p = {abs(P):.2e})\n")


def case_a_resonance():
    I, m, k, g, dE = 100.0, 1.0, 4.0, 0.1, 0.01
    w0 = math.sqrt(k / m)
    Om = w0
    p = (I, m, k, g, dE)
    # A = 0: x(0)=0, x'(0)=d/dt[-g/(2mω0) t cos ω0 t]_{t=0} = -g/(2mω0)
    tr = rk4((0.0, Om, 0.0, -g / (2 * m * w0)), p, 20.0, 1e-3)
    print("=== Trường hợp A(b): cộng hưởng Ω = ω0 ===")
    print(f"I={I}, m={m}, k={k}, g={g}, ΔE={dE}, Ω=ω0={w0}")
    print(f"Tốc độ tăng biên độ g/(2mω0) = {g/(2*m*w0):.5f} /đơn vị t")
    # so tại các đỉnh cos ω0t = ±1, tức t_n = nπ/ω0 (lấy đỉnh gần nhất ≤ T)
    for T in (5, 10, 15, 20):
        tn = math.floor(T * w0 / math.pi) * math.pi / w0
        amp = abs(min(tr, key=lambda r: abs(r[0] - tn))[3])
        print(f"  đỉnh t={tn:7.4f}: |x| số = {amp:.5f}   giải tích g·t/(2mω0) = {g*tn/(2*m*w0):.5f}")
    err = max(abs(x - (-g / (2 * m * w0)) * t * math.cos(w0 * t)) for t, _, _, x, _ in tr)
    print(f"max|x_số − x_cộng hưởng| trên [0,20] = {err:.2e}")
    print(f"Trôi năng lượng tương đối max = {energy_drift(tr, p):.2e}\n")


def case_b(g):
    I, m, k, dE = 1.0, 1.0, 1.0, 1.0
    wth2, w02 = dE / I, k / m
    disc = math.sqrt((wth2 - w02)**2 + 4 * g**2 / (I * m))
    wp2, wm2 = 0.5 * (wth2 + w02 + disc), 0.5 * (wth2 + w02 - disc)
    stable = dE * k > g**2
    print(f"=== Trường hợp B: dao động nhỏ, g={g} ===")
    print(f"I={I}, m={m}, k={k}, ΔE={dE} → ωθ = {math.sqrt(wth2):.4f}, ω0 = {math.sqrt(w02):.4f}")
    print(f"ΔE·k = {dE*k}, g² = {g**2:.4f} → {'ỔN ĐỊNH' if stable else 'KHÔNG ổn định'}")
    print(f"ω+² = {wp2:.6f}, ω−² = {wm2:.6f}")
    p = (I, m, k, g, dE)
    if stable:
        wp, wm = math.sqrt(wp2), math.sqrt(wm2)
        Tb = 2 * math.pi / (wp - wm)
        print(f"ω+ = {wp:.6f}, ω− = {wm:.6f}, chu kỳ phách 2π/(ω+−ω−) = {Tb:.4f}")
        tr = rk4((0.01, 0.0, 0.0, 0.0), p, Tb, 1e-3)
        # năng lượng dao động riêng của x: ½mv² + ½kx²
        ex = [(r[0], 0.5 * m * r[4]**2 + 0.5 * k * r[3]**2) for r in tr]
        tmax, emax = max(ex, key=lambda z: z[1])
        e_tot = 0.5 * dE * 0.01**2
        print(f"θ(0)=0.01, x(0)=0: năng lượng x cực đại tại t = {tmax:.3f} (dự đoán T_phách/2 = {Tb/2:.3f}),"
              f" chiếm {emax/e_tot:.1%} năng lượng ban đầu")
        print(f"Trôi năng lượng tương đối max = {energy_drift(tr, p):.2e}\n")
    else:
        lam = math.sqrt(-wm2)
        tr = rk4((1e-6, 0.0, 0.0, 0.0), p, 20.0, 1e-3)
        r5 = next(r for r in tr if r[0] >= 10.0)
        r10 = next(r for r in tr if r[0] >= 20.0 - 1e-9)
        lam_num = math.log(abs(r10[1]) / abs(r5[1])) / 10.0
        print(f"Tốc độ tăng giải tích λ = √(−ω−²) = {lam:.5f}; đo từ θ(t) trên [10,20]: {lam_num:.5f}")
        th_end = max(abs(r[1]) for r in tr)
        print(f"θ(0)=1e-6 → |θ| max trên [0,20] = {th_end:.3f} rad (hệ rời khỏi (0,0))\n")


if __name__ == "__main__":
    case_a()
    case_a_resonance()
    case_b(0.2)
    case_b(1.2)
