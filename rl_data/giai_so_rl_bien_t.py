"""Học tăng cường (RL) có biến thời gian t của Học thuyết E cho hệ E_kép (HocThuyetE (8).pdf).

Ba cách đưa t vào phương trình RL, mỗi cách giải số riêng:
  B  Cầu Euler–Lagrange ↔ Bellman: cực trị ∫e^{−ρt} L dt cho phương trình EL có thêm số hạng −ρ ∂L/∂q̇ (dạng
     Caldirola–Kanai với γ = −ρ). Kiểm bằng tích phân biến phân rời rạc (DEL) so với RK4, đo bậc hội tụ và dE/dt = 2ρT.
  D  t liên tục (đồng hồ Page–Wootters, nhịp τ = T0/25): điều khiển tối ưu ngẫu nhiên horizon hữu hạn trên đa tạp
     ràng buộc (pha E_past). Bellman V_k(s) = max_u { r(s,u) + E V_{k+1}(s') }, gãy |λ| ≥ R ⟹ hấp thụ, phạt −P.
     So với chính sách dừng có chiết khấu γ = e^{−τ/τ_dw} (τ_dw = thời gian lưu Poisson của Học thuyết E).
     Đánh giá cả hai bằng Monte Carlo trên SDE gốc (bước nhỏ, kiểm λ từng bước).
  S  t bán-Markov (SMDP): hành động = tỉ số lái f chọn ở mỗi lần tái khóa; thời gian lưu τ_c và năng lượng Carnot
     lấy từ mô phỏng Langevin đầy đủ (rl_data/chu_trinh_f.json). Bellman có chiết khấu e^{−ρτ}, horizon hữu hạn theo t,
     Q-learning bán-Markov không mô hình, và dạng đóng E[e^{−ρτ}] từ luật T0 + Exp(μ) và Geom(p) của bài báo.

Đơn vị: ħω0 = 1 (10 MeV), ω0 = 1, ℓ0 = 1, T0 = 2π. Thông số Bảng I của bài báo.
Chạy: python3 giai_so_rl_bien_t.py [B|D|S] …   (B < 5 s, D ≈ 1–2 phút, S ≈ 1 phút; S cần chạy giai_so_rl_thu_du_lieu.py trước)
"""
import json
import math
import sys

import numpy as np

T0 = 2 * math.pi
I, M, K, DE, KAP, G, R = 1.0, 1.0, 1.0, 1.0, 1.0 / 3.0, 0.1, 3.0
GTH, GX, KT = 0.1, 0.70, 0.1            # γθ = 1/τ, γx, k_BT = 1 MeV
TAU_DW = 5.25 * T0                      # thời gian lưu dư trung bình (bài báo: 4.5–6.0 T0)
TICK = T0 / 25                          # nhịp đồng hồ xung (Page–Wootters)
SN = 0.8                                # năng lượng tách nucleon S_n = 8 MeV


# ======================================================================= B. cầu EL ↔ Bellman
def el_rhs(s, rho):
    """EL chiết khấu của L hệ E_kép: Iθ̈ = ρIθ̇ − ΔE sinθ + gx cosθ,  mẍ = ρmẋ − kx + g sinθ."""
    th, w, x, v = s
    return np.array([w, rho * w + (-DE * math.sin(th) + G * x * math.cos(th)) / I,
                     v, rho * v + (-K * x + G * math.sin(th)) / M])


def rk4(f, s, h):
    k1 = f(s); k2 = f(s + h / 2 * k1); k3 = f(s + h / 2 * k2); k4 = f(s + h * k3)
    return s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def grad_pot(q):
    th, x = q
    return np.array([DE * math.sin(th) - G * x * math.cos(th), K * x - G * math.sin(th)])


def del_solve(q0, qd0, rho, h, T):
    """Tích phân biến phân: S_d = Σ_k e^{−ρ t_k} h [½ M_q ((q_{k+1}−q_k)/h)² − Π(q_k)].
    ∂S_d/∂q_k = 0 ⟹ e^{−ρt_{k−1}} M(q_k−q_{k−1})/h − e^{−ρt_k} M(q_{k+1}−q_k)/h − e^{−ρt_k} h ∇Π(q_k) = 0 (hiển)."""
    Mq = np.array([I, M])
    n = int(round(T / h))
    q_prev = np.array(q0, float)
    # khởi động bậc hai: q1 = q0 + h q̇0 + ½h² q̈0
    acc0 = el_rhs(np.array([q0[0], qd0[0], q0[1], qd0[1]]), rho)[[1, 3]]
    q = q_prev + h * np.array(qd0) + 0.5 * h * h * acc0
    e = math.exp(rho * h)
    for _ in range(n - 1):
        q_next = q + e * (q - q_prev) - h * h * grad_pot(q) / Mq
        q_prev, q = q, q_next
    return q


def part_b():
    print("=== B. Cực trị ∫e^{−ρt}L dt ⟺ EL + (−ρ ∂L/∂q̇): DEL rời rạc so với RK4 ===")
    q0, qd0, T = (0.3, 0.0), (0.0, 0.2), 4 * T0
    for rho in (0.0, 0.05, 0.10):
        s = np.array([q0[0], qd0[0], q0[1], qd0[1]])
        hr = T / 40000
        for _ in range(40000):
            s = rk4(lambda z: el_rhs(z, rho), s, hr)
        ref = s[[0, 2]]
        errs = []
        for n in (1250, 2500, 5000):              # h = T/n đúng điểm cuối
            h = T / n
            errs.append(np.max(np.abs(del_solve(q0, qd0, rho, h, T) - ref)))
        order = [math.log2(errs[i] / errs[i + 1]) for i in range(2)]
        print(f"ρ={rho:.2f}: |q_DEL − q_RK4| tại t=4T0, h=T/1250,/2500,/5000: "
              + " ".join(f"{e:.2e}" for e in errs) + f"  bậc {order[0]:.2f}, {order[1]:.2f}")
    # dE/dt = 2ρT  (động năng T) và VdP nhỏ-tốc-độ ≡ chiết khấu ρ = γθ/I trên phương trình θ
    rho, h = 0.10, 1e-3
    s = np.array([0.3, 0.0, 0.0, 0.2])
    E = lambda z: 0.5 * I * z[1]**2 + 0.5 * M * z[3]**2 - DE * math.cos(z[0]) + 0.5 * K * z[2]**2 - G * z[2] * math.sin(z[0])
    Tk = lambda z: 0.5 * I * z[1]**2 + 0.5 * M * z[3]**2
    E0, integ = E(s), 0.0
    for _ in range(int(round(2 * T0 / h))):
        s1 = rk4(lambda z: el_rhs(z, rho), s, h)
        integ += h * (Tk(s) + Tk(s1))           # ∫2ρT dt theo hình thang (ρ nhân sau)
        s = s1
    print(f"ρ=0.10, t∈[0,2T0]: ΔE = {E(s) - E0:.6f}, ∫2ρT dt = {rho * integ:.6f}"
          f"  (lệch tương đối {abs(E(s) - E0 - rho * integ) / abs(E(s) - E0):.1e})")
    print(f"VdP: Q_θ = γθ(1−θ̇²/Ω0²)θ̇ ⟹ ρ_eff(θ̇) = γθ(1−θ̇²/Ω0²)/I; θ̇→0: ρ_eff = {GTH / I:.3f} ω0, θ̇ = Ω0: ρ_eff = 0")
    print(f"Chiết khấu từ thời gian lưu Poisson: γ_RL = e^(−τ/τ_dw) = {math.exp(-TICK / TAU_DW):.5f} mỗi nhịp T0/25; "
          f"horizon hiệu dụng 1/(1−γ) = {1 / (1 - math.exp(-TICK / TAU_DW)):.1f} nhịp = {TAU_DW / T0:.2f} T0")


# ======================================================================= D. HJB/Bellman có t trên đa tạp ràng buộc
def locked_acc(th, w, u):
    """θ̈ và λ trên đa tạp x = κ sinθ, có mô-men điều khiển u, cản γθ, γx (Rayleigh chiếu lên θ)."""
    s, c = np.sin(th), np.cos(th)
    Mth = I + M * KAP**2 * c * c
    a = (u - GTH * w - GX * KAP**2 * c * c * w + M * KAP**2 * s * c * w * w
         - DE * s - (K * KAP**2 - 2 * G * KAP) * s * c) / Mth
    xdd = KAP * (c * a - s * w * w)
    lam = M * xdd + K * KAP * s - G * s + GX * KAP * c * w
    return a, lam


def crit_speed(u=0.0, gam=False):
    """θ̇ nhỏ nhất để max_θ |λ| ≥ R (không cản nếu gam=False) — kiểm với 3.205 ω0 của bài báo."""
    global GTH, GX
    keep = (GTH, GX)
    if not gam:
        GTH, GX = 0.0, 0.0
    th = np.linspace(-math.pi, math.pi, 20001)
    lo, hi = 0.0, 5.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        _, lam = locked_acc(th, mid, u)
        lo, hi = (lo, mid) if np.max(np.abs(lam)) >= R else (mid, hi)
    GTH, GX = keep
    return hi


NTH, NW, WLO, WHI = 72, 169, -0.5, 3.7
ACTIONS = np.array([-0.2, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6])
GH_X = np.array([0.0, 1.355626179974266, -1.355626179974266, 2.856970013872806, -2.856970013872806])
GH_W = np.array([0.533333333333333, 0.222075922005613, 0.222075922005613, 0.011257411327721, 0.011257411327721])
NSUB = 5


def tick_step(th, w, u):
    """Một nhịp τ (RK4 NSUB bước con) tất định; trả về trạng thái cuối, Δθ, công điều khiển, cờ gãy dọc đường."""
    h = TICK / NSUB
    broke = np.zeros(np.shape(th), bool)
    work = np.zeros(np.shape(th))
    th0 = th
    for _ in range(NSUB):
        a1, l1 = locked_acc(th, w, u)
        broke |= np.abs(l1) >= R
        k1 = (w, a1)
        a2, _ = locked_acc(th + h / 2 * k1[0], w + h / 2 * k1[1], u)
        k2 = (w + h / 2 * a1, a2)
        a3, _ = locked_acc(th + h / 2 * k2[0], w + h / 2 * k2[1], u)
        k3 = (w + h / 2 * a2, a3)
        a4, _ = locked_acc(th + h * k3[0], w + h * k3[1], u)
        k4 = (w + h * a3, a4)
        dth = h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        work += u * dth
        th = th + dth
        w = w + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    _, l_end = locked_acc(th, w, u)
    broke |= np.abs(l_end) >= R
    return th, w, th - th0, work, broke


def noise_sigma(th):
    """Độ lệch chuẩn cú hích θ̇ sau một nhịp: (σθ² + κ²cos²θ σx²) τ / M(θ)²."""
    c = np.cos(th)
    Mth = I + M * KAP**2 * c * c
    return np.sqrt((2 * GTH * KT + KAP**2 * c * c * 2 * GX * KT) * TICK) / Mth


def build_model(c_work):
    """Mô hình chuyển cho mọi (trạng thái lưới, hành động, nút Gauss–Hermite j).
    Thứ tự trong một nhịp: cú hích nhiễu θ̇ → θ̇ + x_j σ(θ) TRƯỚC, rồi dòng tất định có kiểm |λ| ≥ R
    (nhiễu đặt sau sẽ bỏ sót các lần vượt ngưỡng do nhiễu trong nhịp — đã đo: lưới 0.9% so với SDE 98%).
    Nội suy song tuyến (θ tuần hoàn) cho V tại trạng thái cuối."""
    thg = np.linspace(-math.pi, math.pi, NTH, endpoint=False)
    wg = np.linspace(WLO, WHI, NW)
    TH, W = np.meshgrid(thg, wg, indexing='ij')
    TH, W = TH.ravel(), W.ravel()
    ns, na, nj = TH.size, ACTIONS.size, GH_X.size
    idx = np.zeros((na, ns, nj, 4), np.int64)
    wt = np.zeros((na, ns, nj, 4))
    rew = np.zeros((na, ns, nj))
    brk = np.zeros((na, ns, nj), bool)
    dth_g, dw_g = 2 * math.pi / NTH, (WHI - WLO) / (NW - 1)
    sig = noise_sigma(TH)
    for ia, u in enumerate(ACTIONS):
        for j in range(nj):
            th1, w1, dth, work, broke = tick_step(TH, W + GH_X[j] * sig, u)
            rew[ia, :, j] = dth / (2 * math.pi) - c_work * np.maximum(work, 0.0)
            brk[ia, :, j] = broke
            w1 = np.clip(w1, WLO, WHI)
            ft = (np.mod(th1 + math.pi, 2 * math.pi)) / dth_g
            i0 = np.floor(ft).astype(np.int64) % NTH; i1 = (i0 + 1) % NTH; at = ft - np.floor(ft)
            fw = (w1 - WLO) / dw_g
            j0 = np.clip(np.floor(fw).astype(np.int64), 0, NW - 2); aw = fw - j0
            idx[ia, :, j] = np.stack([i0 * NW + j0, i0 * NW + j0 + 1, i1 * NW + j0, i1 * NW + j0 + 1], -1)
            wt[ia, :, j] = np.stack([(1 - at) * (1 - aw), (1 - at) * aw, at * (1 - aw), at * aw], -1)
    return dict(TH=TH, W=W, idx=idx, wt=wt, rew=rew, brk=brk, thg=thg, wg=wg)


def backup(mdl, V, P, gamma=1.0):
    """Q(s,a) = Σ_j w_j [gãy_j ? −P : r_j + γ V(s'_j)];  gãy ⟹ hấp thụ (kết thúc pha E_past)."""
    cont = (V[mdl['idx']] * mdl['wt']).sum(axis=3)
    return (GH_W * np.where(mdl['brk'], -P, mdl['rew'] + gamma * cont)).sum(axis=2)


def solve_finite(mdl, P, nk):
    """Quy nạp lùi theo đồng hồ k = 0..nk−1: chính sách phụ thuộc t."""
    V = np.zeros(mdl['TH'].size)
    pol = np.zeros((nk, V.size), np.int64)
    for k in range(nk - 1, -1, -1):
        Q = backup(mdl, V, P)
        pol[k] = Q.argmax(0)
        V = Q.max(0)
    return V, pol


def solve_discounted(mdl, P, gamma, tol=1e-9):
    V = np.zeros(mdl['TH'].size)
    for it in range(20000):
        Q = backup(mdl, V, P, gamma)
        Vn = Q.max(0)
        if np.max(np.abs(Vn - V)) < tol:
            break
        V = Vn
    return Q.argmax(0), it


def eval_policy_grid(mdl, P, pol_of_k, nk):
    """Đánh giá chính xác trên lưới (cùng mô hình) — trả V_0 và xác suất gãy trước hạn."""
    V = np.zeros(mdl['TH'].size)
    B = np.zeros(mdl['TH'].size)          # xác suất gãy trước hạn
    ar = np.arange(V.size)
    for k in range(nk - 1, -1, -1):
        a = pol_of_k(k)
        i, w = mdl['idx'][a, ar], mdl['wt'][a, ar]
        br = mdl['brk'][a, ar]
        V = (GH_W * np.where(br, -P, mdl['rew'][a, ar] + (V[i] * w).sum(axis=2))).sum(axis=1)
        B = (GH_W * np.where(br, 1.0, (B[i] * w).sum(axis=2))).sum(axis=1)
    return V, B


def mc_sde(mdl, pol_of_k, nk, s0, P, c_work, n=4000, sub=10, seed=7):
    """Monte Carlo trên SDE gốc: bước h = τ/sub, Euler–Maruyama, kiểm |λ| ≥ R mỗi bước, hành động giữ trong nhịp."""
    rng = np.random.default_rng(seed)
    th = np.full(n, s0[0]); w = np.full(n, s0[1])
    alive = np.ones(n, bool); ret = np.zeros(n); tbr = np.full(n, np.nan)
    h = TICK / sub
    dth_g, dw_g = 2 * math.pi / NTH, (WHI - WLO) / (NW - 1)
    for k in range(nk):
        it = np.round(np.mod(th + math.pi, 2 * math.pi) / dth_g).astype(np.int64) % NTH
        iw = np.clip(np.round((w - WLO) / dw_g).astype(np.int64), 0, NW - 1)
        u = ACTIONS[pol_of_k(k)[it * NW + iw]]
        for _ in range(sub):
            a, lam = locked_acc(th, w, u)
            newly = alive & (np.abs(lam) >= R)
            ret[newly] -= P; tbr[newly] = k * TICK; alive &= ~newly
            c = np.cos(th)
            Mth = I + M * KAP**2 * c * c
            dth = w * h
            ret[alive] += dth[alive] / (2 * math.pi) - c_work * np.maximum(u * dth, 0)[alive]
            kick = (math.sqrt(2 * GTH * KT * h) * rng.standard_normal(n)
                    + KAP * c * math.sqrt(2 * GX * KT * h) * rng.standard_normal(n)) / Mth
            th = np.where(alive, th + dth, th)
            w = np.where(alive, w + a * h + kick, w)
    return ret.mean(), ret.std() / math.sqrt(n), 1 - alive.mean(), np.nanmean(tbr) / T0 if (~alive).any() else float('nan')


def part_d(c_work=0.15):
    print("=== D. Bellman có t (đồng hồ τ = T0/25) trên đa tạp E_past, điều khiển mô-men u ===")
    from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params
    base, Rr = params(1.0, 3.0, I=1.0, g=0.1)
    print(f"kiểm λ: tiệm cận √(R/mκ) = {math.sqrt(R / (M * KAP)):.4f}; θ̇_gãy theo quỹ đạo (θ̇ tại θ=0, mã Part I) = "
          f"{omega_crit_general(base, Rr):.4f} (bài báo 3.205)")
    print(f"       ngưỡng tĩnh min{{θ̇ : max_θ|λ(θ,θ̇)| ≥ R}}: không cản {crit_speed():.4f}; có cản {crit_speed(gam=True):.4f};"
          f" u=0.4: {crit_speed(0.4, True):.4f}")
    mdl = build_model(c_work)
    gamma = math.exp(-TICK / TAU_DW)
    for P, Hk in ((2.0, 20), (10.0, 20)):
        nk = int(Hk * 25)
        V0, pol = solve_finite(mdl, P, nk)
        pol_s, it = solve_discounted(mdl, P, gamma)
        Vt, Bt = eval_policy_grid(mdl, P, lambda k: pol[k], nk)
        Vs, Bs = eval_policy_grid(mdl, P, lambda k: pol_s, nk)
        s0 = (0.0, 2.15)                               # θ̇ ngay sau tái khóa (bài báo: θ̇_end ≈ 2.15)
        i0 = int(np.argmin(np.abs(mdl['thg'] - s0[0]))) * NW + int(np.argmin(np.abs(mdl['wg'] - s0[1])))
        print(f"\n-- phạt gãy P = {P:g} vòng, horizon {Hk} T0 = {nk} nhịp, giá công c = {c_work}, "
              f"γ = {gamma:.5f} (VI hội tụ {it} vòng) --")
        print(f"lưới  s0=(0, 2.15): V_t-phụ-thuộc = {Vt[i0]:.3f}, P(gãy) = {Bt[i0]:.3f} | "
              f"V_dừng = {Vs[i0]:.3f}, P(gãy) = {Bs[i0]:.3f}")
        mt = mc_sde(mdl, lambda k: pol[k], nk, s0, P, c_work)
        ms = mc_sde(mdl, lambda k: pol_s, nk, s0, P, c_work)
        print(f"SDE   s0=(0, 2.15): V_t-phụ-thuộc = {mt[0]:.3f} ± {mt[1]:.3f}, P(gãy) = {mt[2]:.3f} | "
              f"V_dừng = {ms[0]:.3f} ± {ms[1]:.3f}, P(gãy) = {ms[2]:.3f}")
        # mép an toàn: θ̇ lớn nhất (tại θ = 0) mà chính sách còn chọn u > 0, theo thời gian còn lại
        ith = int(np.argmin(np.abs(mdl['thg'])))
        line = []
        for left in (20, 10, 5, 2, 1, 0.4):
            k = nk - int(round(left * 25))
            a = pol[k].reshape(NTH, NW)[ith]
            ws = mdl['wg'][ACTIONS[a] > 0]
            line.append(f"{left:g}T0→{ws.max() if ws.size else float('nan'):.2f}")
        a = pol_s.reshape(NTH, NW)[ith]
        ws = mdl['wg'][ACTIONS[a] > 0]
        print("θ̇ lớn nhất còn đẩy (u>0) tại θ=0 theo thời gian còn lại: " + ", ".join(line)
              + f" | dừng: {ws.max() if ws.size else float('nan'):.2f}")
        ubar = [float(np.mean(ACTIONS[pol[nk - int(round(left * 25))]].reshape(NTH, NW)[:, (mdl['wg'] > 1.5) & (mdl['wg'] < 3.0)]))
                for left in (20, 5, 1, 0.4)]
        print("⟨u⟩ trên dải θ̇∈(1.5,3) theo thời gian còn lại 20/5/1/0.4 T0: " + " ".join(f"{x:.3f}" for x in ubar))


# ======================================================================= S. SMDP trên chu trình ba pha
def load_cycles(path='rl_data/chu_trinh_f.json', warm=5 * T0):
    import giai_so_burst_g as BG
    data = json.load(open(path))
    out = {}
    for tr in data:
        ev = [tuple(e) for e in tr['events']]
        dwell, bursts = BG.analyse(ev, warm)
        # ghép mỗi đợt với pha lưu ngay trước nó (cùng thứ tự xuất hiện trong analyse)
        cyc = [(d + b['dur'], b['E'], b['n'], d) for d, b in zip(dwell, bursts)]
        o = out.setdefault(tr['f'], dict(cyc=[], W=0.0, T=0.0))
        o['cyc'] += cyc; o['W'] += tr['W']; o['T'] += tr['T']
    for f, o in out.items():
        o['PW'] = o['W'] / o['T']                     # công suất lái trung bình
        o['cyc'] = np.array(o['cyc'])
    return out


def cyc_reward(o, c_work, rho):
    """R chiết khấu của một chu trình: −c P_W ∫_0^τ e^{−ρt}dt + E_Carnot e^{−ρτ}."""
    tau, E = o['cyc'][:, 0], o['cyc'][:, 1]
    if rho == 0:
        return -c_work * o['PW'] * tau + E
    return -c_work * o['PW'] * (1 - np.exp(-rho * tau)) / rho + E * np.exp(-rho * tau)


def closed_form_disc(rho, mu, p, delta):
    """E[e^{−ρτ}], τ = T0 + Exp(μ) + N δ, N ~ Geom(p) trên {1,2,…}."""
    z = math.exp(-rho * delta)
    return math.exp(-rho * T0) / (1 + rho * mu) * p * z / (1 - (1 - p) * z)


def part_s(c_work=0.15):
    print("=== S. Bellman bán-Markov: hành động f, thời gian lưu τ_c từ Langevin đầy đủ ===")
    D = load_cycles()
    fs = sorted(D)
    print(" f    chu-trình  ⟨τ_c⟩/T0  ⟨dwell−T0⟩/T0  CV_dư   ⟨N⟩    ⟨E_C⟩    P(E≥Sn)  P_W     gain(c)")
    for f in fs:
        o = D[f]; c = o['cyc']
        res = c[:, 3] - T0
        gain = cyc_reward(o, c_work, 0).mean() / c[:, 0].mean()
        print(f"{f:.2f}  {len(c):5d}     {c[:, 0].mean() / T0:6.2f}    {res.mean() / T0:6.2f}       "
              f"{res.std() / res.mean():.2f}   {c[:, 2].mean():5.2f}  {c[:, 1].mean():.3f}   "
              f"{(c[:, 1] >= SN).mean():.3f}   {o['PW']:.4f}  {gain * T0:+.4f}/T0")
    # dạng đóng từ luật thời gian của Học thuyết E so với thực nghiệm
    if 1.05 in D:
        c = D[1.05]['cyc']
        mu = (c[:, 3] - T0).mean(); p = 1 / c[:, 2].mean(); delta = (c[:, 0] - c[:, 3]).sum() / c[:, 2].sum()
        print(f"\nf=1.05: μ = {mu / T0:.2f} T0, p = {p:.3f}, δ = {delta / T0:.3f} T0/cú va")
        for rho in (0.01, 0.03, 0.1):
            print(f"  ρ = {rho:.2f}: E[e^(−ρτ)] thực nghiệm = {np.exp(-rho * c[:, 0]).mean():.4f}, "
                  f"dạng đóng T0+Exp(μ)+Geom(p)·δ = {closed_form_disc(rho, mu, p, delta):.4f}")
    # chính sách dừng theo ρ (một trạng thái tái phát ⟹ V_f = E[R_f]/(1 − E[e^{−ρτ_f}]))
    print("\nf* tối ưu theo ρ (dừng, vô hạn):")
    for rho in (0.0, 0.003, 0.01, 0.03, 0.1, 0.3):
        if rho == 0:
            vals = {f: cyc_reward(D[f], c_work, 0).mean() / D[f]['cyc'][:, 0].mean() for f in fs}
            tag = "gain/T0"
            vals = {f: v * T0 for f, v in vals.items()}
        else:
            vals = {f: cyc_reward(D[f], c_work, rho).mean() / (1 - np.exp(-rho * D[f]['cyc'][:, 0]).mean()) for f in fs}
            tag = "V"
        best = max(vals, key=vals.get)
        print(f"  ρ = {rho:<5}: f* = {best:.2f}  ({tag}: " + ", ".join(f"{f:.2f}:{v:+.3f}" for f, v in vals.items()) + ")")
    # horizon hữu hạn theo t: V(t) = max_f E[R_f(τ ∧ (H−t)) + V(t+τ)]
    H = 60 * T0; dt = TICK; nt = int(round(H / dt))
    V = np.zeros(nt + 1); polt = np.zeros(nt + 1, int)
    for i in range(nt - 1, -1, -1):
        left = H - i * dt
        best, arg = -1e18, 0
        for a, f in enumerate(fs):
            o = D[f]; tau, E = o['cyc'][:, 0], o['cyc'][:, 1]
            done = tau <= left
            j = np.minimum(i + np.round(tau / dt).astype(int), nt)
            val = np.where(done, -c_work * o['PW'] * tau + E + V[j], -c_work * o['PW'] * left).mean()
            if val > best:
                best, arg = val, a
        V[i], polt[i] = best, arg
    seg, prev = [], None
    for i in range(nt + 1):
        if polt[i] != prev:
            seg.append(f"t={i * dt / T0:.1f}T0→f={fs[polt[i]]:.2f}"); prev = polt[i]
    print(f"\nHorizon 60 T0, ρ = 0: V(0) = {V[0]:.3f}; chính sách theo t: " + "; ".join(seg))
    # Q-learning bán-Markov không mô hình (môi trường = lấy mẫu chu trình thực nghiệm)
    rng = np.random.default_rng(3)
    for rho in (0.01, 0.1):
        Q = np.zeros(len(fs)); n = np.zeros(len(fs))
        for ep in range(200000):
            a = rng.integers(len(fs)) if rng.random() < 0.2 else int(Q.argmax())
            o = D[fs[a]]; c = o['cyc'][rng.integers(len(o['cyc']))]
            tau, E = c[0], c[1]
            Rr = -c_work * o['PW'] * (1 - math.exp(-rho * tau)) / rho + E * math.exp(-rho * tau)
            n[a] += 1
            Q[a] += (1 / n[a]**0.6) * (Rr + math.exp(-rho * tau) * Q.max() - Q[a])
        exact = {f: cyc_reward(D[f], c_work, rho).mean() / (1 - np.exp(-rho * D[f]['cyc'][:, 0]).mean()) for f in fs}
        print(f"Q-learning SMDP ρ={rho}: f* học = {fs[int(Q.argmax())]:.2f}, Q* = {Q.max():+.3f};"
              f" chính xác f* = {max(exact, key=exact.get):.2f}, V* = {max(exact.values()):+.3f}")


if __name__ == '__main__':
    for p in (sys.argv[1:] or ['B', 'D', 'S']):
        {'B': part_b, 'D': part_d, 'S': part_s}[p.upper()]()
        print()
