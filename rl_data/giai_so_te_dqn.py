"""TE-DQN: Deep Q-Network có biến thời gian t của Học thuyết E, giải số trên hệ E_kép đủ ba pha.

Môi trường (EKepEnv, vector hoá numpy):
  E_past    : đa tạp ràng buộc x = κ sinθ, phương trình rút gọn (6) + cản γθ, γx + Langevin k_BT, mô-men điều khiển u.
              Nhịp quyết định τ = T0/25 (đồng hồ xung Page–Wootters), 10 bước con Euler–Maruyama, kiểm |λ| ≥ R mỗi bước con.
  E_present : khi gãy, đợt Carnot N ~ Geom(p = 0.195) cú va, mỗi cú cách δ = 0.36 T0 (luật thời gian của bài báo);
              trong đợt động cơ mất đồng bộ — không tính vòng quay.
  E_future  : tái khóa với θ̇ = 2.15 (θ̇_end của bài báo), θ ~ U(−π, π); đồng hồ đã nhảy (1 + 9N) nhịp.
  Phần thưởng mỗi nhịp: Δθ/2π − c·max(u Δθ, 0), c = 0.15. Horizon H = 20 T0 = 500 nhịp đồng hồ.

Chuẩn chính xác: quy hoạch động trên lưới (mô hình nhiễu-trước của giai_so_rl_bien_t) có thêm nhánh đợt Carnot.

Thuật toán (cùng mạng 5/6 → 128 → 128 → 7, Double DQN, replay, mạng đích mềm, Adam, Huber — viết tay bằng numpy):
  V  DQN thường      : không có t, γ = 0.99241 cố định mỗi chuyển, hết horizon coi là kết thúc.
  S  DQN bán-Markov  : không có t, γ^{Δt/τ} theo thời gian trôi thật Δt (gồm cả đợt), hết horizon là cắt (vẫn bootstrap).
  T  TE-DQN          : có t_còn/H trong trạng thái, ρ = 0 (đúng mục tiêu horizon hữu hạn), đồng hồ nhảy theo Δt.
  G  TE-DQN-γ       : có t_còn/H, chiết khấu e^{−ρΔt} (ρ = 1/τ_dw), H là kết thúc thật.
  H  TE-DQN + HJ     : T + định hình thế năng Hamilton–Jacobi Φ(s,t) = ζ θ̇ (H − t)/2π (bất biến chính sách, Φ = 0 ở cuối).
Chạy: python3 giai_so_te_dqn.py DP | train V|S|T|H seed | eval     (kết quả → rl_data/te_dqn_*.json)
"""
import json
import math
import os
import sys
import time

import numpy as np

import giai_so_rl_bien_t as G

T0, TICK = G.T0, G.TICK
K = 500                                  # horizon = 20 T0 theo đồng hồ
P_GEOM, D_TICKS = 0.195, 9               # Geom(p) và δ = 0.36 T0 = 9 nhịp
W_RELOCK = 2.15
C_WORK = 0.15
RHO = 1.0 / G.TAU_DW                     # chiết khấu từ thời gian lưu Poisson
GAMMA_TICK = math.exp(-RHO * TICK)
ACTIONS = G.ACTIONS
NA = ACTIONS.size
NMAX = 80                                # cắt chuỗi Geom: (1−p)^80 ≈ 3e−8
PN = P_GEOM * (1 - P_GEOM) ** np.arange(NMAX)          # P(N = n), n = 1..NMAX


# ======================================================================= môi trường
class EKepEnv:
    def __init__(self, n, seed, sub=10):
        self.n, self.sub, self.rng = n, sub, np.random.default_rng(seed)
        self.reset_all()

    def _relock(self, mask):
        self.th[mask] = self.rng.uniform(-math.pi, math.pi, mask.sum())
        self.w[mask] = W_RELOCK

    def reset_all(self):
        self.th = np.zeros(self.n); self.w = np.zeros(self.n); self.k = np.zeros(self.n, np.int64)
        self._relock(np.ones(self.n, bool))
        self.nbreak = np.zeros(self.n, np.int64)

    def step(self, a):
        """Một nhịp đồng hồ. Trả r, Δk (số nhịp trôi, gồm cả đợt), cờ gãy, cờ hết horizon."""
        u = ACTIONS[a] if a.dtype.kind == 'i' else a
        h = TICK / self.sub
        th, w = self.th.copy(), self.w.copy()
        alive = np.ones(self.n, bool); r = np.zeros(self.n)
        for _ in range(self.sub):
            acc, lam = G.locked_acc(th, w, u)
            alive &= np.abs(lam) < G.R
            c = np.cos(th)
            Mth = G.I + G.M * G.KAP**2 * c * c
            dth = w * h
            r += np.where(alive, dth / (2 * math.pi) - C_WORK * np.maximum(u * dth, 0), 0.0)
            kick = (math.sqrt(2 * G.GTH * G.KT * h) * self.rng.standard_normal(self.n)
                    + G.KAP * c * math.sqrt(2 * G.GX * G.KT * h) * self.rng.standard_normal(self.n)) / Mth
            th = np.where(alive, th + dth, th)
            w = np.where(alive, w + acc * h + kick, w)
        broke = ~alive
        dk = np.ones(self.n, np.int64)
        if broke.any():                                       # E_present: đợt Carnot, rồi E_future: tái khóa
            N = self.rng.geometric(P_GEOM, broke.sum())
            dk[broke] += D_TICKS * N
            self.nbreak[broke] += 1
        self.th, self.w = np.mod(th + math.pi, 2 * math.pi) - math.pi, w
        self._relock(broke)
        self.k = self.k + dk
        end = self.k >= K
        return r, dk, broke, end

    def lam(self):
        return G.locked_acc(self.th, self.w, 0.0)[1]


def features(th, w, lam, k, with_t):
    f = [np.cos(th), np.sin(th), w / 3.0, lam / G.R]
    if with_t:
        f.append(np.maximum(K - k, 0) / K)
    return np.stack(f, -1)


# ======================================================================= quy hoạch động chính xác (chuẩn)
def dp_solve(mdl=None):
    """Horizon hữu hạn có nhánh đợt Carnot: Q_k = Σ_j w_j [r_j + (gãy ? B(k+1) : V_{k+1}(s'_j))],
    B(m) = Σ_n P(N=n) V̄_relock(m + 9n), V̄_relock(k) = ⟨V_k(θ, 2.15)⟩_θ, V_{k≥K} = 0."""
    mdl = mdl or G.build_model(C_WORK)
    jr = int(np.argmin(np.abs(mdl['wg'] - W_RELOCK)))
    ns = mdl['TH'].size
    Vrel = np.zeros(K + D_TICKS * NMAX + 2)
    V = np.zeros(ns)
    pol = np.zeros((K, ns), np.int8)
    offs = D_TICKS * np.arange(1, NMAX + 1)
    for k in range(K - 1, -1, -1):
        B = float((PN * Vrel[k + 1 + offs]).sum())
        cont = (V[mdl['idx']] * mdl['wt']).sum(axis=3)
        Q = (G.GH_W * (mdl['rew'] + np.where(mdl['brk'], B, cont))).sum(axis=2)
        pol[k] = Q.argmax(0)
        V = Q.max(0)
        Vrel[k] = V.reshape(G.NTH, G.NW)[:, jr].mean()
    return V, pol, mdl, Vrel[0]


def dp_stationary(mdl, tol=1e-7):
    """Chính sách dừng có chiết khấu SMDP: chuyển thường γ, đợt e^{−ρ(1+9N)τ}."""
    jr = int(np.argmin(np.abs(mdl['wg'] - W_RELOCK)))
    disc_b = float((PN * np.exp(-RHO * TICK * (1 + D_TICKS * np.arange(1, NMAX + 1)))).sum())
    V = np.zeros(mdl['TH'].size)
    for it in range(30000):
        vr = V.reshape(G.NTH, G.NW)[:, jr].mean()
        cont = (V[mdl['idx']] * mdl['wt']).sum(axis=3)
        Q = (G.GH_W * (mdl['rew'] + np.where(mdl['brk'], disc_b * vr, GAMMA_TICK * cont))).sum(axis=2)
        Vn = Q.max(0)
        if np.max(np.abs(Vn - V)) < tol:
            break
        V = Vn
    return Q.argmax(0).astype(np.int8)


def grid_lookup(th, w):
    it = np.round(np.mod(th + math.pi, 2 * math.pi) / (2 * math.pi / G.NTH)).astype(np.int64) % G.NTH
    iw = np.clip(np.round((w - G.WLO) / ((G.WHI - G.WLO) / (G.NW - 1))).astype(np.int64), 0, G.NW - 1)
    return it * G.NW + iw


# ======================================================================= đánh giá chung (số ngẫu nhiên chung)
def evaluate(policy, n=2000, seed=12345):
    """policy(env) → mảng hành động (chỉ số) hoặc mô-men liên tục. Trả trung bình lợi tức, SE, số lần gãy/tập."""
    env = EKepEnv(n, seed)
    ret = np.zeros(n); done = np.zeros(n, bool)
    while not done.all():
        a = policy(env)
        r, dk, broke, end = env.step(a)
        ret += np.where(done, 0.0, r)
        done |= end
    return float(ret.mean()), float(ret.std() / math.sqrt(n)), float(env.nbreak.mean())


# ======================================================================= mạng nơ-ron numpy
class MLP:
    def __init__(self, din, dout, rng, hid=128):
        he = lambda a, b: rng.standard_normal((a, b)) * math.sqrt(2.0 / a)
        self.p = [he(din, hid), np.zeros(hid), he(hid, hid), np.zeros(hid), he(hid, dout) * 0.1, np.zeros(dout)]
        self.m = [np.zeros_like(x) for x in self.p]; self.v = [np.zeros_like(x) for x in self.p]; self.t = 0

    def forward(self, x, cache=False):
        W1, b1, W2, b2, W3, b3 = self.p
        z1 = x @ W1 + b1; a1 = np.maximum(z1, 0)
        z2 = a1 @ W2 + b2; a2 = np.maximum(z2, 0)
        out = a2 @ W3 + b3
        if cache:
            self._c = (x, z1, a1, z2, a2)
        return out

    def backward_adam(self, gout, lr):
        x, z1, a1, z2, a2 = self._c
        W1, b1, W2, b2, W3, b3 = self.p
        gW3 = a2.T @ gout; gb3 = gout.sum(0)
        g2 = (gout @ W3.T) * (z2 > 0)
        gW2 = a1.T @ g2; gb2 = g2.sum(0)
        g1 = (g2 @ W2.T) * (z1 > 0)
        gW1 = x.T @ g1; gb1 = g1.sum(0)
        grads = [gW1, gb1, gW2, gb2, gW3, gb3]
        gn = math.sqrt(sum(float((g * g).sum()) for g in grads))
        if gn > 10.0:                                        # cắt chuẩn gradient
            grads = [g * (10.0 / gn) for g in grads]
        self.t += 1
        for i, g in enumerate(grads):
            self.m[i] = 0.9 * self.m[i] + 0.1 * g
            self.v[i] = 0.999 * self.v[i] + 0.001 * g * g
            mh = self.m[i] / (1 - 0.9**self.t); vh = self.v[i] / (1 - 0.999**self.t)
            self.p[i] -= lr * mh / (np.sqrt(vh) + 1e-8)

    def soft_update(self, src, tau):
        for i in range(len(self.p)):
            self.p[i] = (1 - tau) * self.p[i] + tau * src.p[i]

    def copy_from(self, src):
        self.p = [x.copy() for x in src.p]


VARIANTS = {
    'V': dict(with_t=False, smdp=False, finite=False, hj=False, name='DQN thường'),
    'S': dict(with_t=False, smdp=True, finite=False, hj=False, name='DQN bán-Markov γ^Δt'),
    'T': dict(with_t=True, smdp=True, finite=True, hj=False, name='TE-DQN'),
    'H': dict(with_t=True, smdp=True, finite=True, hj=True, name='TE-DQN + Hamilton–Jacobi'),
    'G': dict(with_t=True, smdp=True, finite=False, hj=False, end_terminal=True, name='TE-DQN-γ (t + e^{−ρΔt})'),
}
ZETA = 0.6                                # ≈ (lợi tức DP 34.4/500 nhịp) / (θ̇ τ/2π tại θ̇ ≈ 2.8)


def phi(w, k):
    """Thế Hamilton–Jacobi: với rô-to tự do S = −h t + p θ (sách, (7.16)), giá trị còn lại ≈ θ̇ (H−t)/2π vòng."""
    return ZETA * w * np.maximum(K - k, 0) * TICK / (2 * math.pi)


def train(var, seed, iters=8000, nenv=32, batch=256, updates=2, lr=5e-4, tau_soft=0.005, buf=200000):
    cfg = VARIANTS[var]
    rng = np.random.default_rng(seed)
    din = 5 if cfg['with_t'] else 4
    net, tgt = MLP(din, NA, rng), MLP(din, NA, rng); tgt.copy_from(net)
    env = EKepEnv(nenv, seed + 99)
    S = np.zeros((buf, din)); A = np.zeros(buf, np.int64); Rw = np.zeros(buf); S2 = np.zeros((buf, din))
    Dsc = np.zeros(buf); ptr = size = 0
    t0 = time.time()
    s = features(env.th, env.w, env.lam(), env.k, cfg['with_t'])
    best, best_p = -1e9, None
    pol_of = lambda e: net.forward(features(e.th, e.w, e.lam(), e.k, cfg['with_t'])).argmax(1)
    for it in range(iters):
        eps = max(0.05, 1.0 - it / (0.5 * iters))
        lr_t = lr * (1.0 - 0.8 * it / iters)                           # giảm tuyến tính lr → 0.2 lr
        if it >= iters // 2 and (it % 2500 == 0 or it == iters - 1):    # kiểm định (seed riêng, không dùng tập đánh giá)
            vm = evaluate(pol_of, n=300, seed=777 + seed)[0]
            if vm > best:
                best, best_p = vm, [x.copy() for x in net.p]
        greedy = net.forward(s).argmax(1)
        a = np.where(rng.random(nenv) < eps, rng.integers(NA, size=nenv), greedy)
        w_old, k_old = env.w.copy(), env.k.copy()
        r, dk, broke, end = env.step(a)
        s2 = features(env.th, env.w, env.lam(), env.k, cfg['with_t'])
        if cfg['hj']:
            r = r + phi(env.w, env.k) * (~end) - phi(w_old, k_old)     # F = Φ(s′,t′) − Φ(s,t), Φ(cuối) = 0
        if cfg['finite']:
            disc = np.where(end, 0.0, 1.0)                              # ρ = 0, kết thúc thật ở H
        elif cfg['smdp'] and cfg.get('end_terminal'):
            disc = np.where(end, 0.0, np.exp(-RHO * TICK * dk))         # t trong trạng thái ⟹ H là kết thúc thật
        elif cfg['smdp']:
            disc = np.exp(-RHO * TICK * dk)                             # hết horizon = cắt, vẫn bootstrap
        else:
            disc = np.where(end, 0.0, GAMMA_TICK)                       # cách làm thường: γ cố định, cắt = kết thúc
        idx = (ptr + np.arange(nenv)) % buf
        S[idx], A[idx], Rw[idx], S2[idx], Dsc[idx] = s, a, r, s2, disc
        ptr = (ptr + nenv) % buf; size = min(size + nenv, buf)
        if end.any():                                                   # tập mới cho môi trường đã hết giờ
            env.k[end] = 0; env.nbreak[end] = 0
            env._relock(end)
            s2 = features(env.th, env.w, env.lam(), env.k, cfg['with_t'])
        s = s2
        if size >= 5000:
            for _ in range(updates):
                j = rng.integers(size, size=batch)
                a_star = net.forward(S2[j]).argmax(1)                   # Double DQN
                y = Rw[j] + Dsc[j] * tgt.forward(S2[j])[np.arange(batch), a_star]
                q = net.forward(S[j], cache=True)
                err = q[np.arange(batch), A[j]] - y
                g = np.zeros_like(q)
                g[np.arange(batch), A[j]] = np.clip(err, -1, 1) / batch     # Huber
                net.backward_adam(g, lr_t)
                tgt.soft_update(net, tau_soft)
    net.p = best_p
    m, se, nb = evaluate(pol_of)
    out = dict(var=var, name=cfg['name'], seed=seed, ret=m, se=se, breaks=nb, val=best, sec=time.time() - t0,
               iters=iters, transitions=iters * nenv)
    os.makedirs('rl_data', exist_ok=True)
    json.dump(out, open(f'rl_data/te_dqn_{var}_{seed}.json', 'w'))
    np.savez(f'rl_data/te_dqn_{var}_{seed}.npz', *net.p)
    print(json.dumps(out, ensure_ascii=False), flush=True)


def run_dp():
    t0 = time.time()
    V, pol, mdl, vrel0 = dp_solve()
    pol_s = dp_stationary(mdl)
    print(f"DP: V̄_relock(0) trên lưới = {vrel0:.3f}  ({time.time() - t0:.0f} s)")
    res = {}
    res['DP_t'] = evaluate(lambda e: pol[np.minimum(e.k, K - 1), grid_lookup(e.th, e.w)])
    res['DP_dung'] = evaluate(lambda e: pol_s[grid_lookup(e.th, e.w)])
    res['u0'] = evaluate(lambda e: np.zeros(e.n, np.int64) + 1)
    best = None
    for ia in range(NA):
        v = evaluate(lambda e, ia=ia: np.full(e.n, ia, np.int64), n=500)
        best = (v, ia) if best is None or v[0] > best[0][0] else best
    res[f'u_const_{ACTIONS[best[1]]:g}'] = evaluate(lambda e: np.full(e.n, best[1], np.int64))
    Om = 1.05 * 3.2048
    res['VdP_bai_bao'] = evaluate(lambda e: G.GTH * e.w * (1 - e.w**2 / Om**2) + G.GTH * e.w)
    for k_, v in res.items():
        print(f"{k_:16s} lợi tức = {v[0]:.3f} ± {v[1]:.3f}, số lần gãy/tập = {v[2]:.3f}")
    json.dump(dict(vrel0=vrel0, **{k_: list(v) for k_, v in res.items()}), open('rl_data/te_dqn_dp.json', 'w'))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'DP':
        run_dp()
    elif cmd == 'train':
        train(sys.argv[2], int(sys.argv[3]), iters=int(sys.argv[4]) if len(sys.argv) > 4 else 8000)
    elif cmd == 'eval':
        import glob
        rows = [json.load(open(f)) for f in sorted(glob.glob('rl_data/te_dqn_[VSTHG]_*.json'))]
        for r in rows:
            print(f"{r['name']:28s} seed {r['seed']}: {r['ret']:.3f} ± {r['se']:.3f}, gãy/tập {r['breaks']:.3f}, {r['sec']:.0f} s")
