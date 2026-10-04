"""Giải số bản thiết kế E-DQN (4 trụ cột) trên môi trường E_kép đủ ba pha của giai_so_te_dqn.

Bản thiết kế:
  Trụ 1  chiết khấu vật lý γ(Δt) = e^{−ρΔt}, ρ = 0.03 (≈ 1/τ_dw = 0.0303); tùy chọn chiết khấu Van der Pol
         ρ_eff(s) = (γθ/I)(1 − θ̇²/Ω0²)  (kẹp ≥ 0 — với θ̇ > Ω0 công thức cho ρ < 0, tức γ > 1, Bellman phân kỳ).
  Trụ 2  kênh đồng hồ [T_rem/H, sin(2πt/T0), cos(2πt/T0)] nối vào (cosθ, sinθ, θ̇/3).
  Trụ 3  chuyển bán-Markov (s, a, R_cum, Δt, s′): đồng hồ nhảy 1 + 9N nhịp khi gãy (đợt Geom), đích R + e^{−ρΔt} max Q̄.
  Trụ 4  mặt nạ hành động theo phản lực λ(s,u) (hàm chắn): loại u nếu đường đi tất định trong nhịp — xuất phát
         từ θ̇ + mσ(θ) (nhiễu-trước) — chạm |λ| ≥ R; m = 0 (λ tất định) hoặc m = 1.
  Mạng: tanh, đầu dueling Q = V + A − mean(A), Adam lr = 1e−3, Huber, cắt gradient 10.

Biến thể (bỏ từng trụ để đo đóng góp):
  E   đầy đủ: trụ 1 (ρ = 0.03) + trụ 2 + trụ 3 + trụ 4 (m = 1)
  E0  như E nhưng mặt nạ m = 0 (chỉ λ tất định, không tính nhiễu)
  Em  E bỏ trụ 4 (không mặt nạ)
  Ec  E bỏ sin/cos pha đồng hồ (chỉ giữ T_rem/H)
  Ev  E thay ρ hằng bằng ρ_eff(s) Van der Pol
Chạy (từ hoc-thuyet-e/): PYTHONPATH=.:rl_data python3 rl_data/giai_so_e_dqn.py train E 0 [iters] | eval
"""
import json
import math
import os
import sys
import time

import numpy as np

import giai_so_rl_bien_t as G
import giai_so_te_dqn as TE

K, TICK, NA, ACTIONS = TE.K, TE.TICK, TE.NA, TE.ACTIONS
RHO = 0.03
OMEGA0 = 1.05 * 3.2048
VARIANTS = {
    'E': dict(mask=1.0, clock=True, vdp=False, name='E-DQN đầy đủ'),
    'E0': dict(mask=0.0, clock=True, vdp=False, name='E-DQN, mặt nạ λ tất định'),
    'Em': dict(mask=None, clock=True, vdp=False, name='E-DQN bỏ mặt nạ λ'),
    'Ec': dict(mask=1.0, clock=False, vdp=False, name='E-DQN bỏ sin/cos đồng hồ'),
    'Ev': dict(mask=1.0, clock=True, vdp=True, name='E-DQN chiết khấu Van der Pol ρ_eff(s)'),
}


def feats(env, clock):
    t = env.k * TICK
    f = [np.cos(env.th), np.sin(env.th), env.w / 3.0, np.maximum(K - env.k, 0) / K]
    if clock:
        f += [np.sin(2 * math.pi * t / TE.T0), np.cos(2 * math.pi * t / TE.T0)]
    return np.stack(f, -1)


def action_mask(th, w, m):
    """True = được phép. Dự báo: dòng tất định một nhịp từ (θ, θ̇ + mσ) có chạm |λ| ≥ R không (λ gồm số hạng mκcosθ u/M)."""
    w0 = w + m * G.noise_sigma(th)
    ok = np.stack([~G.tick_step(th, w0, u)[4] for u in ACTIONS], -1)
    ok[~ok.any(1), 0] = True                   # mọi hành động đều gãy ⟹ cho phép phanh u = −0.2
    return ok


def rho_eff(w):
    return np.maximum(G.GTH / G.I * (1 - w * w / OMEGA0**2), 0.0)


class DuelingTanh:
    def __init__(self, din, dout, rng, hid=128):
        xa = lambda a, b: rng.standard_normal((a, b)) * math.sqrt(1.0 / a)
        self.p = [xa(din, hid), np.zeros(hid), xa(hid, hid), np.zeros(hid),
                  xa(hid, 1) * 0.1, np.zeros(1), xa(hid, dout) * 0.1, np.zeros(dout)]
        self.m = [np.zeros_like(x) for x in self.p]; self.v = [np.zeros_like(x) for x in self.p]; self.t = 0

    def forward(self, x, cache=False):
        W1, b1, W2, b2, Wv, bv, Wa, ba = self.p
        h1 = np.tanh(x @ W1 + b1); h2 = np.tanh(h1 @ W2 + b2)
        V = h2 @ Wv + bv; A = h2 @ Wa + ba
        if cache:
            self._c = (x, h1, h2)
        return V + A - A.mean(1, keepdims=True)

    def backward_adam(self, gq, lr):
        x, h1, h2 = self._c
        W1, b1, W2, b2, Wv, bv, Wa, ba = self.p
        gV = gq.sum(1, keepdims=True); gA = gq - gq.mean(1, keepdims=True)
        gWv, gbv, gWa, gba = h2.T @ gV, gV.sum(0), h2.T @ gA, gA.sum(0)
        g2 = (gV @ Wv.T + gA @ Wa.T) * (1 - h2 * h2)
        gW2, gb2 = h1.T @ g2, g2.sum(0)
        g1 = (g2 @ W2.T) * (1 - h1 * h1)
        gW1, gb1 = x.T @ g1, g1.sum(0)
        grads = [gW1, gb1, gW2, gb2, gWv, gbv, gWa, gba]
        gn = math.sqrt(sum(float((g * g).sum()) for g in grads))
        if gn > 10.0:
            grads = [g * (10.0 / gn) for g in grads]
        self.t += 1
        for i, g in enumerate(grads):
            self.m[i] = 0.9 * self.m[i] + 0.1 * g
            self.v[i] = 0.999 * self.v[i] + 0.001 * g * g
            self.p[i] -= lr * (self.m[i] / (1 - 0.9**self.t)) / (np.sqrt(self.v[i] / (1 - 0.999**self.t)) + 1e-8)

    def soft_update(self, src, tau):
        for i in range(len(self.p)):
            self.p[i] = (1 - tau) * self.p[i] + tau * src.p[i]


def greedy(net, x, ok):
    q = net.forward(x)
    return np.where(ok, q, -1e9).argmax(1)


def make_policy(net, cfg):
    def pol(e):
        ok = action_mask(e.th, e.w, cfg['mask']) if cfg['mask'] is not None else np.ones((e.n, NA), bool)
        return greedy(net, feats(e, cfg['clock']), ok)
    return pol


def train(var, seed, iters=30000, nenv=32, batch=256, updates=2, lr=1e-3, tau_soft=0.005, buf=200000):
    cfg = VARIANTS[var]
    rng = np.random.default_rng(seed)
    din = 6 if cfg['clock'] else 4
    net, tgt = DuelingTanh(din, NA, rng), DuelingTanh(din, NA, rng)
    tgt.p = [x.copy() for x in net.p]
    env = TE.EKepEnv(nenv, seed + 99)
    S = np.zeros((buf, din)); A = np.zeros(buf, np.int64); Rw = np.zeros(buf); S2 = np.zeros((buf, din))
    Dsc = np.zeros(buf); OK2 = np.ones((buf, NA), bool); ptr = size = 0
    allow = lambda e: action_mask(e.th, e.w, cfg['mask']) if cfg['mask'] is not None else np.ones((e.n, NA), bool)
    s, ok = feats(env, cfg['clock']), allow(env)
    best, best_p, curve, t0 = -1e9, None, [], time.time()
    pol = make_policy(net, cfg)
    for it in range(iters):
        eps = max(0.05, 1.0 - it / (0.5 * iters))
        lr_t = lr * (1.0 - 0.8 * it / iters)
        if it % 2500 == 0 or it == iters - 1:                          # đường học (kiểm định seed riêng)
            vm = TE.evaluate(pol, n=300, seed=777 + seed)[0]
            curve.append((it * nenv, vm))
            if it >= iters // 2 and vm > best:
                best, best_p = vm, [x.copy() for x in net.p]
        rnd = np.array([rng.choice(np.flatnonzero(o)) for o in ok])
        a = np.where(rng.random(nenv) < eps, rnd, greedy(net, s, ok))
        w_old = env.w.copy()
        r, dk, broke, end = env.step(a)
        rho = rho_eff(w_old) if cfg['vdp'] else RHO
        disc = np.where(end, 0.0, np.exp(-rho * TICK * dk))           # trụ 1 + 3: e^{−ρΔt}, Δt gồm cả đợt Carnot
        s2, ok2 = feats(env, cfg['clock']), allow(env)
        idx = (ptr + np.arange(nenv)) % buf
        S[idx], A[idx], Rw[idx], S2[idx], Dsc[idx], OK2[idx] = s, a, r, s2, disc, ok2
        ptr = (ptr + nenv) % buf; size = min(size + nenv, buf)
        if end.any():
            env.k[end] = 0; env.nbreak[end] = 0; env._relock(end)
            s2, ok2 = feats(env, cfg['clock']), allow(env)
        s, ok = s2, ok2
        if size >= 5000:
            for _ in range(updates):
                j = rng.integers(size, size=batch)
                a_star = np.where(OK2[j], net.forward(S2[j]), -1e9).argmax(1)       # Double DQN trong tập cho phép
                y = Rw[j] + Dsc[j] * tgt.forward(S2[j])[np.arange(batch), a_star]
                q = net.forward(S[j], cache=True)
                err = q[np.arange(batch), A[j]] - y
                g = np.zeros_like(q); g[np.arange(batch), A[j]] = np.clip(err, -1, 1) / batch
                net.backward_adam(g, lr_t)
                tgt.soft_update(net, tau_soft)
    net.p = best_p
    m_, se, nb = TE.evaluate(make_policy(net, cfg))
    # hành vi theo hạn chót: ⟨u⟩ của chính sách trên dải θ̇ ∈ (2.0, 2.9), θ đều, theo T_rem
    th = np.linspace(-math.pi, math.pi, 36, endpoint=False)
    TH, W = np.meshgrid(th, np.linspace(2.0, 2.9, 10)); TH, W = TH.ravel(), W.ravel()
    ubar = {}
    for left in (20, 5, 1, 0.4):
        class _E:
            pass
        e = _E(); e.th, e.w, e.n = TH, W, TH.size; e.k = np.full(TH.size, K - int(round(left * 25)))
        ubar[left] = float(ACTIONS[make_policy(net, cfg)(e)].mean())
    out = dict(var=var, name=cfg['name'], seed=seed, ret=m_, se=se, breaks=nb, val=best, curve=curve,
               ubar=ubar, sec=time.time() - t0, transitions=iters * nenv)
    json.dump(out, open(f'rl_data/e_dqn_{var}_{seed}.json', 'w'), ensure_ascii=False)
    print(json.dumps({k_: v for k_, v in out.items() if k_ != 'curve'}, ensure_ascii=False), flush=True)


def summarize():
    import glob
    rows = {}
    for f in sorted(glob.glob('rl_data/e_dqn_*.json')) + sorted(glob.glob('rl_data/te_dqn_[VSTH]_*.json')):
        r = json.load(open(f)); rows.setdefault((r['var'], r['name']), []).append(r)
    dp = json.load(open('rl_data/te_dqn_dp.json'))
    print(f"chuẩn DP phụ thuộc t: {dp['DP_t'][0]:.3f}; DP dừng: {dp['DP_dung'][0]:.3f}")
    for (v, name), rs in rows.items():
        rets = [r['ret'] for r in rs]
        print(f"{v:3s} {name:40s} n={len(rs)} lợi tức {np.mean(rets):.3f} (min {min(rets):.3f}, max {max(rets):.3f}) "
              f"= {100 * np.mean(rets) / dp['DP_t'][0]:.1f}% DP; gãy/tập {np.mean([r['breaks'] for r in rs]):.3f}"
              + (f"; ⟨u⟩@T_rem 20/5/1/0.4 T0: " + " ".join(f"{np.mean([r['ubar'][k_] for r in rs]):.3f}" for k_ in ('20', '5', '1', '0.4')) if 'ubar' in rs[0] else ""))
        if 'curve' in rs[0]:
            c = np.mean([[p[1] for p in r['curve']] for r in rs], 0)
            print("     đường học (chuyển → lợi tức kiểm định): " + " ".join(f"{p[0] // 1000}k:{x:.1f}" for p, x in zip(rs[0]['curve'], c)))


if __name__ == '__main__':
    if sys.argv[1] == 'train':
        train(sys.argv[2], int(sys.argv[3]), iters=int(sys.argv[4]) if len(sys.argv) > 4 else 30000)
    else:
        summarize()
