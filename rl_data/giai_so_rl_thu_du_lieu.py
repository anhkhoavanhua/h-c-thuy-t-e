"""Thu dữ liệu chu trình E_past → E_present → E_future cho RL bán-Markov (SMDP) theo tỉ số lái f.

Mô phỏng Langevin đầy đủ của Part I (giai_so_dong_co_langevin.run): γθ = 1/τ = 0.1, γx = 0.70, k_BT = 0.1 (1 MeV),
Ω0 = f·θ̇_gãy. Mỗi quỹ đạo lưu chuỗi sự kiện ('gãy', t, θ̇) / ('khóa', t, θ̇, ΔE_Carnot), công lái W và độ dài T.
Chạy: python3 giai_so_rl_thu_du_lieu.py   (8 lõi, ≈ 30–60 phút; tiếp tục được nếu bị ngắt) → rl_data/chu_trinh_f.json
"""
import json
import os
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L

F_LIST = (0.95, 1.05, 1.13)
N_TRAJ, T_LEN = 4, 120 * L.T0


def _one(a):
    f, seed = a
    r = L.run(0.1, f * L.WC, 0.1, T_LEN, 2e-3, seed, gx=0.70, w_start=0.82 * L.WC)
    return dict(f=f, seed=seed, T=T_LEN, W=r['W'], carnot=r['carnot'], events=r['events'])


if __name__ == '__main__':
    os.makedirs('rl_data', exist_ok=True)
    part = 'rl_data/chu_trinh_f.part.jsonl'                # lưu từng quỹ đạo để không mất khi bị ngắt
    done = set()
    if os.path.exists(part):
        done = {(d['f'], d['seed']) for d in map(json.loads, open(part))}
    jobs = [(f, 1000 * i + s) for i, f in enumerate(F_LIST) for s in range(N_TRAJ) if (f, 1000 * i + s) not in done]
    with Pool(12) as p, open(part, 'a') as fh:
        for r in p.imap_unordered(_one, jobs):
            fh.write(json.dumps(r) + '\n'); fh.flush()
            print('xong f =', r['f'], 'seed =', r['seed'], file=sys.stderr, flush=True)
    out = [json.loads(l) for l in open(part)]
    with open('rl_data/chu_trinh_f.json', 'w') as fh:
        json.dump(out, fh)
    print('đã lưu', len(out), 'quỹ đạo', file=sys.stderr)
