"""Xung đột Hiện tại ↔ Tương lai trong hạt nhân "tán xạ" tương lai — và nguồn năng lượng cần thiết.

Mỗi va chạm Carnot = một lần xung đột Hiện tại ↔ Tương lai (xung-dot-hien-tai-tuong-lai-kiem-chung.md).
Năng lượng Carnot tích lũy E_exc(k) = Σ_{i≤k} ΔE_i; khi E_exc ≥ S_n, neutron có thể thoát (Part I, V.A):
phổ ε ~ (E_exc − S_n)·Beta(2, n−1), n = 2k ± 1 (Ericson).

  A  Theo nguồn năng lượng (tỉ số lái f, dữ liệu Langevin đầy đủ rl_data/chu_trinh_f.json, k_BT = 1 MeV):
     số lần xung đột mỗi đợt, P(đợt đạt S_n) với S_n THẬT (IAEA) của hạt nhân hợp phần ⁶⁴Cu, ⁶⁶Cu, ⁶⁰Ni, ⁶²Ni, ⁷⁴Ge
     và S_n = 8 MeV của bài báo; lần xung đột thứ mấy thì vượt S_n; tần suất "tán xạ" (đợt có phát) mỗi T0.
  B  Tương lai tách nhánh tại mỗi lần xung đột: kết cục {Tương lai thắng (tái khóa giữ), Hiện tại thắng và chưa đủ S_n,
     Hiện tại thắng và đã vượt S_n (tán xạ)} — entropy H_k (nat) theo chỉ số xung đột k.
  C  Phổ năng lượng tán xạ ra (neutron) và giá năng lượng: Carnot / T0 chia cho nat tương lai / T0.
Chạy: python3 giai_so_tan_xa_tuong_lai.py   (≈ 10 giây)
"""
import csv
import json
import math
import os

import numpy as np

import giai_so_burst_g as BG

T0 = 2 * math.pi
HW = 10.0                                        # ħω0 = 10 MeV
HERE = os.path.dirname(os.path.abspath(__file__))


def sn_real():
    out = {}
    for d in csv.DictReader(open(os.path.join(HERE, 'rl_data', 'iaea_ground_states.csv'))):
        key = f"{d['symbol']}{int(d['z']) + int(d['n'])}"
        if key in ('Cu64', 'Cu66', 'Ni60', 'Ni62', 'Ge74') and d['sn'].strip():
            out[key] = float(d['sn']) / 1000.0
    return out


def bursts_by_f():
    data = json.load(open(os.path.join(HERE, 'rl_data', 'chu_trinh_f.json')))
    res = {}
    for tr in data:
        ev = [tuple(e) for e in tr['events']]
        _, bursts = BG.analyse(ev, 5 * T0)
        o = res.setdefault(tr['f'], dict(b=[], T=0.0))
        o['b'] += [b['pulses'] for b in bursts]; o['T'] += tr['T'] - 5 * T0
    return res


def main():
    SN = sn_real(); SN['bài báo'] = 8.0
    res = bursts_by_f()
    rng = np.random.default_rng(17)
    print("S_n thật (IAEA, MeV): " + ", ".join(f"{k}: {v:.2f}" for k, v in SN.items() if k != 'bài báo'))
    print("\n=== A. Nguồn năng lượng (tỉ số lái f) → tán xạ tương lai ===")
    print("  f      đợt   ⟨N⟩ xung đột   ⟨E_exc⟩ đợt (MeV)   " + "   ".join(f"P(≥S_n {k})" for k in SN) + "   tán xạ/T0 (S_n = 8)")
    summ = {}
    for f in sorted(res):
        bs = res[f]['b']
        Eb = np.array([sum(p) * HW for p in bs]); Nb = np.array([len(p) for p in bs])
        P = {k: float((Eb >= v).mean()) for k, v in SN.items()}
        rate = P['bài báo'] * len(bs) / (res[f]['T'] / T0)
        summ[f] = (len(bs), Nb.mean(), Eb.mean(), P, rate)
        print(f"  {f:.2f}   {len(bs):3d}     {Nb.mean():5.2f}          {Eb.mean():6.2f}          "
              + "      ".join(f"{P[k]:.2f}" for k in SN) + f"          {rate:.3f}")
    print("\n  Lần xung đột thứ mấy thì E_exc vượt S_n = 8 MeV (trong các đợt có vượt):")
    for f in sorted(res):
        ks = []
        for p in res[f]['b']:
            c = np.cumsum(np.array(p) * HW)
            hit = np.nonzero(c >= 8.0)[0]
            if hit.size:
                ks.append(hit[0] + 1)
        ks = np.array(ks)
        print(f"    f = {f:.2f}: trung vị k = {np.median(ks):.0f}, TB {ks.mean():.2f}, phân bố k=1..6+: "
              + ", ".join(f"{(ks == k).mean():.2f}" for k in range(1, 6)) + f", {(ks >= 6).mean():.2f}")

    print("\n=== B. Tương lai tách nhánh tại lần xung đột thứ k (f = 1.05, S_n = 8 MeV) ===")
    print("   k   đợt tới k   P(Tương lai thắng)   P(Hiện tại thắng, < S_n)   P(Hiện tại thắng, ≥ S_n: tán xạ)   H_k (nat)")
    bs = res[1.05]['b']
    for k in range(1, 9):
        reach = [p for p in bs if len(p) >= k]
        if len(reach) < 8:
            break
        win = np.mean([len(p) == k for p in reach])
        cum = [sum(p[:k]) * HW for p in reach]
        sc = np.mean([(len(p) > k) and (c >= 8.0) for p, c in zip(reach, cum)])
        ns = 1 - win - sc
        q = np.array([win, ns, sc]); q = q[q > 0]
        print(f"  {k:2d}     {len(reach):3d}          {win:.2f}                  {ns:.2f}                       {sc:.2f}"
              f"                     {-(q * np.log(q)).sum():.2f}")

    print("\n=== C. Phổ năng lượng tán xạ ra (neutron, Ericson Beta(2, n−1), n = 2k±1) và giá năng lượng ===")
    for f in sorted(res):
        eps = []
        for p in res[f]['b']:
            E = sum(p) * HW
            if E >= 8.0:
                n = 2 * len(p) + rng.choice([-1, 1])
                n = max(n, 2)
                eps.append((E - 8.0) * rng.beta(2, n - 1))
        eps = np.array(eps)
        if eps.size:
            print(f"  f = {f:.2f}: {eps.size} lần tán xạ; ⟨ε⟩ = {eps.mean():.2f} MeV, trung vị {np.median(eps):.2f},"
                  f" 95% < {np.percentile(eps, 95):.2f} MeV, P(ε > 3 MeV) = {(eps > 3).mean():.2f}")
    # giá năng lượng cho mỗi nat tương lai từ xung đột Hiện tại ↔ Tương lai (lượng tử, e-khuc mục D + e-nang-luong mục A)
    nat_per_C, E_per_C = 3.79, 2.77 * HW
    print(f"  lượng tử (Part III): mỗi lần xung đột Hiện tại ↔ Tương lai mở {nat_per_C:.2f} nat tương lai,"
          f" tốn ≈ {E_per_C:.1f} MeV ⟹ {E_per_C / nat_per_C:.1f} MeV mỗi nat tương lai bị tán xạ")


if __name__ == '__main__':
    main()
