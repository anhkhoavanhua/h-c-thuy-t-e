"""Xung đột Hiện tại ↔ Tương lai (HocThuyetE, Part I): mỗi va chạm Carnot = một lần Tương lai thử tái khóa.
Tương lai thắng ⟹ đợt kết thúc (khóa giữ ≥ T0); Hiện tại thắng ⟹ gãy lại, va chạm tiếp. N ~ Geom(p) là thống kê xung đột.

Đo trên mô phỏng Langevin đầy đủ (giai_so_dong_co_langevin.run, γx = 0.70):
  1. Hàm nguy cơ h_k = P(Tương lai thắng ở va chạm k | đợt đã tới va chạm k): hằng số ⟺ xung đột không nhớ ⟺ Geom.
  2. Xóa trí nhớ: tương quan giữa các chu trình liên tiếp (N_i, N_{i+1}), (lưu_i, lưu_{i+1}), (N_i, lưu_{i+1}).
  3. Đổi cán cân: tắt lái trong pha tự do (nghiêng về Tương lai); thay k_BT, f.
  (Cấm tái khóa ⟹ không có va chạm ⟹ 0 bản ghi Carnot: theo cấu trúc mô hình, không cần chạy.)
Dữ liệu: rl_data/chu_trinh_f.json (f = 0.95/1.05/1.13, k_BT = 1 MeV, 4 × 120 T0) + các lượt mới chạy ở đây.
Chạy: python3 giai_so_xung_dot_hien_tai_tuong_lai.py   (≈ 3–5 phút, 8 lõi)
"""
import json
import math
import os
from multiprocessing import Pool

import numpy as np

import giai_so_burst_g as BG
import giai_so_dong_co_langevin as L

T0 = L.T0
HERE = os.path.dirname(os.path.abspath(__file__))


def _run(a):
    f, kT, seed, T, drive_free = a
    r = L.run(0.1, f * L.WC, kT, T, 2e-3, seed, gx=0.70, w_start=0.82 * L.WC, drive_free=drive_free)
    return dict(f=f, kT=kT, drive_free=drive_free, T=T, events=r['events'], carnot=r['carnot'])


def cycles(tr, warm=5 * T0):
    ev = [tuple(e) for e in tr['events']]
    dwell, bursts = BG.analyse(ev, warm)
    return dwell, bursts


def stats(trs, label):
    Ns, durs, dw_seq, N_seq, pres_time, Ttot, nrec = [], [], [], [], 0.0, 0.0, 0
    for tr in trs:
        dwell, bursts = cycles(tr)
        n = [b['n'] for b in bursts]
        Ns += n; durs += [b['dur'] for b in bursts]
        dw_seq.append(dwell[:len(bursts)]); N_seq.append(n)
        pres_time += sum(b['dur'] for b in bursts); Ttot += tr['T'] - 5 * T0
        nrec += sum(1 for e in tr['events'] if e[0] == 'khóa')
    Ns = np.array(Ns)
    if Ns.size < 5:
        print(f"  {label}: chỉ {Ns.size} đợt — không đủ thống kê"); return None
    # hàm nguy cơ theo chỉ số va chạm
    haz = []
    for k in range(1, 9):
        reach = (Ns >= k).sum(); end = (Ns == k).sum()
        if reach >= 8:
            haz.append((k, end / reach, reach))
    pML = 1 / Ns.mean()
    # tương quan giữa các chu trình liên tiếp
    c_nn, c_dd, c_nd = [], [], []
    for d, n in zip(dw_seq, N_seq):
        d = np.array(d); n = np.array(n)
        m = min(len(d), len(n))
        if m >= 3:
            c_nn += list(zip(n[:m - 1], n[1:m])); c_dd += list(zip(d[:m - 1], d[1:m])); c_nd += list(zip(n[:m - 1], d[1:m]))
    corr = lambda pairs: np.corrcoef(np.array(pairs).T)[0, 1] if len(pairs) > 5 else float('nan')
    # kiểm Geom: χ² của tần suất N so với Geom(p_ML)
    kmax = 12
    obs = np.array([(Ns == k).sum() for k in range(1, kmax)] + [(Ns >= kmax).sum()])
    expp = np.array([pML * (1 - pML)**(k - 1) for k in range(1, kmax)] + [(1 - pML)**(kmax - 1)]) * Ns.size
    ok = expp >= 3
    chi2 = float((((obs - expp)**2) / expp)[ok].sum()); dof = int(ok.sum()) - 2
    print(f"  {label}: {Ns.size} đợt, ⟨N⟩ = {Ns.mean():.2f}, p = 1/⟨N⟩ = {pML:.3f}, χ²(Geom) = {chi2:.1f} / {dof} bậc tự do;"
          f" Hiện tại chiếm {100 * pres_time / Ttot:.1f}% thời gian; bản ghi Carnot {nrec / (sum(t['T'] for t in trs) / T0):.2f} / T0")
    print("     nguy cơ h_k = P(Tương lai thắng ở va chạm k): " + ", ".join(f"k={k}: {h:.2f} (n={n})" for k, h, n in haz))
    print(f"     xóa trí nhớ: corr(N_i, N_i+1) = {corr(c_nn):+.2f}, corr(lưu_i, lưu_i+1) = {corr(c_dd):+.2f},"
          f" corr(N_i, lưu_i+1) = {corr(c_nd):+.2f}  ({len(c_nn)} cặp)")
    return dict(N=Ns.mean(), p=pML, haz=haz, chi2=chi2, dof=dof, pres=pres_time / Ttot)


def main():
    path = os.path.join(HERE, 'rl_data', 'chu_trinh_f.json')
    data = json.load(open(path)) if os.path.exists(path) else []
    print("=== 1–2. Xung đột Hiện tại ↔ Tương lai trong dữ liệu sẵn có (k_BT = 1 MeV, lái bật cả khi tự do) ===")
    for f in sorted({d['f'] for d in data}):
        stats([d for d in data if d['f'] == f], f"f = {f:.2f}")
    print("\n=== 3. Đổi cán cân xung đột (f = 1.05, 4 × 120 T0 mỗi trường hợp) ===")
    jobs = [(1.05, 0.1, 500 + s, 120 * T0, False) for s in range(4)] + \
           [(1.05, 0.05, 600 + s, 120 * T0, True) for s in range(4)]
    with Pool(8) as p:
        res = p.map(_run, jobs)
    stats([r for r in res if not r['drive_free']], "tắt lái khi tự do (nghiêng về Tương lai)")
    stats([r for r in res if r['drive_free'] and r['kT'] == 0.05], "k_BT = 0.5 MeV")

if __name__ == '__main__':
    main()
