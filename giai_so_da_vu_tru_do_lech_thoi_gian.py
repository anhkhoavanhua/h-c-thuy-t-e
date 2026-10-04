"""Độ lệch thời gian chuẩn của đa vũ trụ E_dual (mỗi vũ trụ một động cơ; WdW: trạng thái chung dừng, mọi thước đo
là tỉ số nhịp bản ghi giữa các vũ trụ — không đồng hồ ngoài).

  A  Cổ điển (Part I, Langevin đầy đủ): nhịp thời gian r(f) = số lần gãy / T0 theo tỉ số lái f trong CỬA SỔ VẬN HÀNH
     của bài báo f ∈ (0.90, 1.139), k_BT = 1 MeV. Đa vũ trụ = f phân bố đều trong cửa sổ ⟹
     σ_t = độ lệch chuẩn của ln(r/⟨r⟩), tỉ số giãn lớn nhất r_max/r_min, hệ số biến thiên CV.
  B  Lượng tử (Part III, S 36 chiều): r(a) chính xác (trạng thái dừng); hằng số nội tại (bão hòa, nửa bão hòa).
  C  Đa vũ trụ ghép trường trung bình (lượng tử): N động cơ, a_i log-đều trong [0.3, 3] (±nửa thập phân quanh E_dual);
     mỗi động cơ thêm tiêu tán sang biển chung √κ·L_C và nhận bơm từ biển chung tỉ lệ κ·J̄ (J̄ = nhịp trung bình,
     tự hợp) ⟹ σ_t(κ).
Chạy: python3 giai_so_da_vu_tru_do_lech_thoi_gian.py [A|B|C]   (A ≈ 5–10 phút 8 lõi; B, C < 1 phút)
"""
import json
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))
F_WIN = (0.90, 1.139)


def sigma_stats(fs, rs, fgrid):
    """Nội suy ln r(f) tuyến tính từng đoạn trên lưới đều của cửa sổ ⟹ thống kê đa vũ trụ."""
    lr = np.interp(fgrid, fs, np.log(rs))
    r = np.exp(lr)
    ref = r.mean()
    s = np.log(r / ref)
    return dict(sigma=float(s.std()), cv=float(r.std() / r.mean()), ratio=float(r.max() / r.min()),
                mean=float(ref), lo=float(r.min()), hi=float(r.max()))


def _cl(a):
    import giai_so_dong_co_langevin as L
    f, seed, T = a
    r = L.run(0.1, f * L.WC, 0.1, T, 2e-3, seed, gx=0.70, w_start=0.82 * L.WC)
    return f, sum(e[0] == 'gãy' for e in r['events']), sum(e[0] == 'khóa' for e in r['events']), T


def part_a():
    print("=== A. Cổ điển: nhịp thời gian trong cửa sổ vận hành f ∈ (0.90, 1.139) của bài báo ===")
    data = {}
    path = os.path.join(HERE, 'rl_data', 'chu_trinh_f.json')
    if os.path.exists(path):
        for tr in json.load(open(path)):
            d = data.setdefault(tr['f'], [0, 0, 0.0])
            d[0] += sum(e[0] == 'gãy' for e in tr['events']); d[1] += sum(e[0] == 'khóa' for e in tr['events']); d[2] += tr['T']
    new = [f for f in (0.90, 1.00, 1.10, 1.139) if f not in data]
    if new:
        with Pool(8) as p:
            res = p.map(_cl, [(f, 7000 + 10 * i + s, 100 * T0) for i, f in enumerate(new) for s in range(4)])
        for f, nb, nk, T in res:
            d = data.setdefault(f, [0, 0, 0.0]); d[0] += nb; d[1] += nk; d[2] += T
    fs = np.array(sorted(data)); rb = np.array([data[f][0] / (data[f][2] / T0) for f in fs])
    rk = np.array([data[f][1] / (data[f][2] / T0) for f in fs])
    print("   f       gãy / T0   bản ghi Carnot / T0   (số lần gãy)")
    for f, a, b in zip(fs, rb, rk):
        print(f"  {f:.3f}    {a:7.3f}     {b:7.3f}            ({data[f][0]})")
    fgrid = np.linspace(*F_WIN, 2001)
    for nm, r in (("theo lần gãy", rb), ("theo bản ghi Carnot", rk)):
        m = r > 0
        st = sigma_stats(fs[m], r[m], fgrid)
        print(f"  đa vũ trụ (f đều trong cửa sổ), đồng hồ {nm}: σ_t = std ln(r/⟨r⟩) = {st['sigma']:.3f}"
              f"  (≈ ×{math.exp(st['sigma']):.2f}), CV = {st['cv']:.3f}, r_max/r_min = {st['ratio']:.1f}")
    # bootstrap sai số (Poisson trên số đếm)
    rng = np.random.default_rng(3); sig = []
    for _ in range(2000):
        cnt = rng.poisson([data[f][1] for f in fs]); rr = cnt / np.array([data[f][2] / T0 for f in fs])
        if (rr > 0).all():
            sig.append(sigma_stats(fs, rr, fgrid)['sigma'])
    print(f"  bootstrap Poisson (bản ghi Carnot): σ_t = {np.mean(sig):.3f} ± {np.std(sig):.3f}\n")


def steady_rate(H, Lc, Lb):
    from giai_so_xung_dot_mo_bang import steady, rate
    r = steady(H, [Lc, Lb])
    return rate(Lc, r) * T0, float(np.real(np.trace(H @ r)))


def part_b():
    from giai_so_wdw_pw_e7 import channels
    ls, H, (Lc, Lb) = channels()
    E0 = np.linalg.eigvalsh(H)[0]
    print("=== B. Lượng tử: nhịp thời gian r(a) chính xác của một động cơ (Carnot b = 1) ===")
    print("   a        ΔE (ħω0)   r = J_C (/T0)   r/r(1)")
    r1 = steady_rate(H, Lc, Lb)[0]
    rows = []
    for a in (0.001, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100, 1000):
        r, e = steady_rate(H, Lc, math.sqrt(a) * Lb)
        rows.append((a, e - E0, r)); print(f"  {a:7.3f}   {e - E0:7.3f}    {r:9.4f}     {r / r1:7.3f}")
    rmax = rows[-1][2]
    half = min(rows, key=lambda x: abs(x[2] - rmax / 2))
    print(f"  bão hòa r_max ≈ {rmax:.3f}/T0 (a = 1000); nửa bão hòa ở a ≈ {half[0]};"
          f" giãn thời gian tối đa của vũ trụ 'chạy hết ga' so với E_dual: ×{rmax / r1:.2f};"
          f" so với a → 0: không chặn (r → 0, đóng băng)\n")
    return rows


def part_c(N=12):
    from giai_so_wdw_pw_e7 import channels
    from giai_so_xung_dot_mo_bang import steady, rate
    ls, H, (Lc, Lb) = channels()
    a = np.exp(np.linspace(math.log(0.3), math.log(3.0), N))
    print(f"=== C. Đa vũ trụ ghép trường trung bình: {N} động cơ, a log-đều trong [0.3, 3] ===")
    print("   κ       σ_t = std ln(r/⟨r⟩)   r_max/r_min   ⟨r⟩ (/T0)   J̄ tự hợp   vòng lặp")
    for kap in (0.0, 0.3, 1.0, 3.0, 10.0):
        Jbar = 0.5
        for it in range(60):
            rs = []
            for ai in a:
                Ls = [math.sqrt(1 + kap) * Lc, math.sqrt(ai + kap * Jbar) * Lb]
                rho = steady(H, Ls)
                rs.append(rate(Lc, rho) * T0 * (1 + kap))         # mọi lần tiêu tán (riêng + sang biển chung)
            rs = np.array(rs); new = rs.mean()
            if abs(new - Jbar) < 1e-4:
                break
            Jbar = 0.5 * Jbar + 0.5 * new
        s = np.log(rs / rs.mean())
        print(f"  {kap:5.1f}    {s.std():8.4f}               {rs.max() / rs.min():6.2f}      {rs.mean():.4f}     {Jbar:.4f}     {it + 1}")
    print()


if __name__ == '__main__':
    for p in (sys.argv[1:] or ['B', 'C', 'A']):
        {'A': part_a, 'B': part_b, 'C': part_c}[p.upper()]()
