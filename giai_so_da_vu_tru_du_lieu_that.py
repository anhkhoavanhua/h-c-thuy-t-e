"""Độ lệch thời gian chuẩn σ_t với DỮ LIỆU THẬT: mỗi hạt nhân (hoặc mỗi thí nghiệm) là một động cơ E_dual.
(Không có dữ liệu về các vũ trụ khác; đây là "đa hạt nhân" — miền của chính Học thuyết E.)

Ánh xạ:
  thang riêng ħω0(A) = 41 A^(−1/3) MeV (hệ thống học vỏ, Bohr–Mottelson) ⟹ T0(A) = 2πħ/ħω0(A)
  E năng lượng ε (đơn vị ħω0(A)) ⟹ cường độ lái a với ΔE(a) = ε (đường cong chính xác của động cơ Part III)
  nhịp không thứ nguyên r(a) (bản ghi Carnot / T0) ⟹ nhịp vật lý r_phys = r(a)/T0(A)
  σ_t = std ln(r_phys/⟨r_phys⟩)  (chỉ tỉ số giữa các đồng hồ — trạng thái chung vẫn tĩnh, WdW)

  I   IAEA (nds.iaea.org/relnsd, ground_states, tải 2026-10-02): hạt nhân bền (có độ phổ cập tự nhiên),
      ε = S_n/ħω0(A) — năng lượng hạt nhân giữ được trước khi neutron thoát ("thoát ly khỏi 0").
  II  EXFOR (rl_data/exfor_data, 7 bộ đo phổ neutron Cu/Ni/Ge/⁶⁰Ni(p,xn) 14–22 MeV):
      ε = E* / ħω0(A), E* = E_tới + S_tách của hạt nhân hợp phần (như giai_so_a9_exfor).
Chạy: python3 giai_so_da_vu_tru_du_lieu_that.py   (≈ 30 giây)
"""
import csv
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HBAR = 6.582119569e-22          # MeV·s


def curve():
    """ΔE(a), r(a) chính xác (Part III, Carnot b = 1) trên lưới log a."""
    from giai_so_wdw_pw_e7 import channels
    from giai_so_xung_dot_mo_bang import steady, rate
    ls, H, (Lc, Lb) = channels()
    E0 = np.linalg.eigvalsh(H)[0]
    As = np.logspace(-3, 3, 49)
    dE, r = [], []
    for a in As:
        rho = steady(H, [Lc, math.sqrt(a) * Lb])
        dE.append(float(np.real(np.trace(H @ rho))) - E0); r.append(rate(Lc, rho) * 2 * math.pi)
    return As, np.array(dE), np.array(r)


def rate_from_eps(eps, As, dE, r):
    """ε → a qua ΔE(a) (đơn điệu) → r. Trả r và cờ ε nằm ngoài dải mô hình."""
    out = (eps < dE[0]) | (eps > dE[-1])
    la = np.interp(np.clip(eps, dE[0], dE[-1]), dE, np.log(As))
    return np.interp(la, np.log(As), r), out


def stats(rphys):
    s = np.log(rphys / rphys.mean())
    return s.std(), rphys.max() / rphys.min()


def iaea(As, dE, r):
    path = os.path.join(HERE, 'rl_data', 'iaea_ground_states.csv')
    rows = []
    for d in csv.DictReader(open(path)):
        try:
            Z, N = int(d['z']), int(d['n'])
            ab = d['abundance'].strip(); sn = d['sn'].strip()
            if not ab or not sn or float(ab) <= 0 or Z < 1:
                continue
            rows.append((Z, N, Z + N, float(sn) / 1000.0, d['symbol']))
        except ValueError:
            continue
    A = np.array([x[2] for x in rows], float); Sn = np.array([x[3] for x in rows])
    hw = 41.0 * A**(-1 / 3)
    eps = Sn / hw
    rr, out = rate_from_eps(eps, As, dE, r)
    T0 = 2 * math.pi * HBAR / hw
    rphys = rr / T0
    print(f"=== I. IAEA: {len(rows)} hạt nhân có độ phổ cập tự nhiên (bền/sống rất lâu), A = {int(A.min())}–{int(A.max())} ===")
    print(f"  ħω0(A) = {hw.min():.1f}–{hw.max():.1f} MeV;  S_n = {Sn.min():.2f}–{Sn.max():.2f} MeV;"
          f"  ε = S_n/ħω0 = {eps.min():.2f}–{eps.max():.2f} (dải mô hình ΔE {dE[0]:.2f}–{dE[-1]:.2f}; ngoài dải: {out.sum()})")
    for nm, m in (("mọi hạt nhân bền", np.ones_like(A, bool)), ("A = 40–100 (miền A ≈ 70 của học thuyết)", (A >= 40) & (A <= 100)),
                  ("A ≥ 20", A >= 20)):
        s_all, ratio = stats(rphys[m])
        s_T0, _ = stats(1 / T0[m])
        s_eng, _ = stats(rr[m])
        print(f"  {nm:42s} n = {m.sum():3d}: σ_t = {s_all:.3f} (×{math.exp(s_all):.2f}), r_max/r_min = {ratio:.1f}"
              f"  | chỉ thang T0(A): {s_T0:.3f}; chỉ trạng thái động cơ r(ε): {s_eng:.3f}")
    i_fast, i_slow = int(np.argmax(rphys)), int(np.argmin(rphys))
    f = lambda i: f"{rows[i][4]}-{int(A[i])} (S_n = {Sn[i]:.2f} MeV, ε = {eps[i]:.2f})"
    print(f"  đồng hồ nhanh nhất: {f(i_fast)};  chậm nhất: {f(i_slow)}")
    ib = int(np.argmin(np.abs(A - 70)))
    print(f"  tham chiếu A ≈ 70 ({rows[ib][4]}-{int(A[ib])}): nhịp vật lý {rphys[ib]:.3e} bản ghi/s\n")
    return rphys, A


def exfor(As, dE, r):
    import giai_so_a9_exfor as X
    X.DIR = os.path.join(HERE, 'rl_data', 'exfor_data')            # thư mục EXFOR đã được dời vào rl_data/
    sets = X.datasets()
    # S_tách của hạt nhân hợp phần (MeV): n + ⁶³,⁶⁵Cu → S_n(⁶⁴,⁶⁶Cu) ≈ 7.9/7.1 — dùng giá trị cột cuối của bộ dữ liệu
    print(f"=== II. EXFOR: {len(sets)} bộ đo phổ neutron thật (miền A ≈ 60–74) ===")
    print("  bộ                                        A      E_tới   S_tách   E* (MeV)   ħω0(A)   ε = E*/ħω0   r(ε) /T0    nhịp vật lý (/s)")
    rp = []
    for nm, A, En, _, es, s, Ssep in sets:
        Estar = En + Ssep
        hw = 41.0 * A**(-1 / 3)
        eps = Estar / hw
        rr, out = rate_from_eps(np.array([eps]), As, dE, r)
        T0 = 2 * math.pi * HBAR / hw
        rp.append(rr[0] / T0)
        print(f"  {nm:40s} {A:5.1f}  {En:6.1f}   {Ssep:5.1f}    {Estar:6.1f}    {hw:6.2f}    {eps:7.2f}{' *' if out[0] else '  '}"
              f"    {rr[0]:7.3f}     {rr[0] / T0:.3e}")
    rp = np.array(rp)
    s, ratio = stats(rp)
    print(f"  σ_t (7 bộ EXFOR) = {s:.3f} (×{math.exp(s):.2f}), r_max/r_min = {ratio:.2f}\n")
    return rp


def main():
    As, dE, r = curve()
    print(f"Đường cong động cơ (Part III): ΔE {dE[0]:.2f} → {dE[-1]:.2f} ħω0 khi a = 10⁻³ → 10³; r {r[0]:.4f} → {r[-1]:.3f} /T0\n")
    rI, AI = iaea(As, dE, r)
    rX = exfor(As, dE, r)


if __name__ == '__main__':
    main()
