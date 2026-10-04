"""Giải số luận điểm: "xung đột Quá khứ ↔ Hiện tại làm thời gian không còn đóng băng" (HocThuyetE.pdf … (8)).

Phát biểu kiểm được: trong trạng thái dừng toàn cục (WdW), thời gian quan hệ chỉ "chạy" khi có dòng chu trình bất khả
nghịch J > 0 (bản ghi Carnot mới sinh ra không ngừng). J > 0 cần CẢ lực đẩy lên ngưỡng gãy (Quá khứ → Hiện tại: bơm/
động cơ) LẪN tiêu tán Carnot (Hiện tại → Tương lai/Quá khứ). Nếu một bên thắng tuyệt đối: J = 0 ⟹ đóng băng.

  A  Lượng tử (Part III, S 36 chiều, kênh Carnot L_C và bơm L_B; giai_so_wdw_pw_e7.channels):
     nhân hệ số a cho bơm, b cho Carnot; trạng thái dừng ρ_ss (nhân Liouvillian); dòng bản ghi J_C = Tr L_C ρ L_C†,
     J_B; tỉ lệ chiếm Quá khứ (E < E_br = 4.706 ħω0) và Hiện tại (E ≥ E_br); đóng băng sau khi tắt bơm.
  B  Đối chứng cân bằng: không bơm, bể Carnot ở k_BT > 0 có cân bằng chi tiết (Davies có kênh lên).
     Bản ghi vẫn sinh ra (J_xuống = J_lên > 0) nhưng dòng ròng mỗi kênh = 0 và tỉ lệ đảo chiều chuỗi bản ghi = 1:
     có "tích tắc" nhưng không có chiều — đóng băng theo nghĩa mũi tên.
  C  Cổ điển (Part I, Langevin đầy đủ): tần suất gãy theo k_BT — ở k_BT = 0 điểm hút nằm dưới ngưỡng 13.7%
     (Quá khứ thắng) ⟹ không gãy, không bản ghi; và theo tỉ số lái f (rl_data/chu_trinh_f.json).
Chạy: python3 giai_so_xung_dot_mo_bang.py [A|B|C]   (A, B < 1 phút; C ≈ 2–4 phút, 8 lõi)
"""
import json
import math
import os
import sys

import numpy as np

T0 = 2 * math.pi
E_BR = 4.706


def liou(H, Ls):
    d = H.shape[0]; I = np.eye(d); S = sum(L.conj().T @ L for L in Ls)
    M = -1j * (np.kron(H, I) - np.kron(I, H.T)) - 0.5 * (np.kron(S, I) + np.kron(I, S.T))
    for L in Ls:
        M = M + np.kron(L, L.conj())
    return M


def steady(H, Ls):
    d = H.shape[0]
    M = liou(H, Ls)
    b = np.zeros(d * d, complex); M[0, :] = np.eye(d).reshape(-1); b[0] = 1
    r = np.linalg.solve(M, b).reshape(d, d)
    return (r + r.conj().T) / 2


def rate(L, r):
    return float(np.real(np.trace(L @ r @ L.conj().T)))


def part_a():
    from giai_so_wdw_pw_e7 import channels
    ls, H, (Lc, Lb) = channels()
    E, V = np.linalg.eigh(H)
    Ppast = V[:, E < E_BR] @ V[:, E < E_BR].conj().T
    print("=== A. Lượng tử Part III: dòng bản ghi theo cường độ bơm a (Quá khứ → Hiện tại) và Carnot b ===")
    print(f"S {H.shape[0]} chiều; mức dưới ngưỡng gãy E_br = {E_BR}: {int((E < E_BR).sum())} / {E.size}")
    av = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0); bv = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0)
    J = np.zeros((len(av), len(bv))); PP = np.zeros_like(J)
    for i, a in enumerate(av):
        for j, b in enumerate(bv):
            if a == 0 and b == 0:                                  # chỉ H: không có trạng thái dừng duy nhất, J = 0
                J[i, j], PP[i, j] = 0.0, float('nan'); continue
            Ls = [L for L in (math.sqrt(b) * Lc, math.sqrt(a) * Lb) if np.abs(L).max() > 0]
            r = steady(H, Ls)
            J[i, j] = rate(math.sqrt(b) * Lc, r) * T0              # bản ghi Carnot mỗi T0
            PP[i, j] = float(np.real(np.trace(Ppast @ r)))
    print("  J_C (bản ghi Carnot / T0) ở trạng thái dừng; hàng a (bơm), cột b (Carnot):")
    print("   a\\b " + "".join(f"{b:>9}" for b in bv))
    for i, a in enumerate(av):
        print(f"  {a:5}" + "".join(f"{J[i, j]:9.3f}" for j in range(len(bv))))
    print("  P(Quá khứ) = Tr(P_{E<E_br} ρ_ss):")
    for i, a in enumerate(av):
        print(f"  {a:5}" + "".join(f"{PP[i, j]:9.3f}" for j in range(len(bv))))
    # J so với mức "đồng hiện diện" P_QK·P_HT
    m = (np.array(av)[:, None] > 0) & (np.array(bv)[None, :] > 0) & ~np.isnan(PP)
    coex = (PP * (1 - PP))[m]; jj = J[m]
    print(f"  tương quan(J_C, P_QK·P_HT) trên các ô a, b > 0: {np.corrcoef(coex, jj)[0, 1]:+.3f}")
    # đóng băng sau khi tắt bơm (từ trạng thái dừng a = b = 1)
    r0 = steady(H, [Lc, Lb])
    d = H.shape[0]
    M = liou(H, [Lc])
    w, Vv = np.linalg.eig(M); c = np.linalg.solve(Vv, r0.reshape(-1))
    print("  tắt bơm tại t = 0 (giữ Carnot): J_C(t) =", ", ".join(
        f"{t / T0:g} T0: {rate(Lc, (Vv @ (np.exp(w * t) * c)).reshape(d, d)) * T0:.4f}" for t in (0, T0, 5 * T0, 20 * T0, 100 * T0)))
    rs = steady(H, [Lb])
    print(f"  chỉ bơm (không Carnot): J_C = 0 (không có kênh), J_B dừng = {rate(Lb, rs) * T0:.4f} / T0;"
          f"  P(Quá khứ) = {np.real(np.trace(Ppast @ rs)):.3f}")
    rc = steady(H, [Lc])
    print(f"  chỉ Carnot: ρ_ss → cơ bản, P(Quá khứ) = {np.real(np.trace(Ppast @ rc)):.3f}, J_C = {rate(Lc, rc) * T0:.2e} / T0")
    print(f"  xung đột hai trạng thái đích: ‖ρ_ss(chỉ bơm) − ρ_ss(chỉ Carnot)‖_tr = "
          f"{0.5 * np.abs(np.linalg.eigvalsh(rs - rc)).sum():.3f}\n")
    return J, PP


def part_b():
    from giai_so_wdw_pw_e7 import channels
    from giai_so_part2_mo import davies
    ls, H, (Lc, Lb) = channels()
    u = None
    import giai_so_u_hermite_drude as UH
    _, _, _, _, H2, u = UH.basis(lmax=4, nmax=3)
    print("=== B. Đối chứng cân bằng: không bơm, bể Carnot ở k_BT có cân bằng chi tiết ===")
    print("  k_BT (ħω0)   J_xuống/T0   J_lên/T0   dòng ròng/T0   tỉ lệ lên/xuống   | E_dual (bơm + Carnot T=0): J_C/T0, đảo chiều")
    rE = steady(H, [Lc, Lb])
    jE = rate(Lc, rE) * T0
    for kT in (0.1, 0.3, 1.0, 3.0):
        Ls, Ee, Ve = davies(H, u, 0.5, kT=kT)
        r = steady(H, Ls)
        jd, ju = rate(Ls[0], r) * T0, rate(Ls[1], r) * T0
        print(f"   {kT:5.1f}       {jd:8.4f}    {ju:8.4f}    {jd - ju:+.2e}      {ju / jd:.4f}          |"
              f"  {jE:.3f}, lên/xuống = 0 (không có kênh lên)")
    print("  ⟹ cân bằng: có bản ghi (tích tắc) nhưng dòng ròng = 0, chuỗi bản ghi đảo chiều thời gian vẫn hợp lệ\n")


def _cl(args):
    import giai_so_dong_co_langevin as L
    f, kT, seed, T = args
    r = L.run(0.1, f * L.WC, kT, T, 2e-3, seed, gx=0.70, w_start=0.82 * L.WC)
    ev = r['events']
    return f, kT, sum(e[0] == 'gãy' for e in ev), r['carnot'], T


def part_c(T=40 * T0):
    from multiprocessing import Pool
    print("=== C. Cổ điển Part I: số lần gãy (Quá khứ → Hiện tại) theo nhiệt độ, f = 1.05 ===")
    jobs = [(1.05, kT, s, T) for kT in (0.0, 0.02, 0.05, 0.1) for s in range(2)]
    with Pool(8) as p:
        res = p.map(_cl, jobs)
    agg = {}
    for f, kT, nb, ec, TT in res:
        a = agg.setdefault(kT, [0, 0.0, 0.0]); a[0] += nb; a[1] += ec; a[2] += TT
    print("  k_BT (MeV)   lần gãy / T0   Carnot tiêu tán / T0 (ħω0)")
    for kT, (nb, ec, TT) in sorted(agg.items()):
        print(f"   {kT * 10:5.1f}       {nb / (TT / T0):.4f}         {ec / (TT / T0):.4f}")
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rl_data', 'chu_trinh_f.json')
    if os.path.exists(path):
        print("  theo tỉ số lái f (k_BT = 1 MeV, dữ liệu rl_data, 4 × 120 T0 mỗi f):")
        data = json.load(open(path))
        by = {}
        for tr in data:
            nb = sum(e[0] == 'gãy' for e in tr['events'])
            a = by.setdefault(tr['f'], [0, 0.0, 0.0]); a[0] += nb; a[1] += tr['carnot']; a[2] += tr['T']
        for f, (nb, ec, TT) in sorted(by.items()):
            print(f"   f = {f:.2f}: gãy / T0 = {nb / (TT / T0):.3f}, Carnot / T0 = {ec / (TT / T0):.3f}")
    print()


if __name__ == '__main__':
    for p in (sys.argv[1:] or ['A', 'B', 'C']):
        {'A': part_a, 'B': part_b, 'C': part_c}[p.upper()]()
