""""E năng lượng": năng lượng để thoát khỏi "0" (chân không / trạng thái cơ bản, nơi thời gian đóng băng, h = 0).

  A  Lượng tử (Part III, S 36 chiều; Carnot b = 1, bơm a quét): trạng thái dừng ρ_ss →
     năng lượng vượt chân không ΔE = ⟨H⟩ − E_0, công suất bơm vào P = Tr(H·D_B[ρ]), dòng bản ghi J_C,
     độ mở tương lai h (nhánh), giá năng lượng của thời gian P/h (ħω0 mỗi nat) so với Landauer k_BT (k_BT = 1 MeV).
     Chân không: a = 0 ⟹ ρ = cơ bản, năng lượng điểm không vẫn có nhưng h = 0.
  B  Cổ điển (Part I): hàng rào kích hoạt từ tần suất gãy đo theo k_BT (xung_dot_mo_bang_ket_qua.txt, phần C):
     rate = A e^{−Δ/k_BT} ⟹ Δ; so với Δ_eff ≈ 0.4 MeV của bài báo và khe tất định ½M θ̇_br²(1 − M_det).
Chạy: python3 giai_so_e_nang_luong.py   (≈ 1–2 phút)
"""
import math
import os
import re

import numpy as np

T0 = 2 * math.pi
HERE = os.path.dirname(os.path.abspath(__file__))


def dissip(L, r):
    return L @ r @ L.conj().T - 0.5 * (L.conj().T @ L @ r + r @ L.conj().T @ L)


def part_a():
    from giai_so_wdw_pw_e7 import channels, TAU
    from giai_so_kim_chi_gio import kraus, psi0
    from giai_so_xung_dot_mo_bang import steady, rate
    from giai_so_e_khuc_tuong_lai_mo import openness
    ls, H, (Lc, Lb) = channels()
    E, V = np.linalg.eigh(H); E0 = E[0]
    per = T0 / TAU
    rng = np.random.default_rng(21)
    print("=== A. Lượng tử: năng lượng vượt chân không và giá năng lượng của thời gian (Carnot b = 1) ===")
    print(f"  E_0 (chân không của S, gồm điểm không) = {E0:.3f} ħω0; khe tới mức kích thích đầu = {E[1] - E0:.3f} ħω0;"
          f" ngưỡng gãy E_br = 4.706")
    print("  a (bơm)   ΔE = ⟨H⟩−E0   P_vào (ħω0/T0)   J_C (/T0)   h (nat/T0)   P/h (ħω0/nat)   P/J_C (ħω0/bản ghi)")
    rows = []
    for a in (0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
        if a == 0:
            r = np.outer(V[:, 0], V[:, 0].conj()); P = 0.0; J = 0.0; h = 0.0
        else:
            Lb_a = math.sqrt(a) * Lb
            r = steady(H, [Lc, Lb_a])
            P = float(np.real(np.trace(H @ dissip(Lb_a, r)))) * T0
            J = rate(Lc, r) * T0
            h = float(openness(kraus(H, [Lc, Lb_a], TAU), psi0(), 1200, 150, rng, late=400).mean() * per)
        dE = float(np.real(np.trace(H @ r))) - E0
        ph = P / h if h > 1e-9 else float('nan'); pj = P / J if J > 1e-12 else float('nan')
        rows.append((a, dE, P, J, h, ph, pj))
        print(f"  {a:6.2f}    {dE:8.3f}      {P:9.4f}       {J:8.4f}    {h:8.3f}     {ph:8.3f}        {pj:8.3f}", flush=True)
    kT = 0.1
    print(f"  Landauer: tối thiểu k_BT = {kT} ħω0 (1 MeV) mỗi nat ghi lại; tỉ số P/h ÷ k_BT ở E_dual (a = 1): "
          f"{[r for r in rows if r[0] == 1.0][0][5] / kT:.1f}×\n")
    return rows


def part_b():
    print("=== B. Cổ điển (Part I): hàng rào kích hoạt từ tần suất gãy theo nhiệt độ ===")
    path = os.path.join(HERE, 'xung_dot_mo_bang_ket_qua.txt')
    txt = open(path).read()
    pts = [(float(a) / 10, float(b)) for a, b in re.findall(r"^\s+(\d\.\d)\s+(\d\.\d{4})\s+\d\.\d{4}\s*$", txt, re.M)]
    pts = [(k, r) for k, r in pts if k > 0 and r > 0]          # k_BT (ħω0), lần gãy / T0
    x = np.array([1 / k for k, _ in pts]); y = np.log([r for _, r in pts])
    slope, icpt = np.polyfit(x, y, 1)
    D = -slope
    print("  dữ liệu (k_BT MeV → gãy/T0): " + ", ".join(f"{10 * k:.1f} → {r:.3f}" for k, r in pts))
    print(f"  khớp Arrhenius rate = A e^(−Δ/k_BT): Δ = {D:.4f} ħω0 = {10 * D:.2f} MeV, A = {math.exp(icpt):.2f} /T0"
          f"   (bài báo: Δ_eff ≈ 0.4 MeV)")
    M = 1 + (1 / 3)**2 * 0.5                                   # quán tính hiệu dụng trung bình I + mκ²⟨cos²⟩
    gap = 0.5 * M * 3.2048**2 * (1 - 0.863)
    print(f"  khe tất định (λ ∝ θ̇²): ½⟨M⟩θ̇_br²(1 − M_det) = {gap:.3f} ħω0 = {10 * gap:.1f} MeV — lớn hơn Δ khớp"
          f" {gap / D:.0f}×: nhiễu không cần vượt cả khe, chỉ cần đẩy pha θ tới góc |λ| lớn nhất\n")
    return D, gap


if __name__ == '__main__':
    part_a()
    part_b()
