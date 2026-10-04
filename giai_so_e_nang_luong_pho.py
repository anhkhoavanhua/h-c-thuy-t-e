"""E năng lượng là nhiều năng lượng khác nhau cấu thành + khối vũ trụ |Ψ_N⟩⟩ lớn dần theo số bản ghi.

  A  Thành phần năng lượng của động cơ thoát khỏi 0 (Part III, S 36 chiều, trạng thái dừng E_dual):
     ⟨H⟩ − ⟨H⟩_chân không tách thành động năng rô-to ℓ²/2I, thế −ΔE cosθ, dao động tử (n+½)ħω0, ghép −g x sinθ.
  B  Phổ lượng tử năng lượng: các vạch ω = E_j − E_i mà Carnot nhả ra / bơm nhận vào, tốc độ
     R_ij = |⟨i|L|j⟩|² p_j (gần đúng thế tục, p_j = ⟨j|ρ_ss|j⟩); phân loại theo Δ⟨|ℓ|⟩, Δ⟨n⟩ của trạng thái riêng.
     Kiểm: tổng công suất các vạch so với công suất chính xác Tr(H·D[ρ]).
  C  Khối vũ trụ lớn dần: |Ψ_N⟩⟩ (lịch sử N nhịp) chứa e^{H_N} nhánh hiệu dụng, H_N = ⟨s_N⟩ = entropy chuỗi bản ghi.
     Đo H_N theo t cho: PW đóng / chỉ Carnot / chỉ bơm / cân bằng nhiệt / E_dual.
Chạy: python3 giai_so_e_nang_luong_pho.py   (≈ 1 phút)
"""
import math
from collections import defaultdict

import numpy as np

T0 = 2 * math.pi


def setup():
    from giai_so_wdw_pw_e7 import channels, LMAX, NMAX
    import giai_so_u_hermite_drude as UH
    ls, H, (Lc, Lb) = channels()
    ls_, cos_t, sin_t, pth, H0, u = UH.basis(lmax=LMAX, nmax=NMAX)
    nr, no = len(ls), NMAX + 1
    KE = np.kron(np.diag(ls.astype(float)**2 / 2.0), np.eye(no))
    POT = -np.kron(cos_t, np.eye(no))
    OSC = np.kron(np.eye(nr), np.diag(np.arange(no) + 0.5))
    CPL = H - KE - POT - OSC
    Lop = np.kron(np.diag(np.abs(ls).astype(float)), np.eye(no))
    Nop = np.kron(np.eye(nr), np.diag(np.arange(no).astype(float)))
    return H, Lc, Lb, dict(KE=KE, POT=POT, OSC=OSC, CPL=CPL), Lop, Nop, u


def part_a(H, Lc, Lb, parts):
    from giai_so_xung_dot_mo_bang import steady
    E, V = np.linalg.eigh(H)
    vac = np.outer(V[:, 0], V[:, 0].conj())
    r = steady(H, [Lc, Lb])
    ev = lambda O, x: float(np.real(np.trace(O @ x)))
    print("=== A. Thành phần năng lượng: E_dual (trạng thái dừng) so với chân không (trạng thái cơ bản) ===")
    print("  thành phần                      chân không   E_dual    vượt chân không (ħω0)   (MeV)")
    names = dict(KE="động năng rô-to ℓ²/2I", POT="thế pha −ΔE cosθ", OSC="dao động tử (n+½)ħω0", CPL="ghép −g x sinθ")
    tot0 = tot1 = 0.0
    for k, nm in names.items():
        a, b = ev(parts[k], vac), ev(parts[k], r); tot0 += a; tot1 += b
        print(f"  {nm:30s}  {a:9.4f}  {b:9.4f}   {b - a:+9.4f}               {10 * (b - a):+7.2f}")
    print(f"  {'tổng ⟨H⟩':30s}  {tot0:9.4f}  {tot1:9.4f}   {tot1 - tot0:+9.4f}               {10 * (tot1 - tot0):+7.2f}")
    print(f"  (điểm không trong chân không: dao động tử ½ħω0 = 5 MeV có mặt nhưng h = 0)\n")
    return r


def part_b(H, Lc, Lb, r, Lop, Nop):
    E, V = np.linalg.eigh(H)
    p = np.real(np.diag(V.conj().T @ r @ V))
    lab_l = np.real(np.einsum('ij,ik,kj->j', V.conj(), Lop, V))
    lab_n = np.real(np.einsum('ij,ik,kj->j', V.conj(), Nop, V))
    print("=== B. Phổ các lượng tử năng lượng (vạch ω = E_j − E_i) ===")
    for nm, L, sign in (("Carnot nhả ra (Hiện tại → Tương lai)", Lc, +1), ("bơm nhận vào (Quá khứ → Hiện tại)", Lb, -1)):
        Le = V.conj().T @ L @ V
        R = np.abs(Le)**2 * p[None, :]                              # R[i, j]: j → i
        dE = E[None, :] - E[:, None]                                 # E_j − E_i (năng lượng nhả ra)
        lines = defaultdict(float); kinds = defaultdict(float); tot_r = tot_p = 0.0
        for i in range(len(E)):
            for j in range(len(E)):
                if R[i, j] < 1e-9 or i == j:
                    continue
                w = sign * dE[i, j]
                lines[round(w, 2)] += R[i, j] * T0
                dl = round(lab_l[j] - lab_l[i]); dn = round(lab_n[j] - lab_n[i])
                key = ("ℓ đổi" if dl != 0 else "") + ("+" if dl and dn else "") + ("n đổi" if dn != 0 else "")
                kinds[key or "lai (ℓ, n trộn)"] += R[i, j] * w * T0
                tot_r += R[i, j] * T0; tot_p += R[i, j] * w * T0
        exact = sign * float(np.real(np.trace(H @ (L @ r @ L.conj().T - 0.5 * (L.conj().T @ L @ r + r @ L.conj().T @ L))))) * T0
        exact = -exact if sign > 0 else -exact
        top = sorted(lines.items(), key=lambda x: -x[1] * abs(x[0]))[:10]
        n_lines = sum(1 for v in lines.values() if v > 1e-4)
        print(f"  {nm}: {n_lines} vạch có tốc độ > 1e−4/T0; tổng {tot_r:.3f} lượng tử/T0, công suất {tot_p:.4f} ħω0/T0"
              f" (Lindblad chính xác {abs(exact):.4f})")
        print("    10 vạch mang nhiều năng lượng nhất (ω ħω0 → tốc độ /T0, công suất ħω0/T0):")
        print("     " + "; ".join(f"{w:.2f} → {v:.3f} ({w * v:.3f})" for w, v in top))
        ws = np.array([w for w, v in lines.items() if v > 1e-4])
        print(f"    dải năng lượng: {ws.min():.2f} … {ws.max():.2f} ħω0 ({10 * ws.min():.1f} … {10 * ws.max():.1f} MeV)")
        print("    công suất theo loại chuyển: " + ", ".join(f"{k}: {v:.3f} ({100 * v / tot_p:.0f}%)"
                                                         for k, v in sorted(kinds.items(), key=lambda x: -x[1])))
    print()


def part_c(H, Lc, Lb, u):
    from giai_so_wdw_pw_e7 import TAU
    from giai_so_kim_chi_gio import kraus, psi0
    from giai_so_part2_mo import davies
    rng = np.random.default_rng(31)
    per = T0 / TAU
    regs = [("PW đóng", [0 * Lc]), ("chỉ Carnot (Quá khứ thắng)", [Lc]), ("chỉ bơm (Hiện tại thắng)", [Lb]),
            ("cân bằng nhiệt k_BT = 1", davies(H, u, 0.5, kT=1.0)[0]), ("E_dual (ba pha xung đột)", [Lc, Lb])]
    marks = [int(x * per) for x in (1, 5, 10, 20, 40)]
    print("=== C. Khối vũ trụ lớn dần: số nhánh hiệu dụng của |Ψ_N⟩⟩ = e^{H_N}, H_N = ⟨−ln p(r_1…r_N)⟩ ===")
    print("  chế độ                          " + "".join(f"  H @ {m / per:>4.0f} T0" for m in marks) + "   tăng 20→40 T0 (nat/T0)")
    for nm, Ls in regs:
        Ks = kraus(H, Ls, TAU)
        M = 200; psis = np.tile(psi0().astype(complex), (M, 1)); s = np.zeros(M); out = {}
        for k in range(1, marks[-1] + 1):
            outs = np.stack([psis @ K.T for K in Ks], 0)
            pp = np.real(np.einsum('jmd,jmd->jm', outs.conj(), outs)); pp = np.clip(pp, 1e-300, None); pp /= pp.sum(0, keepdims=True)
            u_ = rng.random(M); c = np.cumsum(pp, 0); j = (u_[None, :] > c).sum(0).clip(0, len(Ks) - 1)
            s += -np.log(pp[j, np.arange(M)])
            psis = outs[j, np.arange(M)] / np.sqrt(pp[j, np.arange(M)])[:, None]
            if k in marks:
                out[k] = s.mean()
        rate = (out[marks[-1]] - out[marks[-2]]) / ((marks[-1] - marks[-2]) / per)
        print(f"  {nm:32s}" + "".join(f"  {out[m]:10.2f}" for m in marks) + f"   {rate:8.3f}")
    print("  (số nhánh hiệu dụng = e^H; ví dụ H = 100 nat ⟹ e^100 ≈ 2.7·10^43 nhánh)\n")


if __name__ == '__main__':
    H, Lc, Lb, parts, Lop, Nop, u = setup()
    r = part_a(H, Lc, Lb, parts)
    part_b(H, Lc, Lb, r, Lop, Nop)
    part_c(H, Lc, Lb, u)
