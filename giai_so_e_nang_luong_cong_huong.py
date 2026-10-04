"""E năng lượng bên ngoài có cấu trúc: bể Carnot với mật độ phổ có đỉnh J(ω) — cộng hưởng với bậc thang rô-to?

Part III (S 36 chiều). Kênh Carnot chuẩn là bể phẳng: L_C = √0.5 · V triu(V†uV) V†. Bể có cấu trúc:
  L_C(ω_p) = √0.5 · V [triu(V†uV) ⊙ √f(ω_ij)] V†,  f(ω) = ε + (1−ε)·γ²/((ω−ω_p)² + γ²),  γ = 0.15 ħω0, ε = 0.05
(chiều cao đỉnh = bể phẳng, sàn 5% để không có trạng thái tối hoàn toàn). Bơm giữ nguyên (a = 1).
Quét ω_p: dòng bản ghi J_C, công suất nhả ra, độ mở tương lai h, năng lượng vượt chân không, vạch nhả mạnh nhất.
Thêm: bể HAI đỉnh (2.46 + 3.54) và bể đỉnh lệch (5.0) để so.
Chạy: python3 giai_so_e_nang_luong_cong_huong.py   (≈ 1–2 phút)
"""
import math

import numpy as np

T0 = 2 * math.pi
GAM, EPS = 0.15, 0.05


def structured(H, u, peaks, gam=None):
    gam = GAM if gam is None else gam
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    w = E[None, :] - E[:, None]                      # ω_ij = E_j − E_i (nhả ra khi j → i, i < j)
    f = np.full_like(w, EPS)
    for wp in peaks:
        f = f + (1 - EPS) * gam**2 / ((w - wp)**2 + gam**2) / len(peaks) * (1 if peaks else 0)
    if not peaks:
        f = np.ones_like(w)
    return V @ (math.sqrt(0.5) * np.triu(ue, 1) * np.sqrt(f)) @ V.conj().T


def lines(H, L, r):
    E, V = np.linalg.eigh(H)
    p = np.real(np.diag(V.conj().T @ r @ V))
    Le = V.conj().T @ L @ V
    R = np.abs(Le)**2 * p[None, :]
    dE = E[None, :] - E[:, None]
    acc = {}
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            if R[i, j] > 1e-9:
                k = round(dE[i, j], 2); acc[k] = acc.get(k, 0.0) + R[i, j] * T0
    return acc


def main():
    from giai_so_wdw_pw_e7 import channels, TAU, LMAX, NMAX
    from giai_so_kim_chi_gio import kraus, psi0
    from giai_so_xung_dot_mo_bang import steady, rate
    from giai_so_e_khuc_tuong_lai_mo import openness
    import giai_so_u_hermite_drude as UH
    ls, H, (Lc0, Lb) = channels()
    u = UH.basis(lmax=LMAX, nmax=NMAX)[5]
    E, V = np.linalg.eigh(H); E0 = E[0]
    per = T0 / TAU
    rng = np.random.default_rng(41)
    print(f"Kiểm: bể phẳng dựng lại = L_C gốc? ‖L − L_C‖ = {np.abs(structured(H, u, []) - Lc0).max():.1e}")
    cases = [("bể phẳng (chuẩn)", [])] + [(f"đỉnh ω_p = {wp:.2f}", [wp]) for wp in
                                          (0.5, 1.0, 1.21, 1.75, 2.0, 2.46, 3.0, 3.54, 4.2, 5.0, 6.0)] + \
            [("hai đỉnh 2.46 + 3.54", [2.46, 3.54])]
    print("bể                          J_C (/T0)   P_nhả (ħω0/T0)   h (nat/T0)   ΔE (ħω0)   vạch mạnh nhất (ħω0: công suất)")
    rows = []
    for nm, pk in cases:
        Lc = structured(H, u, pk)
        r = steady(H, [Lc, Lb])
        J = rate(Lc, r) * T0
        D = Lc @ r @ Lc.conj().T - 0.5 * (Lc.conj().T @ Lc @ r + r @ Lc.conj().T @ Lc)
        P = -float(np.real(np.trace(H @ D))) * T0
        h = float(openness(kraus(H, [Lc, Lb], TAU), psi0(), 1000, 120, rng, late=300).mean() * per)
        dE = float(np.real(np.trace(H @ r))) - E0
        ln = lines(H, Lc, r)
        top = sorted(ln.items(), key=lambda kv: -kv[0] * kv[1])[:3]
        rows.append((nm, J, P, h, dE))
        print(f"{nm:28s} {J:8.3f}     {P:8.3f}          {h:7.3f}     {dE:7.3f}    "
              + ", ".join(f"{w:.2f}: {w * v:.3f}" for w, v in top), flush=True)
    base = rows[0]
    best = max(rows[1:], key=lambda x: x[1])
    worst = min(rows[1:], key=lambda x: x[1])
    print(f"\nmạnh nhất: {best[0]} — J_C ×{best[1] / base[1]:.2f}, h ×{best[3] / base[3]:.2f} so với bể phẳng;"
          f" yếu nhất: {worst[0]} — J_C ×{worst[1] / base[1]:.2f}, h ×{worst[3] / base[3]:.2f}")


def pump_structured(H, Lb, peaks, gam=GAM):
    """Nguồn bơm có cấu trúc: chuyển j → i hấp thụ ω = E_i − E_j; nhân biên độ với √f(ω) (ω ≤ 0 nhận sàn ε)."""
    E, V = np.linalg.eigh(H)
    Le = V.conj().T @ Lb @ V
    w = E[:, None] - E[None, :]
    f = np.full_like(w, EPS)
    for wp in peaks:
        f = f + (1 - EPS) * gam**2 / ((w - wp)**2 + gam**2) / len(peaks)
    return V @ (Le * np.sqrt(f)) @ V.conj().T


def measure(H, Lc, Lb, rng, E0):
    from giai_so_wdw_pw_e7 import TAU
    from giai_so_kim_chi_gio import kraus, psi0
    from giai_so_xung_dot_mo_bang import steady, rate
    from giai_so_e_khuc_tuong_lai_mo import openness
    r = steady(H, [Lc, Lb])
    J = rate(Lc, r) * T0
    D = Lb @ r @ Lb.conj().T - 0.5 * (Lb.conj().T @ Lb @ r + r @ Lb.conj().T @ Lb)
    Pin = float(np.real(np.trace(H @ D))) * T0
    h = float(openness(kraus(H, [Lc, Lb], TAU), psi0(), 1000, 120, rng, late=300).mean() * (T0 / TAU))
    return J, Pin, h, float(np.real(np.trace(H @ r))) - E0


def extra():
    from giai_so_wdw_pw_e7 import channels, LMAX, NMAX
    import giai_so_u_hermite_drude as UH
    ls, H, (Lc0, Lb0) = channels()
    u = UH.basis(lmax=LMAX, nmax=NMAX)[5]
    E0 = np.linalg.eigvalsh(H)[0]
    rng = np.random.default_rng(43)
    print("\n=== Quét độ rộng đỉnh γ (bể Carnot, bơm phẳng) ===")
    print("  γ (ħω0)   một đỉnh 3.54: J_C, h        hai đỉnh 2.46+3.54: J_C, h      (bể phẳng: 0.494, 5.23)")
    for g in (0.05, 0.15, 0.3, 0.6, 1.0, 2.0, 4.0):
        a = measure(H, structured(H, u, [3.54], g), Lb0, rng, E0)
        b = measure(H, structured(H, u, [2.46, 3.54], g), Lb0, rng, E0)
        print(f"  {g:5.2f}     {a[0]:.3f}, {a[2]:.2f}                  {b[0]:.3f}, {b[2]:.2f}", flush=True)
    print("\n=== Nguồn bơm cũng có cấu trúc (γ = 0.15) ===")
    print("  bơm                        bể Carnot                 J_C (/T0)  P_vào (ħω0/T0)  h (nat/T0)  ΔE (ħω0)")
    pk = [2.46, 3.54]
    combos = [("phẳng", None, "phẳng", None), ("đỉnh 2.46+3.54", pk, "phẳng", None),
              ("phẳng", None, "đỉnh 2.46+3.54", pk), ("đỉnh 2.46+3.54", pk, "đỉnh 2.46+3.54", pk),
              ("đỉnh 1.0 (lệch)", [1.0], "đỉnh 2.46+3.54", pk), ("đỉnh 2.46+3.54", pk, "đỉnh 1.0 (lệch)", [1.0])]
    for nb, pb, nc, pc in combos:
        Lb = Lb0 if pb is None else pump_structured(H, Lb0, pb)
        Lc = Lc0 if pc is None else structured(H, u, pc)
        J, P, h, dE = measure(H, Lc, Lb, rng, E0)
        print(f"  {nb:24s}   {nc:24s}  {J:7.3f}    {P:8.3f}       {h:6.2f}     {dE:6.2f}", flush=True)


if __name__ == '__main__':
    import sys
    if 'extra' in sys.argv:
        extra()
    else:
        main()
