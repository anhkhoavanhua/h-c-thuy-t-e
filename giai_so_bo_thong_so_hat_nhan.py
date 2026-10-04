"""Giải số "bộ thông số hoàn thiện" đề xuất: ħω0 = 10 MeV, χ = 1 (E_kép = ΔE = 10 MeV), I = χħ/ω0,
E_vt = 30 MeV (L = ħc/E_vt), e ∈ [0, 0.2], cửa sổ T_relax ≈ 1e-22 s.

(A) SI: kiểm từng số trong bảng đề xuất (λ de Broglie, T0, τ_relax/T0, I, ℓ0, κ, R, v/c).
(B) Động lực học đầy đủ (đơn vị m=ω0=ħ=ℓ0=1 ⟹ năng lượng tính theo ħω0 = 10 MeV):
    góc giữ tĩnh, θ̇_gãy chính xác, năng lượng quay tới hạn; mô phỏng chu kỳ gãy/khóa để đếm số va chạm
    và thời gian tái đồng bộ, so với các cửa sổ τ_relax, 5T0, 10T0.

Chạy: python3 giai_so_bo_thong_so_hat_nhan.py [A|B]   (B ≈ vài phút)
"""
import math
import sys

from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params, static_hold_angle
from giai_so_va_cham_tai_lien_ket import simulate

HBAR = 1.054571817e-34      # J·s
HBAR_MEV_S = 6.582119569e-22  # MeV·s
HBAR_C = 197.3269804        # MeV·fm
MEV = 1.602176634e-13       # J
C = 2.99792458e8
MPC2 = 938.27208816         # MeV

HW0, CHI, EVT = 10.0, 1.0, 30.0


def part_a():
    print("=== A. Kiểm các số SI của bộ thông số (ħω0=10 MeV, χ=1, E_vt=30 MeV) ===")
    for Ek in (5.0,):
        pc = math.sqrt(2 * MPC2 * Ek)
        print(f"de Broglie nucleon E_k={Ek} MeV: pc={pc:.2f} MeV, ƛ=ħc/pc={HBAR_C/pc:.3f} fm, "
              f"λ=2πƛ={2*math.pi*HBAR_C/pc:.3f} fm   (đề xuất: 20–40 fm)")
    w0 = HW0 * MEV / HBAR
    T0 = 2 * math.pi * HBAR_MEV_S / HW0
    tau = 5e-15 / (0.2 * C)
    print(f"T0 = 2πħ/ħω0 = {T0:.3e} s;  τ_relax = 5 fm/(0.2c) = {tau:.3e} s = {tau/T0:.3f} T0   (đề xuất: 10–30 T0 / ~5 T0)")
    print(f"   (nếu ħω0 = 1 MeV: T0 = {2*math.pi*HBAR_MEV_S:.3e} s, τ/T0 = {tau/(2*math.pi*HBAR_MEV_S):.4f})")
    kF = 1.36  # fm^-1, mật độ hạt nhân bão hòa
    pF = HBAR_C * kF
    print(f"v_F thực (k_F=1.36 fm⁻¹): p_F c={pF:.1f} MeV, v_F/c≈{pF/math.sqrt(pF**2+MPC2**2):.3f}")
    I = CHI * HBAR / w0
    l0 = HBAR_C / math.sqrt(MPC2 * HW0)
    kappa = l0 * CHI * HW0 / EVT
    R_mev_fm = EVT / l0
    print(f"I = χħ/ω0 = {I:.3e} kg·m²   (đề xuất: 6.94e-58)")
    print(f"ℓ0 = ħc/√(mc²·ħω0) = {l0:.4f} fm;  κ = ℓ0 E_kép/E_vt = {kappa:.4f} fm   (đề xuất: ≈0.68 fm)")
    print(f"R = E_vt/ℓ0 = {R_mev_fm:.3f} MeV/fm = {R_mev_fm*MEV/1e-15:.1f} N   (đề xuất: ≈1.5e4 N)")
    print(f"L = ħc/E_vt = {HBAR_C/EVT:.3f} fm;  tỉ số E_vt/(√χ ħω0) = {EVT/(math.sqrt(CHI)*HW0):.3f}")
    print(f"kκ/R = {CHI*(HW0/EVT)**2:.4f};  E_quay tiệm cận = E_vt²/(2ħω0) = {EVT**2/(2*HW0):.2f} MeV "
          f"= {EVT**2/(2*HW0)/MPC2:.4f} m_p c²;  v/c = {math.sqrt(CHI*HW0/MPC2):.4f}")
    print(f"Tỉ số 2–5 với ħω0=10 MeV ⟹ E_vt = 20–50 MeV ⟹ L = {HBAR_C/50:.2f}–{HBAR_C/20:.2f} fm   (đề xuất tiêu đề: 20–80 fm)\n")


def part_b(g=0.0, windows=None):
    eps = EVT / HW0
    p, R = params(CHI, eps, I=CHI, g=g)
    th = static_hold_angle(p, R)
    wc = omega_crit_general(p, R)
    E_at = lambda w: 0.5 * (p['I'] + p['m'] * p['kappa']**2) * w * w - p['dE']   # E tại θ=0, hạt đồng bộ
    Erot = 0.5 * p['I'] * wc**2
    print("=== B. Động lực học đầy đủ (đơn vị: m=ω0=ħ=ℓ0=1; 1 đơn vị năng lượng = ħω0 = 10 MeV) ===")
    print(f"I=χ={p['I']}, ΔE={p['dE']}, κ={p['kappa']:.6f}, R={R}, g={g}, T0=2π={2*math.pi:.4f}")
    print(f"góc giữ tĩnh max = {math.degrees(th):.2f}°;  θ̇_gãy chính xác = {wc:.5f} (tiệm cận {eps:.1f}, "
          f"đề xuất ≈3.2);  ½Iθ̇_gãy² = {Erot:.4f} = {Erot*HW0:.2f} MeV (tiệm cận 45 MeV)")
    tau_units = (5e-15 / (0.2 * C)) / (HBAR_MEV_S / HW0)   # τ_relax theo đơn vị 1/ω0
    windows = windows or [("τ_relax", tau_units), ("5T0", 5 * 2 * math.pi), ("10T0", 10 * 2 * math.pi)]
    print(f"Cửa sổ: " + ", ".join(f"{n} = {v:.3f}" for n, v in windows))
    pp = dict(p, R=R)
    print(f"{'θ̇0/θ̇_gãy':>10} {'E dư (MeV)':>10} {'e':>4} {'gãy':>4} {'va chạm':>8} {'t khóa đầu':>11} "
          f"{'t đồng bộ cuối':>15} {'/T0':>6} " + " ".join(f"{'≤'+n+'?':>10}" for n, _ in windows) + f" {'%E dư mất':>10}")
    for f in (1.05, 1.2, 1.5, 2.0):
        w0 = f * wc
        E0, Ec = E_at(w0), E_at(wc)
        for e in (0.0, 0.1, 0.2):
            r = simulate(pp, w0, t_end=10 * 2 * math.pi, e=e)
            ev = r['events']
            imp = [x for x in ev if x['kind'] != 'gãy']
            nb = sum(x['kind'] == 'gãy' for x in ev)
            t1 = next((x['t'] for x in ev if x['kind'] == 'khóa'), float('nan'))
            t_final = ev[-1]['t'] if (r['locked'] and ev and ev[-1]['kind'] == 'khóa') else float('nan')
            if r['locked'] and not ev:
                t_final = 0.0
            ok = [("có" if t_final <= v else "không") if not math.isnan(t_final) else "không" for _, v in windows]
            print(f"{f:10.2f} {(E0-Ec)*HW0:10.3f} {e:4.1f} {nb:4d} {len(imp):8d} {t1:11.4f} {t_final:15.4f} "
                  f"{t_final/(2*math.pi):6.2f} " + " ".join(f"{o:>10}" for o in ok)
                  + f" {(E0-r['E'])/(E0-Ec):10.1%}", flush=True)


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
