"""Giải số bộ thông số chuẩn hóa đề xuất: I = χħ/ω0 (cộng hưởng Ik = mΔE), L ≈ 1 fm (E_vt ≈ Λ_QCD),
cửa sổ tái đồng bộ T ~ (3…10)·T0.

(A) SI: I, θ̇_gãy/ω0, năng lượng quay tới hạn ½Iθ̇² ≈ E_vt²/(2ħω0), v/c của hạt — với ħω0 = 1 MeV và
    ħω0 thực của nucleon trong hạt nhân (≈ 41·A^(-1/3) MeV).
(B) Đơn vị m=ω0=ħ=ℓ0=1: so I = χ (cộng hưởng) với I = 1 — góc giữ tĩnh, θ̇_gãy chính xác, g = 0 và g = 0.1.
(C) Chu kỳ gãy/khóa ở đúng bộ thông số QCD cộng hưởng (χ=1, ε=197.327), cửa sổ T = 10·T0 = 20π.

Chạy: python3 giai_so_cong_huong_qcd.py [A|B|C]   (C ≈ vài phút)
"""
import math
import sys

from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params, static_hold_angle
from giai_so_va_cham_tai_lien_ket import simulate

HBAR = 1.054571817e-34
MEV = 1.602176634e-13
MP_MEV = 938.27208816
HBAR_C = 197.3269804  # MeV·fm


def part_a():
    print("=== A. Bộ thông số QCD + cộng hưởng, đơn vị thật (proton) ===")
    print(f"{'nguồn ħω0':>16} {'ħω0 MeV':>8} {'E_vt MeV':>9} {'χ':>4} {'I (kg·m²)':>11} {'kκ/R':>9} "
          f"{'θ̇_gãy/ω0':>10} {'E_quay tới hạn':>15} {'/m_p c²':>8} {'v/c lúc gãy':>11}")
    rows = [("giả thiết cũ", 1.0)] + [(f"41·A^-1/3, A={A}", 41.0 * A ** (-1 / 3)) for A in (16, 40, 208)]
    for Ev in (HBAR_C / 1.0, 200.0):
        for name, hw in rows:
            for chi in (1.0,):
                w0 = hw * MEV / HBAR
                I = chi * HBAR / w0
                ratio = chi * (hw / Ev)**2
                wg = Ev / (math.sqrt(chi) * hw)
                Erot = Ev**2 / (2 * hw)                      # ½Iθ̇² với θ̇ = θ̇_gãy (tiệm cận), MeV
                vc = math.sqrt(chi * hw / MP_MEV)            # v_max = κθ̇_gãy
                print(f"{name:>16} {hw:8.3f} {Ev:9.2f} {chi:4.1f} {I:11.3e} {ratio:9.2e} {wg:10.2f} "
                      f"{Erot/1000:12.3f} GeV {Erot/MP_MEV:8.2f} {vc:11.4f}")
    print("(E_quay tới hạn ≈ E_vt²/(2ħω0), không phụ thuộc χ)\n")


def part_b():
    print("=== B. Cộng hưởng I = χ so với I = 1 (m=ω0=ħ=ℓ0=1) ===")
    print(f"{'g':>4} {'χ':>5} {'ε/√χ':>6} {'I':>5} {'θ giữ tĩnh':>11} {'θ̇_gãy':>10} {'tiệm cận':>9} {'sai lệch':>9}")
    for g in (0.0, 0.1):
        for chi in (0.1, 1.0, 10.0):
            for f in (0.5, 1.0, 5.0):
                eps = f * math.sqrt(chi)
                for I in ((chi,) if chi == 1.0 else (chi, 1.0)):
                    p, R = params(chi, eps, I=I, g=g)
                    th = static_hold_angle(p, R)
                    wc = omega_crit_general(p, R)
                    dev = f"{(f - wc) / wc:+9.2%}" if wc > 0 else "        —"
                    tag = "=χ" if I == chi else "=1"
                    print(f"{g:4.1f} {chi:5.1f} {f:6.1f} {tag:>5} {math.degrees(th):10.2f}° {wc:10.5f} {f:9.5f} {dev}")
    print()


def part_c():
    chi, eps = 1.0, HBAR_C / 1.0
    p, R = params(chi, eps, I=chi, g=0.0)
    p = dict(p, R=R)
    wc = omega_crit_general({k: v for k, v in p.items() if k != 'R'}, R)
    T0 = 2 * math.pi
    E_crit = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * wc**2 - p['dE']
    print(f"=== C. Chu kỳ gãy/khóa ở bộ QCD cộng hưởng: χ={chi}, ε=E_vt/ħω0={eps:.3f}, I={p['I']}, "
          f"κ={p['kappa']:.6f}, R={R:.3f}, g=0 ===")
    print(f"θ̇_gãy chính xác = {wc:.4f} (tiệm cận {eps/math.sqrt(chi):.4f}); E tại ngưỡng = {E_crit:.3f}; "
          f"cửa sổ T = 10·T0 = {10*T0:.3f}")
    for e in (0.0, 0.3, 0.6):
        w_start = 1.05 * wc
        E0 = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * w_start**2 - p['dE']
        r = simulate(p, w_start, t_end=10 * T0, e=e)
        ev = r['events']
        nb = sum(x['kind'] == 'gãy' for x in ev)
        imp = [x for x in ev if x['kind'] != 'gãy']
        nk = sum(x['kind'] == 'khóa' for x in imp)
        lost = E0 - r['E']
        per = lost / max(len(imp), 1)
        need = E0 - E_crit
        print(f"e={e}: θ̇(0)={w_start:.3f} (1.05×), E0={E0:.3f}, gãy={nb}, va chạm={len(imp)} (khóa {nk}), "
              f"E mất={lost:.4f} ({lost/need:.2%} phần dư trên ngưỡng {need:.3f}), "
              f"mất TB/va chạm={per:.4f} ⟹ cần ≈{need/per if per>0 else float('inf'):.0f} va chạm, "
              f"cuối: {'ĐỒNG BỘ' if r['locked'] else 'TỰ DO'}, θ̇={r['s'][1]:.3f}", flush=True)
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B", "C"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
