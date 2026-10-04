"""Giải số hai đề xuất mới:
(A) κ = ℓ0·E_kép/E_vt với E_vt = ħc/L — tính bằng đơn vị SI thật cho proton, kiểm các điều kiện
    κ ≤ ℓ0 (đề xuất tự đặt ra), liên kết giữ được ở trạng thái nghỉ (kκ < R), và θ̇_gãy/ω0.
(B) Hệ số hồi phục e ∈ [0,1] khi tái liên kết, với E_kép = ΔE (Rκ = ΔE): số va chạm, năng lượng
    mất, trạng thái cuối — kiểm nhận định "e > 0 làm hạt nảy hỗn loạn, phá đồng bộ".

Chạy: python3 giai_so_thang_vu_tru_va_hoi_phuc.py   (chỉ thư viện chuẩn)
"""
import math

from giai_so_va_cham_tai_lien_ket import omega_crit, simulate

HBAR = 1.054571817e-34      # J·s
C = 2.99792458e8            # m/s
EV = 1.602176634e-19        # J
MP = 1.67262192e-27         # kg (proton)
MPC = 3.0856775814913673e22  # m
H0 = 67.4e3 / MPC           # 1/s (Planck 2018)


def part_a(E_kep_eV=1e6, hw0_eV=1e6):
    w0 = hw0_eV * EV / HBAR
    l0 = math.sqrt(HBAR / (MP * w0))
    k = MP * w0**2
    Ek = E_kep_eV * EV
    print(f"=== A. Thang vi mô–vĩ mô: proton, ħω0 = {hw0_eV/1e6:g} MeV, E_kép = {E_kep_eV/1e6:g} MeV ===")
    print(f"ω0 = {w0:.4e} rad/s, ℓ0 = {l0:.4e} m, k = mω0² = {k:.4e} N/m")
    print(f"Điều kiện κ ≤ ℓ0 ⟺ E_vt ≥ E_kép ⟺ L ≤ ħc/E_kép = {HBAR*C/Ek:.4e} m")
    cases = [
        ("Bán kính Hubble c/H0", C / H0),
        ("Thang năng lượng tối ρΛ^(1/4)≈2.24 meV", HBAR * C / (2.24e-3 * EV)),
        ("L = 1 m", 1.0),
        ("L = ħc/E_kép (biên)", HBAR * C / Ek),
        ("L = 1 fm", 1e-15),
    ]
    print(f"{'Trường hợp':42} {'L (m)':>11} {'E_vt (eV)':>11} {'κ (m)':>11} {'κ/ℓ0':>10} "
          f"{'R (N)':>11} {'kκ/R':>10} {'θ̇_gãy/ω0':>11}")
    for name, L in cases:
        Ev = HBAR * C / L
        kappa = l0 * Ek / Ev
        kappa_formula = L * Ek / (C * math.sqrt(HBAR * MP * w0))   # dạng gộp trong đề xuất
        assert abs(kappa / kappa_formula - 1) < 1e-9
        R = Ev / l0
        ratio_static = k * kappa / R                                 # >1: gãy ngay khi đứng yên
        wg = Ev / math.sqrt(HBAR * w0 * Ek)                           # tiệm cận θ̇_gãy/ω0
        print(f"{name:42} {L:11.3e} {Ev/EV:11.3e} {kappa:11.3e} {kappa/l0:10.2e} "
              f"{R:11.3e} {ratio_static:10.2e} {wg:11.3e}")
    print("(kκ/R = ħω0·E_kép/E_vt²; kκ/R > 1 ⟹ lực đàn hồi tĩnh đã vượt ngưỡng: không thể đồng bộ)\n")


def part_b():
    p = dict(I=1.0, m=1.0, k=4.0, g=0.1, dE=0.5, kappa=0.2, R=2.5)   # Rκ = 0.5 = ΔE
    wc = omega_crit(p)
    E0 = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * 25.0 - p['dE']
    print(f"=== B. Quét hệ số hồi phục e: {', '.join(f'{k}={v}' for k, v in p.items())} (Rκ = {p['R']*p['kappa']} = ΔE) ===")
    print(f"θ'(0)=5, E0={E0:.6f}, θ̇ tới hạn chính xác = {wc:.5f}, t∈[0,200], khóa khi |u| < 1e-3")
    print(f"{'e':>4} {'gãy':>4} {'nảy':>5} {'khóa':>5} {'E mất':>10} {'%E0':>7} {'Σ Carnot':>10} "
          f"{'max|lệch|':>10} {'t khóa cuối':>11} {'trạng thái cuối':>16} {'θ̇ cuối':>9}")
    for e in (0.0, 0.3, 0.6, 0.9, 1.0):
        r = simulate(p, 5.0, t_end=200.0, e=e)
        ev = r['events']
        nb = sum(x['kind'] == 'gãy' for x in ev)
        nn = sum(x['kind'] == 'nảy' for x in ev)
        nk = sum(x['kind'] == 'khóa' for x in ev)
        imp = [x for x in ev if x['kind'] in ('nảy', 'khóa')]
        lost = sum(x['E_before'] - x['E_after'] for x in imp)
        carnot = sum(x['carnot'] for x in imp)
        dev = max((abs(x['E_before'] - x['E_after'] - x['carnot']) for x in imp), default=0.0)
        tl = max((x['t'] for x in ev if x['kind'] == 'khóa'), default=float('nan'))
        print(f"{e:4.1f} {nb:4d} {nn:5d} {nk:5d} {lost:10.6f} {lost/E0:7.2%} {carnot:10.6f} {dev:10.1e} "
              f"{tl:11.4f} {'ĐỒNG BỘ' if r['locked'] else 'TỰ DO':>16} {r['s'][1]:9.5f}")


if __name__ == "__main__":
    part_a()
    part_b()
