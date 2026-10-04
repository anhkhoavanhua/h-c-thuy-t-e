"""Giải số ba đề xuất:
(A) L cỡ hạt nhân (1–200 fm) và tỉ số χ = E_kép/ħω0 — bảng SI cho proton.
(B) Kiểm công thức χ trên toàn bộ động lực học (đơn vị không thứ nguyên m = ω0 = ħ = ℓ0 = 1):
    biên k·κ/R = χ(ħω0/E_vt)² < 1, góc giữ tĩnh lớn nhất, θ̇_gãy chính xác so với tiệm cận E_vt/(√χ ħω0).
(C) Quét mịn hệ số hồi phục e = 0…1 để tìm dải e thực sự tái đồng bộ (khóa lại) trong t ≤ 200.

Chạy: python3 giai_so_chi_va_dai_hoi_phuc.py [A|B|C]   (mặc định chạy cả ba; C ≈ 15–20 phút)
"""
import math
import sys

from giai_so_lien_ket_e_gay import lam_closed, lam_rot_max
from giai_so_va_cham_tai_lien_ket import simulate

HBAR_C_MEV_FM = 197.3269804  # ħc tính bằng MeV·fm


def part_a(hw0=1.0):
    print(f"=== A. L cỡ hạt nhân, proton ħω0 = {hw0} MeV, E_vt = ħc/L ===")
    print(f"{'L (fm)':>8} {'E_vt (MeV)':>11} {'χ':>6} {'E_kép':>7} {'κ/ℓ0':>8} {'kκ/R':>9} "
          f"{'biên √χ·ħω0':>12} {'giữ được?':>10} {'θ̇_gãy/ω0':>10}")
    for L in (0.2, 1.0, 10.0, 100.0, 197.3269804, 200.0):
        Ev = HBAR_C_MEV_FM / L
        for chi in (0.1, 1.0, 10.0):
            Ek = chi * hw0
            ratio = chi * (hw0 / Ev)**2
            print(f"{L:8.2f} {Ev:11.4f} {chi:6.1f} {Ek:7.2f} {Ek/Ev:8.4f} {ratio:9.3e} "
                  f"{math.sqrt(chi)*hw0:12.4f} {'có' if ratio < 1 else 'KHÔNG':>10} {Ev/(math.sqrt(chi)*hw0):10.4f}")
    print()


def params(chi, eps, I=1.0, g=0.0):
    """Đơn vị m=ω0=ħ=ℓ0=1 ⟹ k=1; E_kép = ΔE = χ; E_vt = ε; κ = ℓ0 E_kép/E_vt; R = E_vt/ℓ0."""
    return dict(I=I, m=1.0, k=1.0, g=g, dE=chi, kappa=chi / eps), eps


def static_hold_angle(p, R, n=20000):
    """Góc |θ| lớn nhất (0..π) mà hệ ĐỨNG YÊN vẫn có |λ| < R, quét từ θ=0 ra ngoài."""
    for i in range(1, n + 1):
        th = math.pi * i / n
        if abs(lam_closed(th, 0.0, **p)) >= R:
            return math.pi * (i - 1) / n
    return math.pi


def omega_crit_general(p, R):
    """θ'(0) nhỏ nhất (xuất phát θ=0) để max|λ| trên quỹ đạo đồng bộ chạm R; 0 nếu chưa quay đã gãy."""
    if lam_rot_max(1e-9, p, n=4000) >= R:
        return 0.0
    hi = 1.0
    while lam_rot_max(hi, p, n=4000) < R:
        hi *= 2.0
    lo = 0.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if lam_rot_max(mid, p, n=4000) >= R:
            hi = mid
        else:
            lo = mid
    return hi


def part_b():
    print("=== B. Kiểm công thức χ trên động lực học đầy đủ (m=ω0=ħ=ℓ0=1, I=1, g=0) ===")
    print(f"{'χ':>6} {'ε=E_vt/ħω0':>11} {'ε/√χ':>7} {'kκ/R':>10} {'θ giữ tĩnh max':>15} "
          f"{'θ̇_gãy chính xác':>17} {'tiệm cận ε/√χ':>14} {'sai lệch':>9}")
    for chi in (0.1, 1.0, 10.0):
        for f in (0.5, 1.0, 2.0, 5.0, 20.0):
            eps = f * math.sqrt(chi)
            p, R = params(chi, eps)
            ratio = p['k'] * p['kappa'] / R
            th_s = static_hold_angle(p, R)
            wc = omega_crit_general(p, R)
            asym = eps / math.sqrt(chi)
            dev = f"{(asym - wc) / wc:+9.2%}" if wc > 0 else "      —"
            print(f"{chi:6.1f} {eps:11.4f} {f:7.2f} {ratio:10.3e} {math.degrees(th_s):13.2f}° "
                  f"{wc:17.5f} {asym:14.5f} {dev}")
    print()


def part_c():
    p = dict(I=1.0, m=1.0, k=4.0, g=0.1, dE=0.5, kappa=0.2, R=2.5)   # Rκ = ΔE (như bảng 2.2)
    print(f"=== C. Quét mịn e: {', '.join(f'{k}={v}' for k, v in p.items())}, t∈[0,200], khóa khi |u| < 1e-3 ===")
    for w0 in (5.0, 6.0):
        E0 = 0.5 * (p['I'] + p['m'] * p['kappa']**2) * w0**2 - p['dE']
        print(f"-- θ'(0) = {w0}, E0 = {E0:.6f}")
        print(f"{'e':>4} {'gãy':>4} {'nảy':>5} {'khóa':>5} {'%E0 mất':>8} {'t khóa đầu':>11} {'cuối':>8} {'θ̇ cuối':>9}")
        for e in [i / 10 for i in range(11)]:
            r = simulate(p, w0, t_end=200.0, e=e)
            ev = r['events']
            nb = sum(x['kind'] == 'gãy' for x in ev)
            nn = sum(x['kind'] == 'nảy' for x in ev)
            nk = sum(x['kind'] == 'khóa' for x in ev)
            t1 = next((x['t'] for x in ev if x['kind'] == 'khóa'), float('nan'))
            print(f"{e:4.1f} {nb:4d} {nn:5d} {nk:5d} {(E0 - r['E']) / E0:8.2%} {t1:11.4f} "
                  f"{'ĐỒNG BỘ' if r['locked'] else 'TỰ DO':>8} {r['s'][1]:9.5f}", flush=True)
    print()


if __name__ == "__main__":
    parts = sys.argv[1:] or ["A", "B", "C"]
    for name in parts:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
