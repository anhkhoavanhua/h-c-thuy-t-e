"""Giải số ba đề xuất:
(A) τ từ "tương tác mạnh" α_πNN ≈ 14 thay vào Γ/ω0 = (2/3)α ħω0/(m_p c²); ngưỡng phát pion thật m_π c²;
    so với hằng số ghép giả-vectơ f²/4π ≈ 0.075.
(B) Cơ chế Kramers: đo TỐC ĐỘ GÃY Ở TRẠNG THÁI DỪNG (bỏ giai đoạn chuyển tiếp) theo k_BT, và nhiệt độ
    Kiểm Arrhenius: ln r = const − ΔK/k_BT với ΔK = 1.24.
(C) "γx = const ≈ 0.70 ⟹ η ≈ 0.40 ⟹ ΔE = 16.4 A^-1/3": đo η thật với γx = 0.70, và với bộ tham số
    "cố định tiên nghiệm" nhất quán (E_vt = 30 MeV cố định ⟹ E_vt/ħω0 đổi theo A).

Chạy: python3 giai_so_kramers_A.py [A|B|C]   (B ≈ 8 phút, C ≈ 3 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_chi_va_dai_hoi_phuc import omega_crit_general, params

MPC2, MPIC2 = 938.27208816, 139.57039
GX = 0.70
TAU = 10.0


def part_a():
    print("=== A. τ từ tương tác mạnh ===")
    for name, a in (("α_πNN = g²/4π ≈ 14 (giả-vô hướng)", 14.0), ("f²/4π ≈ 0.075 (giả-vectơ)", 0.075)):
        r = (2 / 3) * a * 10.0 / MPC2
        print(f"{name}: Γ/ω0 = {r:.4f} ⟹ τ = {1/r:.2f}")
    print(f"Ngưỡng phát pion thật: m_π c² = {MPIC2:.2f} MeV ≫ ħω0 = 10 MeV ⟹ không phát được pion thật (tỉ số {MPIC2/10:.1f})")
    print(f"Quan hệ g²/4π = (2m_p/m_π)² f²/4π: (2·938.27/139.57)² × 0.075 = {(2*MPC2/MPIC2)**2*0.075:.2f}\n")


def setup(cfg):
    """Đặt tham số toàn cục của mô-đun L cho một cấu hình (chạy trong từng tiến trình con)."""
    base, R, wc = cfg
    L.BASE, L.R, L.WC = base, R, wc


def _stationary(a):
    cfg, kT, seed, T, warm, dt = a
    setup(cfg)
    base, R, wc = cfg
    w_ss = wc * 0.82                                  # khởi động gần vận tốc dừng dưới ngưỡng
    r = L.run(1.0 / TAU, 1.05 * wc, kT, T, dt, seed, gx=GX, w_start=w_ss, sample_every=None)
    br = [e[1] for e in r['events'] if e[0] == 'gãy' and e[1] > warm]
    lk = [e[1] for e in r['events'] if e[0] == 'khóa']
    # "đợt thoát": lần gãy xảy ra sau khi đã khóa ổn định ≥ T0 (loại các lần gãy lại ngay)
    episodes, last_lock = 0, None
    for e in r['events']:
        if e[0] == 'khóa':
            last_lock = e[1]
        elif e[0] == 'gãy' and e[1] > warm and (last_lock is None or e[1] - last_lock >= L.T0):
            episodes += 1
    return len(br), episodes


def std_cfg():
    base, R = params(1.0, 3.0, I=1.0, g=0.1)
    return base, R, omega_crit_general(base, R)


def part_b(N=16):
    cfg = std_cfg()
    base, R, wc = cfg
    kap = base['kappa']
    w_dung = 1.05 * wc * math.sqrt(1 - GX * kap * kap / (2 / TAU))
    dK = 0.5 * base['I'] * (wc**2 - w_dung**2)
    print(f"=== B. Kramers: tốc độ gãy dừng theo k_BT (τ={TAU}, γx={GX}, θ̇_gãy={wc:.4f}, "
          f"θ̇_dừng ước lượng={w_dung:.4f}, ΔK={dK:.3f}) ===")
    T, warm, dt = 100 * L.T0, 10 * L.T0, 2e-3
    rows = []
    with Pool() as pool:
        for kT in (0.1, 0.2, 0.3, 0.5, 0.8):
            res = pool.map(_stationary, [(cfg, kT, 9000 + i, T, warm, dt) for i in range(N)])
            nb = sum(r[0] for r in res)
            ne = sum(r[1] for r in res)
            span = N * (T - warm)
            rate_e = ne / span
            rows.append((kT, rate_e))
            print(f"k_BT={kT:.1f} ({kT*10:g} MeV): gãy (sau khởi động) = {nb}, đợt thoát = {ne} trên "
                  f"{span/L.T0:.0f} T0 ⟹ r = {rate_e:.5f} /đv t = {rate_e*L.T0:.4f} /T0", flush=True)
    pts = [(1 / k, math.log(r)) for k, r in rows if r > 0]
    if len(pts) >= 2:
        n = len(pts)
        mx, my = sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n
        slope = sum((p[0] - mx) * (p[1] - my) for p in pts) / sum((p[0] - mx)**2 for p in pts)
        print(f"Khớp Arrhenius ln r = a + b/k_BT trên {n} điểm có r>0: b = {slope:.3f} (Kramers dự đoán b = −ΔK = {-dK:.3f})")
    print(f"Kramers với prefactor bất kỳ: r(1 MeV)/r(5 MeV) = exp(−ΔK(1/0.1 − 1/0.5)) = {math.exp(-dK*(10-2)):.2e}\n")


def _eta(a):
    cfg, kT, seed = a
    setup(cfg)
    base, R, wc = cfg
    r = L.run(1.0 / TAU, 1.05 * wc, kT, 10 * L.T0, 2e-3, seed, gx=GX)
    return r['carnot'], sum(e[0] == 'khóa' for e in r['events'])


def part_c(N=32):
    print(f"=== C. η với γx = {GX} (τ = {TAU}, nhiễu 1 MeV, N = {N}) ===")
    nuclei = [("¹⁶O", 16, 6.130), ("⁴⁰Ca", 40, 3.737), ("²⁰⁸Pb", 208, 2.614)]
    print(f"{'hạt nhân':>8} {'giả thiết':>28} {'E_vt/ħω0':>9} {'η':>7} {'ΔE mô hình':>11} {'16.4A^-1/3':>11} {'3⁻':>7} {'lệch':>7}")
    with Pool() as pool:
        for name, A, e3 in nuclei:
            hw = 41 * A ** (-1 / 3)
            for label, ratio in (("tỉ số 3 cố định (đề xuất)", 3.0), ("E_vt = 30 MeV cố định", 30.0 / hw)):
                base, R = params(1.0, ratio, I=1.0, g=0.1)
                cfg = (base, R, omega_crit_general(base, R))
                res = pool.map(_eta, [(cfg, 1.0 / hw, 7000 + i) for i in range(N)])
                eta = sum(r[0] for r in res) / max(sum(r[1] for r in res), 1)
                dE = eta * hw
                print(f"{name:>8} {label:>28} {ratio:9.3f} {eta:7.4f} {dE:8.3f} MeV {16.4*A**(-1/3):8.3f} MeV "
                      f"{e3:7.3f} {(dE-e3)/e3:+7.1%}", flush=True)
    print(f"Đề xuất ghi ⁴⁰Ca 'lệch +20%': 16.4·40^-1/3 = {16.4*40**(-1/3):.3f} MeV so với 3.737 ⟹ {(16.4*40**(-1/3)-3.737)/3.737:+.1%}\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "C", "B"]:
        {"A": part_a, "B": part_b, "C": part_c}[name.upper()]()
