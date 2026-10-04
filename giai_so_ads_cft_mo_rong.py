"""AdS/CFT mở rộng cho E_dual (HocThuyetE (7).pdf, Bảng I). Nền: L_AdS = ħc/ħω0 = 19.73 fm (giai_so_ads_cft.py).
Đơn vị không thứ nguyên của bản thảo: ħ = m = ω0 = 1, ℓ0 = 1, nên L = mc/(ħ/ℓ0)… = c/(ω0 ℓ0) = 9.690 ℓ0.

(F) Ngưỡng gãy R = E_vt/ℓ0 đọc theo AdS: lực hấp dẫn AdS (NR) tại bán kính r là k r ⟹ R = k·r_br.
    Tính lại θ̇_br từ Eq. (7) của bản thảo, rồi thêm hiệu chỉnh AdS chính xác cho thế tĩnh
    V'(x) = k x/√(1 + x²/L²) và ước lượng khối lượng tương đối tính γ(v) m.
(G) Độ sâu xuyên tâm (UV/IR): năng lượng dao động E ↔ bán kính quay đầu r_max = ℓ0 √(2E/ħω0) trong AdS.
(H) Hawking–Page trong AdS4: T(r_h) = (1 + 3r_h²/L²)/(4π r_h), F ∝ r_h(1 − r_h²/L²) ⟹ T_min, T_HP (giải số),
    đổi ra MeV. Rồi quét k_BT bằng mô phỏng Langevin Bảng I (hàm _traj của giai_so_dong_ho_photon.py, không sửa):
    tỉ lệ đợt đồng hồ đọc "Tương lai trước Hiện tại", tỉ lệ bước Φ tăng, sai số t̂.
(K) Kiểm trùng hợp: đổi γx hoặc ΔE (T_HP không đổi) rồi quét lại k_BT.
Chạy: python3 giai_so_ads_cft_mo_rong.py [F|G|H|K]   (F, G tức thì; H ≈ 4 phút, K ≈ 3 phút trên 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

HBARC, MN, HW = 197.3270, 938.919, 10.0
L_FM = HBARC / HW
L0_FM = math.sqrt(HBARC**2 / (MN * HW))
L = L_FM / L0_FM                                  # L theo đơn vị ℓ0
I_, M_, K_, DE, G_, KAP, R_ = 1.0, 1.0, 1.0, 1.0, 0.1, 1 / 3, 3.0


def lam(th, w, ads=False, rel=False):
    """Phản lực ràng buộc Eq. (7); ads: thế tĩnh AdS chính xác kx → V'(x); rel: quán tính dọc γ³m, v = κθ̇cosθ."""
    x = KAP * np.sin(th)
    m = M_ / (1 - (KAP * w * np.cos(th) / L)**2)**1.5 if rel else M_
    kx = K_ * x / np.sqrt(1 + (x / L)**2) if ads else K_ * x
    num = np.sin(th) * (-I_ * m * KAP * w**2 - DE * m * KAP * np.cos(th) + G_ * m * KAP**2 * np.cos(th)**2) \
        + I_ * (kx - G_ * np.sin(th))
    return num / (I_ + m * KAP**2 * np.cos(th)**2)


def energy(th, w, ads, rel):
    """Năng lượng trên đa tạp ràng buộc; x = κ sinθ, v = κθ̇ cosθ, c = L (đơn vị ℓ0 ω0)."""
    x, v = KAP * np.sin(th), KAP * w * np.cos(th)
    kin = M_ * L**2 * (1 / np.sqrt(1 - (v / L)**2) - 1) if rel else 0.5 * M_ * v**2
    pot = K_ * L**2 * (np.sqrt(1 + (x / L)**2) - 1) if ads else 0.5 * K_ * x**2
    return 0.5 * I_ * w**2 + kin - DE * np.cos(th) + pot - G_ * x * np.sin(th)


def lam_orbit_max(w0, ads=False, rel=False, n=4001):
    """max|λ| trên một vòng; θ̇(θ) lấy từ bảo toàn năng lượng với θ̇(0) = w0 (như cách tính 3.205)."""
    th = np.linspace(0, 2 * math.pi, n)
    E = energy(0.0, w0, ads, rel)
    lo, hi = np.zeros(n), np.full(n, w0 * 3)
    for _ in range(60):                                # giải θ̇(θ) ≥ 0 bằng chia đôi (năng lượng tăng theo θ̇)
        mid = (lo + hi) / 2
        big = energy(th, mid, ads, rel) > E
        hi, lo = np.where(big, mid, hi), np.where(big, lo, mid)
    return np.max(np.abs(lam(th, hi, ads, rel)))


def wbr(**kw):
    lo, hi = 2.0, 5.0
    for _ in range(50):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if lam_orbit_max(mid, **kw) < R_ else (lo, mid)
    return hi


def part_f():
    print("=== F. Ngưỡng gãy R là thoát khỏi giếng AdS ===")
    F1 = HW / L0_FM                                   # k ℓ0 = ħω0/ℓ0 (MeV/fm)
    print(f"lực hấp dẫn AdS tại r = ℓ0: k ℓ0 = ħω0/ℓ0 = {F1:.3f} MeV/fm = {F1*160.2177:.1f} N")
    print(f"R = E_vt/ℓ0 = {30/L0_FM*160.2177:.1f} N = {30/HW:.1f} × k ℓ0 ⟹ r_br = 3ℓ0 = 9κ = {3*L0_FM:.3f} fm "
          f"= {3*L0_FM/L_FM:.3f} L = {3*L0_FM/(1.2*70**(1/3)):.2f} R_hạt nhân")
    b0, b1, b2 = wbr(), wbr(ads=True), wbr(ads=True, rel=True)
    print(f"θ̇_br (θ̇ tại θ = 0, θ̇(θ) theo bảo toàn năng lượng): {b0:.4f} ω0 (bản thảo: 3.205; tiệm cận √(R/mκ) = {math.sqrt(R_/KAP):.3f})")
    print(f"  + thế tĩnh AdS chính xác:        {b1:.4f} ω0 (Δ = {b1-b0:+.1e})")
    print(f"  + động học tương đối tính (γ³m):  {b2:.4f} ω0 (Δ = {b2-b0:+.1e}, v_max/c = {KAP*b2/L:.4f})\n")


def part_g():
    print("=== G. Độ sâu xuyên tâm của các mức năng lượng (bán kính quay đầu trong AdS) ===")
    print(f"{'trạng thái':>34} {'E (MeV)':>8} {'r_max (fm)':>11} {'r_max/L':>8}")
    rows = [("điểm không ħω0/2", 5.0), ("hạt bị kéo ở ngưỡng ½m(κθ̇_br)²", 0.5 * (KAP * 3.205)**2 * HW),
            ("một va chạm Carnot ⟨ΔE⟩", 4.3), ("Sn (phát nucleon)", 8.0), ("trần đợt (bảng V.B.3)", 20.4),
            ("E_vt", 30.0)]
    for name, E in rows:
        r = L0_FM * math.sqrt(2 * E / HW)
        print(f"{name:>34} {E:8.2f} {r:11.3f} {r/L_FM:8.3f}")
    print(f"(R_hạt nhân = {1.2*70**(1/3):.2f} fm; mọi mức nằm ở r ≪ L: vùng hồng ngoại sâu, AdS gần như phẳng + thế điều hòa)\n")


def part_h_theory():
    rh = np.linspace(0.05, 5, 200001)                  # theo đơn vị L
    T = (1 + 3 * rh**2) / (4 * math.pi * rh)
    i = np.argmin(T)
    Fz = rh[np.argmin(np.abs(rh * (1 - rh**2)))]
    j = np.argmin(np.abs(rh - Fz))
    print("=== H. Hawking–Page trong AdS4 (T theo đơn vị ħc/L = ħω0) ===")
    print(f"T_min = {T[i]:.5f} (giải tích √3/2π = {math.sqrt(3)/(2*math.pi):.5f}) = {T[i]*HW:.3f} MeV tại r_h = {rh[i]:.4f} L")
    print(f"T_HP  = {T[j]:.5f} (giải tích 1/π = {1/math.pi:.5f}) = {T[j]*HW:.3f} MeV tại r_h = {rh[j]:.4f} L")
    print(f"⟹ k_BT/ħω0 < {T[i]:.3f}: chỉ có AdS nhiệt (pha 'giam'); > {T[j]:.3f}: lỗ đen AdS thắng (pha 'giải giam')\n")


def _analyse(res, T, warm, every):
    from giai_so_dong_ho_photon import phi_at
    T0 = 2 * math.pi
    nb = rev = 0
    up, errs = [], []
    for samp, ev in res:
        for i, e in enumerate(ev):
            if e[0] != 'gãy' or e[1] < warm:
                continue
            prev = next((x for x in reversed(ev[:i]) if x[0] == 'khóa'), None)
            if prev is None or e[1] - prev[1] < T0:
                continue
            end = None
            for j in range(i + 1, len(ev)):
                if ev[j][0] == 'khóa':
                    nxt = next((y for y in ev[j + 1:] if y[0] == 'gãy'), None)
                    if nxt is None or nxt[1] - ev[j][1] >= T0:
                        end = ev[j][1] if nxt is not None else None
                        break
            if end is None or end > T - every:
                continue
            nb += 1
            rev += phi_at(samp, every, end) - phi_at(samp, every, e[1]) <= 0
        k0 = int(warm / every)
        ph = samp[k0:] - samp[k0]
        up.append(float(np.mean(np.diff(ph) > 0)))
        tt = np.arange(len(ph)) * every
        rate = ph[-1] / tt[-1]
        if abs(rate) > 1e-9:
            errs.append(float(np.sqrt(np.mean(((ph / rate - tt) / T0)**2))))
    return nb, rev, float(np.mean(up)), float(np.median(errs)) if errs else float('nan')


def part_h_sim(N=8):
    from giai_so_dong_ho_photon import _traj
    from giai_so_burst_g import cfg_for
    T0 = 2 * math.pi
    cfg = cfg_for(0.1)
    T, warm, every = 50 * T0, 5 * T0, 0.05
    print(f"=== H'. Quét k_BT bằng mô phỏng Langevin Bảng I ({N} quỹ đạo × 50 T0 mỗi mức) ===")
    print(f"{'k_BT (MeV)':>10} {'k_BT/T_HP':>10} {'đợt':>5} {'Tương lai TRƯỚC Hiện tại':>25} {'bước Φ tăng':>12} "
          f"{'sai số t̂ trung vị (T0)':>23}")
    with Pool() as pool:
        for kT in (0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5):
            res = pool.map(_traj, [(cfg, kT, 71000 + i, T, every) for i in range(N)])
            nb, rev, up, err = _analyse(res, T, warm, every)
            print(f"{kT*10:10.1f} {kT*math.pi:10.2f} {nb:5d} {rev:>12d}/{nb:<4d} ({rev/max(nb,1):5.1%}) {up:12.1%} "
                  f"{err:23.2f}", flush=True)
    print()


def _traj_k(a):
    import giai_so_dong_co_langevin as Lm
    from giai_so_kramers_A import setup
    cfg, kT, seed, T, every, gx = a
    setup(cfg)
    r = Lm.run(0.1, 1.05 * cfg[2], kT, T, 2e-3, seed, gx=gx, w_start=0.82 * cfg[2], sample_every=every)
    return np.array(r['samples']), r['events']


def part_k(N=8):
    """Kiểm trùng hợp: đổi tham số KHÔNG đụng ω0 (T_HP giữ nguyên 3.18 MeV), xem điểm vỡ thứ tự có dịch không."""
    from giai_so_chi_va_dai_hoi_phuc import params, omega_crit_general
    T0 = 2 * math.pi
    T, warm, every = 50 * T0, 5 * T0, 0.05
    cases = [("γx = 0.35 (Bảng I còn lại)", 1.0, 3.0, 0.35), ("ΔE = 0.5, E_vt = 1.5 (κ giữ ℓ0/3)", 0.5, 1.5, 0.7)]
    print(f"=== K. Kiểm trùng hợp với T_HP ({N} quỹ đạo × 50 T0 mỗi mức; T_HP = 3.18 MeV không đổi) ===")
    with Pool() as pool:
        for name, dE, evt, gx in cases:
            base, R = params(dE, evt, I=1.0, g=0.1)
            cfg = (base, R, omega_crit_general(base, R))
            print(f"-- {name}: θ̇_br = {cfg[2]:.3f}")
            for kT in (0.1, 0.2, 0.3, 0.4):
                res = pool.map(_traj_k, [(cfg, kT, 81000 + i, T, every, gx) for i in range(N)])
                nb, rev, up, err = _analyse(res, T, warm, every)
                print(f"{kT*10:10.1f} MeV {nb:5d} đợt, Tương lai TRƯỚC Hiện tại {rev:>3d}/{nb:<3d} ({rev/max(nb,1):5.1%}), "
                      f"bước Φ tăng {up:6.1%}, sai số t̂ trung vị {err:6.2f} T0", flush=True)
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["F", "G", "H"]:
        name = name.upper()
        if name == "F":
            part_f()
        elif name == "G":
            part_g()
        elif name == "K":
            part_k()
        else:
            part_h_theory(); part_h_sim()
