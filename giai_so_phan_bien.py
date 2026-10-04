"""Giải số cho bốn phương án trả lời phản biện của Part I.

(A) "Công nguồn lái cân bằng chính xác với tiêu tán Carnot": đo công lái W, tổng Carnot C, biến thiên năng lượng ΔE
    trong Langevin (τ=10, γx=0.7, f=1.05, 1 MeV); phần còn lại W − C − ΔE = (ma sát γx) − (năng lượng nhiễu đưa vào).
    Kèm tất định: ở trạng thái dừng không gãy, công lái chỉ đi vào ma sát.
(B) Hấp thụ mặt Koning–Delaroche (2003), neutron:
    W_D = d1 ΔE² e^{−d2 ΔE}/(ΔE² + d3²), d1 = 16(1 − (N−Z)/A), d2 = 0.0180 + 0.003802/(1+e^{(A−156)/8}), d3 = 11.5 MeV;
    W_V = w1 ΔE²/(ΔE² + w2²); γ_total = 2(W_V + W_D)/ħω0 với ħω0 = 10 MeV. ΔE = E − E_F, E_F = −11.2814 + 0.02646A.
Chạy: python3 giai_so_phan_bien.py [A|B]   (A ≈ 1 phút, 8 lõi)
"""
import math
import sys
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import cfg_for
from giai_so_kramers_A import setup
from giai_so_lien_ket_e_gay import energy_p

T0 = L.T0


def _job(a):
    cfg, kT, seed, T = a
    setup(cfg)
    wc = cfg[2]
    s0 = [0.0, 0.82 * wc, 0.0, L.BASE['kappa'] * 0.82 * wc]
    r = L.run(0.1, 1.05 * wc, kT, T, 2e-3, seed, gx=0.7, w_start=0.82 * wc)
    return r['W'], r['carnot'], r['E'] - energy_p(s0, L.BASE), sum(e[0] == 'gãy' for e in r['events'])


def part_a(N=8):
    cfg = cfg_for(0.1)
    T = 60 * T0
    print("=== A. Ngân sách công của nguồn lái (60 T0 mỗi quỹ đạo) ===")
    with Pool() as pool:
        for kT in (0.0, 0.1):
            res = pool.map(_job, [(cfg, kT, 41000 + i, T) for i in range(N if kT > 0 else 1)])
            W = sum(r[0] for r in res); C = sum(r[1] for r in res); dE = sum(r[2] for r in res)
            nb = sum(r[3] for r in res)
            print(f"k_BT = {kT*10:g} MeV: công lái W = {W:.2f}, Carnot C = {C:.2f} ({C/W:.1%} của W), ΔE hệ = {dE:+.2f}, "
                  f"W − C − ΔE = {W-C-dE:.2f} (ma sát − nhiễu đưa vào), số lần gãy = {nb}")
    print()


def part_b():
    print("=== B. Koning–Delaroche neutron: W_V + W_D và γ_total = 2(W_V+W_D)/ħω0 ===")
    for A, Z in ((56, 26), (70, 32), (208, 82)):
        N = A - Z
        w1, w2 = 12.195 + 0.0167 * A, 73.55 + 0.0795 * A
        d1 = 16.0 * (1 - (N - Z) / A)
        d2 = 0.0180 + 0.003802 / (1 + math.exp((A - 156.0) / 8.0))
        d3 = 11.5
        EF = -11.2814 + 0.02646 * A
        print(f"A={A} (Z={Z}): E_F = {EF:.2f} MeV, d1 = {d1:.2f}, d2 = {d2:.4f}")
        hit = None
        for dE in (2, 4, 5, 6, 7, 8, 10, 12, 15):
            WV = w1 * dE**2 / (dE**2 + w2**2)
            WD = d1 * dE**2 * math.exp(-d2 * dE) / (dE**2 + d3**2)
            g = 2 * (WV + WD) / 10.0
            if hit is None and g >= 0.7:
                hit = dE
            print(f"   ΔE = {dE:4.1f} MeV (E = {EF+dE:+6.2f} MeV): W_V = {WV:.3f}, W_D = {WD:.3f}, γ_total = {g:.3f}")
        print(f"   ⟹ γ_total ≥ 0.70 từ ΔE ≈ {hit} MeV")
    print("Lưu ý: W_D là ĐỘ SÂU của thừa số dạng đạo hàm Woods–Saxon (đỉnh ở bề mặt r ≈ R); với nucleon trong lòng "
          "hạt nhân, trung bình trên hàm sóng nhỏ hơn độ sâu này. E = E_F + ΔE < 0 là vùng liên kết, nằm ngoài miền "
          "tán xạ (E > 0) mà KD03 được khớp.\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["B", "A"]:
        {"A": part_a, "B": part_b}[name.upper()]()
