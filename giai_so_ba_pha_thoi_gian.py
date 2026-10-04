"""Giải số cho bảng ánh xạ ba pha E_quá khứ → E_hiện tại → E_tương lai và ba đoạn văn đề xuất cho bài báo.

(A) Ngắt lái trong toàn đợt (Q_θ = 0 từ lần gãy tới khi một lần khóa kéo dài ≥ T0), có/không tắt nhiễu trong đợt:
    tổng Carnot mỗi đợt, ngân sách quay ½I(θ̇_gãy² − θ̇_cuối²), tỉ lệ đợt ≤ 25 MeV.
(B) Pha E_tương lai: θ̇ ngay sau lần khóa kết thúc đợt có < θ̇_gãy không.
(C) Xung Carnot so với lượng tử ħω0: P(xung ≥ 1ħω0), P(≥ 2ħω0); bước mở kênh k, E_max, ⟨ε⟩ Ericson 2E/(n+1)
    với n = 2k ± 1 so với n = 3 của đề xuất.

Script này KHÔNG sửa giai_so_dong_co_langevin.py; nó dùng vòng lặp riêng với các hàm có sẵn trong đó.
Đơn vị m=ω0=ħ=ℓ0=1, ħω0 = 10 MeV, f = 1.05, τ = 10, γx = 0.70, g = 0.1, S_n = 8 MeV.
Chạy: python3 giai_so_ba_pha_thoi_gian.py   (≈ 8 phút, 8 lõi)
"""
import math
import random
from multiprocessing import Pool

import giai_so_dong_co_langevin as L
from giai_so_burst_g import cfg_for
from giai_so_kramers_A import setup
from giai_so_lien_ket_e_gay import rk4_step
from giai_so_va_cham_tai_lien_ket import impulse

TAU, GX, HW0, SN = 10.0, 0.70, 10.0, 8.0
T0 = L.T0


def run_burst(cfg, kT, seed, T, warm, mode, dt=2e-3):
    """Như L.run, thêm chế độ ngắt theo đợt. mode ∈ {None, 'drive', 'drive+noise'}. Trả về danh sách đợt."""
    setup(cfg)
    B, R, wc = L.BASE, L.R, cfg[2]
    gam, Om = 1 / TAU, 1.05 * wc
    rng = random.Random(seed)
    sig, sigx = math.sqrt(2 * gam * kT), math.sqrt(2 * GX * kT)
    s = [0.0, 0.82 * wc, 0.0, B['kappa'] * 0.82 * wc]
    t, locked, ref, in_burst, t_lock = 0.0, True, 0.0, False, 0.0
    bursts, cur = [], None
    while t < T:
        if in_burst and locked and t - t_lock >= T0:            # đợt khép lại: E_tương lai
            in_burst = False
            if cur is not None:
                bursts.append(cur)
            cur = None
        quiet = mode is not None and in_burst
        g = 0.0 if quiet else gam
        sc = lambda q, h: rk4_step(lambda z: L.f_con(z, g, Om, GX), q, h)
        sf = lambda q, h: rk4_step(lambda z: L.f_free(z, g, Om, GX), q, h)
        hit = lambda q: abs(L.acc_con(q, g, Om, GX)[2]) >= R
        h = dt
        if locked:
            if hit(s):                                           # gãy
                if not in_burst:
                    in_burst = True
                    cur = dict(w_break=s[1], pulses=[], w_end=None) if t > warm else None
                locked, ref = False, 0.0
                continue
            s1 = sc(s, dt)
            if hit(s1):
                h, s1 = L.locate(sc, s, dt, hit)
        else:
            s1 = sf(s, dt)
            f1 = L.fval(s1)
            if ref == 0.0:
                ref = math.copysign(1.0, f1)
            elif f1 * ref <= 0:
                h, s1 = L.locate(sf, s, dt, lambda q: L.fval(q) * ref <= 0)
                s1, Lam, u = impulse(s1, B)
                c = math.cos(s1[0])
                meff = B['m'] * B['I'] / (B['I'] + B['m'] * B['kappa']**2 * c * c)
                if cur is not None:
                    cur['pulses'].append(0.5 * meff * u * u * HW0)
                    cur['w_end'] = s1[1]
                locked, t_lock = True, t + h
        s, t = list(s1), t + h
        if not (quiet and mode == 'drive+noise'):
            k1 = sig * math.sqrt(h) * rng.gauss(0, 1)
            kx = sigx * math.sqrt(h) * rng.gauss(0, 1)
            if locked:
                c = math.cos(s[0])
                s[1] += (k1 + B['kappa'] * c * kx) / (B['I'] + B['m'] * B['kappa']**2 * c * c)
                s[3] = B['kappa'] * c * s[1]
            else:
                s[1] += k1 / B['I']
                s[3] += kx / B['m']
    return [b for b in bursts if b['w_end'] is not None]


def _job(a):
    return run_burst(*a)


def pct(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))]


def main(N=12):
    cfg = cfg_for(0.1)
    wc = cfg[2]
    T, warm = 70 * T0, 10 * T0
    print(f"θ̇_gãy = {wc:.4f};  {N} quỹ đạo × 60 T0 mỗi cấu hình\n")
    with Pool() as pool:
        for mode in (None, 'drive', 'drive+noise'):
            for kT in (0.1, 0.5):
                bs = [b for r in pool.map(_job, [(cfg, kT, 31000 + i, T, warm, mode) for i in range(N)]) for b in r]
                tot = [sum(b['pulses']) for b in bs]
                bud = [0.5 * (b['w_break']**2 - b['w_end']**2) * HW0 for b in bs]
                below = sum(1 for b in bs if b['w_end'] < wc) / len(bs)
                pul = [p for b in bs for p in b['pulses']]
                opened, eps_prop, eps_lo, eps_hi = [], [], [], []
                for b in bs:
                    acc = 0.0
                    for k, p in enumerate(b['pulses'], 1):
                        acc += p
                        if acc >= SN:
                            E = acc - SN
                            opened.append((k, E))
                            eps_lo.append(2 * E / (2 * k))          # n = 2k − 1
                            eps_hi.append(2 * E / (2 * k + 2))      # n = 2k + 1
                            eps_prop.append(2 * E / 4)              # n = 3 (đề xuất)
                            break
                name = {None: "không ngắt (tham chiếu)", 'drive': "ngắt LÁI cả đợt", 'drive+noise': "ngắt LÁI + NHIỄU cả đợt"}[mode]
                print(f"== {name}, k_BT = {kT*HW0:g} MeV: {len(bs)} đợt ==")
                print(f"   (A) tổng Carnot/đợt: TB {sum(tot)/len(tot):.2f}, trung vị {pct(tot,.5):.2f}, 90% {pct(tot,.9):.2f}, "
                      f"max {max(tot):.1f} MeV;  tỉ lệ đợt ≤ 25 MeV: {sum(1 for x in tot if x <= 25)/len(tot):.1%}")
                print(f"       ngân sách ½I(θ̇_gãy²−θ̇_cuối²): TB {sum(bud)/len(bud):.2f}, 10–90% {pct(bud,.1):.2f}–{pct(bud,.9):.2f}, "
                      f"max {max(bud):.1f} MeV;  tỉ lệ ngân sách ≤ 25 MeV: {sum(1 for x in bud if x <= 25)/len(bud):.1%};  "
                      f"tỉ lệ đợt tổng ≤ ngân sách: {sum(1 for x, y in zip(tot, bud) if x <= y)/len(bs):.1%}")
                print(f"   (B) θ̇ sau lần khóa kết thúc đợt < θ̇_gãy: {below:.1%}")
                print(f"   (C) xung: TB {sum(pul)/len(pul)/HW0:.3f} ħω0; P(≥1ħω0) = {sum(1 for p in pul if p >= HW0)/len(pul):.1%}, "
                      f"P(≥2ħω0) = {sum(1 for p in pul if p >= 2*HW0)/len(pul):.1%}")
                if opened:
                    ks = [o[0] for o in opened]
                    Es = [o[1] for o in opened]
                    print(f"       mở kênh {len(opened)}/{len(bs)}; k TB {sum(ks)/len(ks):.2f} (trung vị {pct(ks,.5)}); "
                          f"E_max TB {sum(Es)/len(Es):.2f} (trung vị {pct(Es,.5):.2f}) MeV; "
                          f"⟨ε⟩ n=2k+1…2k−1: {sum(eps_hi)/len(eps_hi):.2f}–{sum(eps_lo)/len(eps_lo):.2f} MeV; "
                          f"nếu cố định n=3: {sum(eps_prop)/len(eps_prop):.2f} MeV", flush=True)
            print()


if __name__ == "__main__":
    main()
