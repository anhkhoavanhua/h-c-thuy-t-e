"""Phát neutron THẬT trong động cơ E_dual (Part I, Langevin đầy đủ) — tương lai sau tán xạ có quay về Quá khứ nhanh hơn?

Thêm bể kích thích nội E_res (MeV) nhận mọi năng lượng Carnot. Bể nung nóng hạt nhân (khí Fermi):
  T_nuc = √(E_res / a), a = A/8 MeV⁻¹ (A = 70) ⟹ nhiễu Langevin dùng k_BT = max(1 MeV, T_nuc).
Phát thật: khi E_res ≥ S_n = 8 MeV, một neutron thoát với ε = (E_res − S_n)·Beta(2, n−1), n = 2k+1
(k = số va chạm kể từ lần phát trước); E_res ← E_res − S_n − ε.

  (a) chuẩn Part I: k_BT = 1 MeV cố định (không phản hồi)
  (b) phản hồi nhiệt, KHÔNG phát: E_res chỉ tăng ⟹ hạt nhân nóng dần
  (c) phản hồi nhiệt + phát thật
Đo: lần gãy / T0, tỉ lệ thời gian ở Quá khứ, ⟨N⟩ va chạm/đợt, nhịp phát, phổ ε, T_nuc, và thời gian từ lúc phát
tới khi tái khóa giữ được (quay về Quá khứ), so với thời gian tương ứng từ lúc E_exc vượt S_n ở (a).
Chạy: python3 giai_so_phat_that.py   (≈ 5–10 phút, 8 lõi)
"""
import math
import random
from multiprocessing import Pool

import numpy as np

import giai_so_dong_co_langevin as L

T0 = L.T0
import os
A_MASS, SN = 70.0, 8.0
KT0 = float(os.environ.get('KT0', '1.0'))  # MeV — nhiệt độ nền của bể
ALEV = A_MASS / 8.0                         # mật độ mức a (MeV⁻¹)
HW = 10.0


def run(mode, seed, T, f=1.05, gam=0.1, gx=0.70, dt=2e-3):
    rng = random.Random(seed); nrng = np.random.default_rng(seed)
    Om = f * L.WC
    s = [0.0, 0.82 * L.WC, 0.0, L.BASE['kappa'] * 0.82 * L.WC]
    t, locked, ref = 0.0, True, 0.0
    Eres, k_since = 0.0, 0
    ev, emits, Tlog = [], [], []
    sc = lambda q, h: L.rk4_step(lambda z: L.f_con(z, gam, Om, gx), q, h)
    sf = lambda q, h: L.rk4_step(lambda z: L.f_free(z, gam, Om, gx), q, h)
    lam_hit = lambda q: abs(L.acc_con(q, gam, Om, gx)[2]) >= L.R
    nxt = 0.0
    while t < T:
        kT_MeV = KT0 if mode == 'a' else max(KT0, math.sqrt(max(Eres, 0.0) / ALEV))
        kT = kT_MeV / HW
        sig, sigx = math.sqrt(2 * gam * kT), math.sqrt(2 * gx * kT)
        if t >= nxt:
            Tlog.append(kT_MeV); nxt += 0.5
        h = dt
        if locked:
            if lam_hit(s):
                ev.append(('gãy', t)); locked, ref = False, 0.0; continue
            s1 = sc(s, dt)
            if lam_hit(s1):
                h, s1 = L.locate(sc, s, dt, lam_hit)
        else:
            s1 = sf(s, dt); f1 = L.fval(s1)
            if ref == 0.0:
                ref = math.copysign(1.0, f1)
            elif f1 * ref <= 0:
                h, s1 = L.locate(sf, s, dt, lambda q: L.fval(q) * ref <= 0)
                s1, Lam, u = L.impulse(s1, L.BASE)
                c = math.cos(s1[0])
                meff = L.BASE['m'] * L.BASE['I'] / (L.BASE['I'] + L.BASE['m'] * L.BASE['kappa']**2 * c * c)
                dE = 0.5 * meff * u * u * HW
                ev.append(('khóa', t + h, dE)); locked = True
                if mode != 'a':
                    Eres += dE; k_since += 1
                    if mode == 'c' and Eres >= SN:
                        n = max(2, 2 * k_since + 1)
                        eps = (Eres - SN) * nrng.beta(2, n - 1)
                        emits.append((t + h, eps, k_since)); Eres -= SN + eps; k_since = 0
        s, t = list(s1), t + h
        kick = sig * math.sqrt(h) * rng.gauss(0, 1); kx = sigx * math.sqrt(h) * rng.gauss(0, 1)
        if locked:
            c = math.cos(s[0])
            s[1] += (kick + L.BASE['kappa'] * c * kx) / (L.BASE['I'] + L.BASE['m'] * L.BASE['kappa']**2 * c * c)
            s[3] = L.BASE['kappa'] * c * s[1]
        else:
            s[1] += kick / L.BASE['I']; s[3] += kx / L.BASE['m']
    return dict(mode=mode, ev=ev, emits=emits, T=T, Tlog=Tlog, Eres=Eres)


def relock_after(ev, times):
    """Thời gian từ mỗi mốc tới lần khóa đầu tiên giữ được ≥ T0 (không gãy lại trong T0) — quay về Quá khứ."""
    out = []
    for t0 in times:
        later = [e for e in ev if e[1] >= t0]
        for i, e in enumerate(later):
            if e[0] == 'khóa':
                nb = next((x for x in later[i + 1:] if x[0] == 'gãy'), None)
                if nb is None or nb[1] - e[1] >= T0:
                    out.append(e[1] - t0); break
    return np.array(out)


def summarize(rs, label, T_each):
    tot = len(rs) * T_each / T0
    nb = sum(sum(1 for e in r['ev'] if e[0] == 'gãy') for r in rs)
    nk = sum(sum(1 for e in r['ev'] if e[0] == 'khóa') for r in rs)
    # tỉ lệ thời gian khóa liên tục ≥ T0 (Quá khứ) — ước theo các khoảng giữa 'khóa' và 'gãy' kế tiếp
    past = 0.0
    for r in rs:
        e = r['ev']
        for a, b in zip(e, e[1:]):
            if a[0] == 'khóa' and b[0] == 'gãy' and b[1] - a[1] >= T0:
                past += b[1] - a[1] - T0
    emits = [x for r in rs for x in r['emits']]
    Tm = np.mean([np.mean(r['Tlog'][len(r['Tlog']) // 5:]) for r in rs])
    line = (f"  {label:34s} gãy/T0 = {nb / tot:.3f}, va chạm/T0 = {nk / tot:.3f}, Quá khứ = {100 * past / (len(rs) * T_each):4.1f}%,"
            f" T_nuc TB = {Tm:.2f} MeV")
    if emits:
        eps = np.array([x[1] for x in emits])
        line += f", phát/T0 = {len(emits) / tot:.3f}, ⟨ε⟩ = {eps.mean():.2f} MeV (trung vị {np.median(eps):.2f})"
    if label.startswith('(b)'):
        line += f", E_res cuối TB = {np.mean([r['Eres'] for r in rs]):.0f} MeV"
    print(line)


def main(N=8, T=100 * T0):
    with Pool(8) as p:
        res = p.starmap(run, [(m, 500 + s, T) for m in 'abc' for s in range(N)])
    print(f"Phát neutron thật: k_BT nền = {KT0} MeV, f = 1.05, A = {A_MASS:.0f}, a = {ALEV:.2f} MeV⁻¹, S_n = {SN} MeV; {N} quỹ đạo × {T / T0:.0f} T0 mỗi biến thể")
    for m, lab in (('a', f'(a) chuẩn, k_BT = {KT0} MeV cố định'), ('b', '(b) nung nóng, KHÔNG phát'), ('c', '(c) nung nóng + phát thật')):
        summarize([r for r in res if r['mode'] == m], lab, T)
    # thời gian quay về Quá khứ sau tán xạ
    ta = []
    for r in [r for r in res if r['mode'] == 'a']:
        acc, marks, in_burst = 0.0, [], False
        e = r['ev']
        for i, x in enumerate(e):
            if x[0] == 'khóa':
                acc += x[2]
                if acc >= SN:
                    marks.append(x[1]); acc = 0.0
        ta.append(relock_after(e, marks))
    tc = [relock_after(r['ev'], [x[0] for x in r['emits']]) for r in res if r['mode'] == 'c']
    ta = np.concatenate(ta); tc = np.concatenate(tc)
    print(f"\n  thời gian từ mốc tán xạ tới khi quay về Quá khứ (khóa giữ ≥ T0):")
    print(f"    (a) mốc = E_exc tích lũy vượt S_n (không phát):  TB {ta.mean() / T0:.2f} T0, trung vị {np.median(ta) / T0:.2f} T0 (n = {ta.size})")
    print(f"    (c) mốc = neutron thoát thật:                    TB {tc.mean() / T0:.2f} T0, trung vị {np.median(tc) / T0:.2f} T0 (n = {tc.size})")
    print(f"    tỉ số (c)/(a) = {tc.mean() / ta.mean():.2f}")


if __name__ == '__main__':
    main()
