"""Hai câu hỏi mở cuối của Part II (E_dual, HocThuyetE (7).pdf).

H1  Sàn χ² đạt được: χ² giữa HAI phép đo độc lập cùng phản ứng Cu(n,xn) 14.1 MeV (Baba 1988 vs Takahashi 1988), cùng
    cách xử lý (tích phân góc, sai số từng bin, chuẩn hóa hình dạng). Nếu hai thí nghiệm lệch nhau quá sai số khai báo
    thì χ²/k ≈ 1 không đạt được với bất kỳ mô hình nào; dùng độ lệch giữa hai phép đo làm sai số hệ thống thực tế.
H2  σ_inv của neutron từ thế quang học toàn cục Koning–Delaroche 2003 (Nucl. Phys. A 713, 231): V_V, W_V, W_D, V_so,
    W_so; Numerov theo (ℓ, j), khớp Riccati–Bessel ⟹ S_ℓj ⟹ σ_R(ε) = π/k² Σ (2j+1)/2 (1 − |S_ℓj|²). Bảng σ_R(A, ε)
    lưu exfor_data/om_sigma_n.npz; giai_so_part2_mo5 dùng bảng này khi đặt biến môi trường E_DUAL_OM=1.
S   Giới hạn cổ điển của số va chạm/đợt: mô hình chặn gãy (giếng ràng buộc lực cực đại R) dưới phép thu nhỏ ħ_eff = s:
    I → I/s, ΔE → ΔE/s, g → g/√s, κ → κ/√s, giếng D → D/s, w → w/√s, bơm Γ0 → Γ0/s (f theo ℓ thực = sℓ);
    năng lượng thực = s × mã. Đo va chạm/đợt, ΔE/va chạm, lưu, CV theo s; so 5.12 của Part I.
Chạy: python3 giai_so_part2_mo6.py [H1|H2|S]
"""
import math
import os
import random
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(HERE, "exfor_data")
T0 = 2 * math.pi


# ======================================================================= H1
def part_h1():
    import giai_so_part2_mo5 as M
    sets = {d[0]: d for d in M.datasets()}
    b = sets["Cu(n,xn) 14.1 — Baba 1988"]; t = sets["Cu(n,xn) 14.1 — Takahashi 1988"]
    lo, hi = 1.0, 10.1
    eb, sb, eb_err = b[5], b[6], b[7]
    et, st, st_err = t[5], t[6], t[7]
    mb = (eb >= lo) & (eb <= hi) & np.isfinite(sb)
    e = eb[mb]
    Sb, Eb = sb[mb], eb_err[mb]
    St = np.interp(e, et, st); Et = np.interp(e, et, st_err)
    w = np.gradient(e)
    Ab, At = np.sum(Sb * w), np.sum(St * w)
    nb, nt = Sb / Ab, St / At
    nbe, nte = Eb / Ab, Et / At
    k = len(e) - 1
    print("=== H1. Sàn χ²: hai phép đo độc lập Cu(n,xn) 14.1 MeV (Baba 1988 vs Takahashi 1988), cửa sổ 1–10.1 MeV ===")
    print(f"tích phân: Baba {Ab:.0f} mb, Takahashi {At:.0f} mb (lệch {100*(At/Ab-1):+.1f}%); ⟨ε⟩ {np.sum(e*nb*w):.2f} vs "
          f"{np.sum(e*nt*w):.2f} MeV; sai số thống kê trung vị/bin: Baba {np.median(nbe/nb):.1%}, Takahashi {np.median(nte/nt):.1%}")
    for extra in (0.0, 0.02, 0.05, 0.08):
        chi = np.sum((nb - nt)**2 / (nbe**2 + nte**2 + (extra * nb)**2 + (extra * nt)**2))
        print(f"   thêm {extra:.0%} không tương quan mỗi phép đo: χ²/k = {chi/k:.2f}  (k = {k})")
    rel = np.abs(nb - nt) / (0.5 * (nb + nt))
    print(f"độ lệch tương đối giữa hai phép đo theo bin: trung vị {np.median(rel):.1%}, 90% {np.percentile(rel, 90):.1%}\n")
    return float(np.median(rel))


# ======================================================================= H2
def kd_neutron(A, Z, E):
    """Thế quang học toàn cục KD 2003 cho neutron: trả (độ sâu, bán kính, độ nhòe) các thành phần."""
    N = A - Z
    Ef = -11.2814 + 0.02646 * A
    f = E - Ef
    v1 = 59.30 - 21.0 * (N - Z) / A - 0.024 * A
    v2 = 0.007228 - 1.48e-6 * A; v3 = 1.994e-5 - 2.0e-8 * A; v4 = 7e-9
    w1 = 12.195 + 0.0167 * A; w2 = 73.55 + 0.0795 * A
    d1 = 16.0 - 16.0 * (N - Z) / A; d2 = 0.0180 + 0.003802 / (1 + math.exp((A - 156.0) / 8.0)); d3 = 11.5
    vso1 = 5.922 + 0.0030 * A; vso2 = 0.0040; wso1 = -3.1; wso2 = 160.0
    VV = v1 * (1 - v2 * f + v3 * f**2 - v4 * f**3)
    WV = w1 * f**2 / (f**2 + w2**2)
    WD = d1 * f**2 * math.exp(-d2 * f) / (f**2 + d3**2)
    Vso = vso1 * math.exp(-vso2 * f); Wso = wso1 * f**2 / (f**2 + wso2**2)
    A3 = A**(1 / 3)
    rV, aV = 1.3039 - 0.4054 / A3, 0.6778 - 1.487e-4 * A
    rD, aD = 1.3424 - 0.01585 * A3, 0.5446 - 1.656e-4 * A
    rso, aso = 1.1854 - 0.647 / A3, 0.59
    return dict(VV=VV, WV=WV, WD=WD, Vso=Vso, Wso=Wso, RV=rV * A3, aV=aV, RD=rD * A3, aD=aD, Rso=rso * A3, aso=aso)


def riccati(l, x):
    """Trả (F_l, G_l): F = x j_l(x) ~ sin(x − lπ/2) (truy hồi lùi Miller), G = −x y_l(x) ~ cos(x − lπ/2) (truy hồi tiến)."""
    n0 = l + int(x) + 30
    jn1, jn = 0.0, 1e-30
    vals = {}
    for n in range(n0, -1, -1):                       # j_{n−1} = (2n+1)/x j_n − j_{n+1}
        jm = (2 * n + 1) / x * jn - jn1
        jn1, jn = jn, jm
        vals[n - 1] = jm
        if abs(jm) > 1e250:
            jn1 /= 1e250; jn /= 1e250; vals = {k: v / 1e250 for k, v in vals.items()}
    scale = (math.sin(x) / x) / vals[0]
    F = x * vals[l] * scale
    y0, y1 = -math.cos(x) / x, -math.cos(x) / x**2 - math.sin(x) / x
    yl = y0 if l == 0 else y1
    for n in range(1, l):
        y0, y1 = y1, (2 * n + 1) / x * y1 - y0
        yl = y1
    return F, -x * yl


def sigma_R(A, Z, E, Rmax=15.0, h=0.05):
    mu = 939.565 * A / (A + 1)
    k = math.sqrt(2 * mu * E) / 197.327
    p = kd_neutron(A, Z, E)
    r = np.arange(h, Rmax + 1e-9, h)
    fV = 1 / (1 + np.exp((r - p["RV"]) / p["aV"]))
    xD = np.exp((r - p["RD"]) / p["aD"]); gD = 4 * xD / (1 + xD)**2          # −4a_D df_D/dr (dạng bề mặt)
    xs = np.exp((r - p["Rso"]) / p["aso"]); dfso = -xs / (1 + xs)**2 / p["aso"]
    Uc = -p["VV"] * fV - 1j * p["WV"] * fV - 1j * p["WD"] * gD
    Uso = 2.0 * (p["Vso"] + 1j * p["Wso"]) * dfso / r                          # λ_π² = 2 fm²; nhân l·σ
    fac = 2 * mu / 197.327**2
    c = h * h / 12
    sig = 0.0
    for l in range(int(k * Rmax) + 4):
        for j in ((0.5,) if l == 0 else (l - 0.5, l + 0.5)):
            lsig = j * (j + 1) - l * (l + 1) - 0.75                             # l·σ = 2 l·s
            q = fac * (E - Uc - Uso * lsig) - l * (l + 1) / r**2
            u = np.zeros(len(r), complex)
            u[0], u[1] = r[0]**(l + 1), r[1]**(l + 1)
            for i in range(1, len(r) - 1):
                u[i + 1] = (2 * u[i] * (1 - 5 * c * q[i]) - u[i - 1] * (1 + c * q[i - 1])) / (1 + c * q[i + 1])
                if abs(u[i + 1]) > 1e150:
                    u[: i + 2] /= 1e150
            F1, G1 = riccati(l, k * r[-3]); F2, G2 = riccati(l, k * r[-1])
            Hp1, Hm1, Hp2, Hm2 = G1 + 1j * F1, G1 - 1j * F1, G2 + 1j * F2, G2 - 1j * F2
            R_ = u[-3] / u[-1]
            S = (Hm1 - R_ * Hm2) / (Hp1 - R_ * Hp2)
            sig += (2 * j + 1) / 2 * (1 - abs(S)**2)
    return math.pi / k**2 * sig


def part_h2():
    print("=== H2. σ_R(ε) neutron từ thế quang học KD 2003 (Numerov ℓ, j) so với Dostrovsky ===")
    import giai_so_part2_mo5 as M
    As = np.arange(50, 81, 2)
    eps = np.concatenate([np.linspace(0.05, 2, 12), np.linspace(2.5, 25, 20)])
    tab = np.zeros((len(As), len(eps)))
    with Pool() as pool:
        for ia, A in enumerate(As):
            Z = int(round(A / (1.98 + 0.0155 * A**(2 / 3))))
            tab[ia] = pool.starmap(sigma_R, [(int(A), Z, float(e)) for e in eps])
    np.savez(os.path.join(DIR, "om_sigma_n.npz"), A=As, eps=eps, sigma=tab)
    print(f"{'A':>4} {'ε (MeV)':>8} {'σ_R KD (mb)':>12} {'Dostrovsky (mb)':>16} {'tỉ số':>6}")
    for A in (58, 64, 74):
        for e in (0.5, 1, 2, 5, 10, 14):
            s_om = np.interp(e, eps, tab[list(As).index(A)]) * 10
            s_d = float(M.sigma_inv("n", np.array([e]), 0, A)[0]) * 10
            print(f"{A:4d} {e:8.1f} {s_om:12.0f} {s_d:16.0f} {s_om/s_d:6.2f}")
    print("kiểm: σ_non đo Cu(n) 14.1 MeV = 1486 mb; KD cho A = 64:",
          f"{np.interp(14.1, eps, tab[list(As).index(64)])*10:.0f} mb\n")


# ======================================================================= S
def scaled_bond_model(s, lmax, nmax, R=3.0, w=1.0, pstar=0.3, Gc=0.5, G0=2.0, l0=4.0):
    import giai_so_u_hermite_drude as UH
    import giai_so_part2_mo as P
    keep = (UH.I_, UH.DE, UH.G, UH.KAP)
    UH.I_, UH.DE, UH.G, UH.KAP = 1.0 / s, 1.0 / s, 0.1 / math.sqrt(s), (1 / 3) / math.sqrt(s)
    try:
        m = P.rotor_osc(lmax, nmax)
        kap = UH.KAP
    finally:
        UH.I_, UH.DE, UH.G, UH.KAP = keep
    d = m["d"]
    Y = np.kron(np.eye(d), m["x"]) - kap * np.kron(m["sin"], np.eye(nmax + 1))
    yv, Q = np.linalg.eigh(Y)
    ws = w / math.sqrt(s); D = R * w * math.sqrt(math.e) / s
    H = m["H"] + (Q * (-D * np.exp(-yv**2 / (2 * ws * ws)))) @ Q.conj().T
    Pout = (Q * (np.abs(yv) > 2 * ws)) @ Q.conj().T
    E, V = np.linalg.eigh(H)
    po = np.real(np.einsum("ij,ik,kj->j", V.conj(), Pout, V))
    ue = V.conj().T @ m["u"] @ V
    Lc = V @ (math.sqrt(Gc) * np.triu(ue, 1) * (po > pstar)[None, :]) @ V.conj().T
    Lp = P.pump(m, G0, l0, s=s)
    return m, H, [Lc, Lp], V


def _s_traj(args):
    s, lmax, nmax, seed, Tmax, dt = args
    m, H, Ls, V = scaled_bond_model(s, lmax, nmax)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = V[:, 0].astype(complex)
    t, r, jumps = 0.0, rng.random(), []
    while t < Tmax:
        psi = U @ psi; t += dt
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            if k == 0:
                Eb = np.real(np.vdot(psi, H @ psi)); lb = np.real(np.vdot(psi, Labs @ psi))
            psi = Ls[k] @ psi; psi /= np.linalg.norm(psi)
            if k == 0:
                jumps.append((t, s * (Eb - np.real(np.vdot(psi, H @ psi))), s * lb))   # đơn vị thực
            r = rng.random()
    return jumps


def part_s(N=16):
    from giai_so_part2_mo3 import burst_stats
    print("=== S. Giới hạn cổ điển của mô hình chặn gãy (ℓ0 = 4, Γ0 = 2, Γc = 0.5, p* = 0.3) theo ħ_eff = s ===")
    print(f"{'s':>5} {'cơ sở':>12} {'đợt':>5} {'va chạm/đợt':>12} {'ΔE/va chạm (ħω0)':>17} {'ΔE/đợt':>8} {'lưu (T0)':>9} "
          f"{'CV':>5} {'ℓ thực khởi phát':>16}")
    cfg = ((1.0, 6, 12, 2400.0, 0.01), (0.5, 11, 20, 1600.0, 0.01), (0.25, 20, 34, 1200.0, 0.02))
    with Pool() as pool:
        for s, lmax, nmax, Tmax, dt in cfg:
            res = pool.map(_s_traj, [(s, lmax, nmax, 14000 + i, Tmax, dt) for i in range(N)])
            st = burst_stats(res, Tmax / 6)
            eb = []
            for jl in res:
                cur, last = 0.0, None
                for tj, e, lb in [j for j in jl if j[0] > Tmax / 6]:
                    if last is not None and tj - last >= T0:
                        eb.append(cur); cur = 0.0
                    cur += e; last = tj
            print(f"{s:5.2f} {f'ℓ≤{lmax},n≤{nmax}':>12} {st['n']:5d} {st['size']:12.2f} {st['eJ']:17.3f} "
                  f"{np.mean(eb) if eb else float('nan'):8.3f} {st['dwell']:9.2f} {st['cv']:5.2f} {st['lon']:16.2f}", flush=True)
    print("(Part I cổ điển: 5.12 va chạm/đợt × 0.40–0.46 ħω0 ≈ 2.2 ħω0/đợt; ngưỡng p_θ = 3.56)\n")


if __name__ == "__main__":
    for a in sys.argv[1:] or ["H1", "H2", "S"]:
        {"H1": part_h1, "H2": part_h2, "S": part_s}[a.upper()]()
