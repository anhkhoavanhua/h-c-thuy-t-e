"""Giải số chín câu hỏi mở của Part II (E_dual, HocThuyetE (7).pdf, mô hình V.C.1).

A1  Gãy liên kết trong phương trình chủ: thêm giếng ràng buộc V_b = −D e^{−y²/2w²}, y = x − κ sinθ, lực cực đại = R
    (w = 1, D = R w √e); MCWF — trung bình theo sự kiện của P_ngoài = P(|y| > 2w) và ⟨|ℓ|⟩ quanh lúc khởi phát đợt Carnot.
A2  Bloch–Redfield không secular (cùng bể phẳng T = 0 với Davies): ρ̇ = −i[H,ρ] − [û, Λρ] + [û, ρΛ†], Λ = (Γc/2) û₋;
    so trạng thái dừng, tính dương, độ ổn định Liouvillian với Davies.
A3  Giới hạn cổ điển: ħ_eff = s. Biến đổi chính tắc ⟹ trong mã (ħ = 1): I → I/s, ΔE → ΔE/s, g → g/√s, κ → κ/√s,
    bơm Γ0 → Γ0/s (giữ mô-men cổ điển), f theo ℓ thực = sℓ; năng lượng thực = s × mã. Đo năng lượng/bước, bước/đợt, công suất.
A4  Độ nhạy theo Γc, Γ0, ℓ0 (Liouvillian chính xác + dòng chuyển dời trong cơ sở riêng).
A5  Bơm Van der Pol lượng tử cho rô-to: L↑ = √γ e^{iθ}√(ℓ+½)₊, L↓ = √(γ/Ω0²) e^{−iθ} ℓ₊^{3/2} ⟹ d⟨ℓ⟩/dt = γ(ℓ+½ − ℓ³/Ω0²)
    (mô-men Rayleigh–Van der Pol của Part I, γ = 1/τ = 0.1, Ω0 = 1.05 θ̇_br).
A6  Davies ở nhiệt độ hữu hạn (cân bằng chi tiết: tốc độ lên = tốc độ xuống × e^{−ω/k_BT}), k_BT = 1 và 5 MeV.
A7  HEOM (Tanimura–Kubo) cho động cơ: bể Drude J(ω) = 2λγ_D ω/(ω² + γ_D²) ghép qua û, k_BT = 0.1, Matsubara + hiệu chỉnh
    đuôi Ishizaki–Tanimura; so với Davies từ cùng J(ω) ở ba cường độ λ; hội tụ theo (K, độ sâu).
A8  Darwinism với 10 bản ghi (S nhỏ ℓ ≤ 3, n ≤ 2): I(R:F) theo |F|, cặp bó sóng ràng buộc vs cặp trạng thái riêng.
A9  Phổ nucleon phát ra: đợt lượng tử (MCWF V.C.1) vs chuỗi cổ điển Part I vs chuẩn bay hơi Weisskopf (a = A/8).
Chạy: python3 giai_so_part2_mo.py A1 … A9
"""
import itertools
import math
import random
import sys
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

import giai_so_u_hermite_drude as UH
from giai_so_part2_master import liouvillian

T0 = 2 * math.pi
HW = 10.0


# ======================================================================= tiện ích chung
def rotor_osc(lmax, nmax):
    ls, cos_t, sin_t, pth, H, u = UH.basis(lmax=lmax, nmax=nmax)
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
    x = (a + a.conj().T) / math.sqrt(2)
    d = len(ls)
    up = np.diag(np.ones(d - 1), -1).astype(complex)
    return dict(ls=ls, d=d, cos=cos_t, sin=sin_t, H=H, u=u, a=a, x=x, up=up, nmax=nmax,
                Lop=np.kron(np.diag(ls.astype(float)), np.eye(nmax + 1)),
                Nn=np.kron(np.eye(d), a.conj().T @ a))


def pump(m, G0=0.2, l0=4.0, s=1.0):
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - (s * m["ls"].astype(float))**2 / l0**2)))
    return math.sqrt(G0 / s) * np.kron(m["up"] @ F, np.eye(m["nmax"] + 1))


def davies(H, u, Gc, kT=0.0):
    """Davies phẳng: xuống Γc|u_ij|², lên Γc|u_ij|² e^{−ω/kT} (cân bằng chi tiết)."""
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    Ls = [V @ (math.sqrt(Gc) * np.triu(ue, 1)) @ V.conj().T]
    if kT > 0:
        w = E[:, None] - E[None, :]                      # E_i − E_j
        B = math.sqrt(Gc) * np.tril(ue, -1) * np.sqrt(np.exp(-np.clip(w, 0, None) / kT))
        Ls.append(V @ B @ V.conj().T)
    return Ls, E, V


def steady(H, Ls):
    D = H.shape[0]
    L = liouvillian(H, Ls)
    b = np.zeros(D * D, complex); L[0, :] = np.eye(D).reshape(-1); b[0] = 1
    r = np.linalg.solve(L, b).reshape(D, D)
    return (r + r.conj().T) / 2


def labels(m, V):
    Ab = np.abs(m["Lop"])
    ll = np.array([float(np.real(V[:, j].conj() @ Ab @ V[:, j])) for j in range(V.shape[1])])
    nn = np.array([float(np.real(V[:, j].conj() @ m["Nn"] @ V[:, j])) for j in range(V.shape[1])])
    return ll, nn


def flux(m, rho, H, u, Gc):
    """Dòng chuyển dời Davies xuống ở trạng thái dừng (gần đúng secular): F_ij = Γc|u_ij|² p_j."""
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    R = Gc * np.abs(np.triu(ue, 1))**2
    p = np.real(np.diag(V.conj().T @ rho @ V))
    Fm = R * p[None, :]
    ll, nn = labels(m, V)
    cls = defaultdict(lambda: [0.0, 0.0])
    for i, j in zip(*np.nonzero(Fm > 1e-12)):
        dl, dn = round(ll[j]) - round(ll[i]), round(nn[j]) - round(nn[i])
        k = "ℓ→ℓ−1" if (dl == 1 and dn == 0) else "n→n−1" if (dl == 0 and dn == 1) else "khác"
        cls[k][0] += Fm[i, j]; cls[k][1] += Fm[i, j] * (E[j] - E[i])
    tot = sum(v[0] for v in cls.values()); et = sum(v[1] for v in cls.values())
    return {k: (v[0] / tot, v[1] / et) for k, v in cls.items()}, et / tot, tot * T0


def summary(m, rho):
    ev = lambda O: float(np.real(np.trace(O @ rho)))
    w = np.linalg.eigvalsh(rho)
    return ev(m["Lop"]), ev(m["H"]) if "Hb" not in m else ev(m["Hb"]), float(np.real(np.trace(rho @ rho))), w.min()


# ======================================================================= A1
def a1_model():
    lmax, nmax, w, R = 6, 12, 1.0, 3.0
    m = rotor_osc(lmax, nmax)
    d = m["d"]
    Y = np.kron(np.eye(d), m["x"]) - UH.KAP * np.kron(m["sin"], np.eye(nmax + 1))
    yv, Q = np.linalg.eigh(Y)
    D = R * w * math.sqrt(math.e)
    Vb = (Q * (-D * np.exp(-yv**2 / (2 * w * w)))) @ Q.conj().T
    Pout = (Q * (np.abs(yv) > 2 * w)) @ Q.conj().T
    H = m["H"] + Vb
    return m, H, Pout


def _a1_traj(args):
    seed, Tmax, dt = args
    m, H, Pout = a1_model()
    Ls = davies(H, m["u"], 0.5)[0] + [pump(m)]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    E0, V0 = np.linalg.eigh(H)
    rng = random.Random(seed)
    psi = V0[:, 0].astype(complex)                       # trạng thái cơ bản: hạt bị giữ trong giếng
    t, r, step = 0.0, rng.random(), 0
    ts, po, lab, jumps = [], [], [], []
    Labs = np.abs(m["Lop"])
    while t < Tmax:
        psi = U @ psi; t += dt; step += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            psi = Ls[k] @ psi; psi /= np.linalg.norm(psi)
            if k == 0:
                jumps.append(t)
            r = rng.random()
        if step % 10 == 0:
            p = psi / np.linalg.norm(psi)
            ts.append(t); po.append(np.real(np.vdot(p, Pout @ p))); lab.append(np.real(np.vdot(p, Labs @ p)))
    return np.array(ts), np.array(po), np.array(lab), np.array(jumps)


def part_a1(N=24, Tmax=600.0, dt=0.01):
    m, H, Pout = a1_model()
    E, V = np.linalg.eigh(H)
    pb = [float(np.real(V[:, j].conj() @ Pout @ V[:, j])) for j in range(8)]
    print("=== A1. Gãy liên kết trong phương trình chủ (giếng ràng buộc lực cực đại R = 3, w = 1; ℓ ≤ 6, n ≤ 12) ===")
    print("8 mức thấp của H + V_b: E = " + ", ".join(f"{e:.2f}" for e in E[:8]) + "; P_ngoài = " + ", ".join(f"{p:.3f}" for p in pb))
    with Pool() as pool:
        res = pool.map(_a1_traj, [(3000 + i, Tmax, dt) for i in range(N)])
    win = np.arange(-T0, T0 + 1e-9, 0.1)
    acc_p, acc_l, onl, n_on = np.zeros(len(win)), np.zeros(len(win)), [], 0
    allp, alll = [], []
    for ts, po, lab, jumps in res:
        keep = ts > Tmax / 5
        allp.append(po[keep]); alll.append(lab[keep])
        js = jumps[jumps > Tmax / 5 + T0]
        last = -1e9
        for tj in js:
            if tj - last >= T0 and tj < ts[-1] - T0:          # khởi phát đợt
                ip = np.interp(tj + win, ts, po); il = np.interp(tj + win, ts, lab)
                acc_p += ip; acc_l += il; n_on += 1; onl.append(np.interp(tj - 0.05, ts, lab))
            last = tj
    acc_p /= n_on; acc_l /= n_on
    gp, gl = np.concatenate(allp).mean(), np.concatenate(alll).mean()
    rate = sum(len(r_[3][r_[3] > Tmax / 5]) for r_ in res) / (N * 0.8 * Tmax / T0)
    print(f"{N} quỹ đạo × {Tmax:g}: tốc độ nhảy Carnot {rate:.3f}/T0, {n_on} lần khởi phát đợt")
    print(f"trung bình toàn thời gian: P_ngoài {gp:.4f}, ⟨|ℓ|⟩ {gl:.3f}")
    print(f"{'t − t_khởi phát (T0)':>21} {'P_ngoài':>9} {'/ TB':>6} {'⟨|ℓ|⟩':>7}")
    for tt in (-1.0, -0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5, 1.0):
        i = int(np.argmin(np.abs(win - tt * T0)))
        print(f"{tt:21.2f} {acc_p[i]:9.4f} {acc_p[i]/gp:6.2f} {acc_l[i]:7.3f}")
    print(f"⟨|ℓ|⟩ ngay trước khởi phát: TB {np.mean(onl):.3f}, trung vị {np.median(onl):.3f} "
          f"(ngưỡng cổ điển θ̇_br = 3.205 với I = 1)\n")


# ======================================================================= A2
def part_a2(lmax=5, nmax=4, Gc=0.5):
    m = rotor_osc(lmax, nmax)
    H, u = m["H"], m["u"]
    Lp = pump(m)
    Ld, E, V = davies(H, u, Gc)
    D = H.shape[0]
    I = np.eye(D)
    Lam = 0.5 * Gc * (V @ np.triu(V.conj().T @ u @ V, 1) @ V.conj().T)
    Lred = liouvillian(H, [Lp])
    uL, Ldu = u @ Lam, Lam.conj().T @ u
    Lred += -np.kron(uL, I) + np.kron(Lam, u.T) + np.kron(u, Lam.conj().T.T) - np.kron(I, Ldu.T)
    b = np.zeros(D * D, complex)
    rows = []
    for name, Lsup in (("Davies (secular, Lindblad)", liouvillian(H, Ld + [Lp])), ("Bloch–Redfield (không secular)", Lred)):
        lam = np.linalg.eigvals(Lsup)
        A_ = Lsup.copy(); A_[0, :] = I.reshape(-1); b[:] = 0; b[0] = 1
        rho = np.linalg.solve(A_, b).reshape(D, D); rho = (rho + rho.conj().T) / 2
        lm, eh, pu, mn = summary(m, rho)
        fr, eJ, rt = flux(m, rho, H, u, Gc)
        rows.append((name, lm, eh, pu, mn, np.max(lam.real), np.sort(-lam.real)[1], fr.get("ℓ→ℓ−1", (0, 0))[0], eJ, rho))
    E_, V_ = np.linalg.eigh(H)
    ue = V_.conj().T @ u @ V_
    gaps = np.abs(E_[:, None] - E_[None, :])
    off = (gaps > 1e-9) & (np.abs(ue) > 1e-6)
    ratio = Gc * np.abs(ue[off])**2 / gaps[off]
    print(f"=== A2. Bloch–Redfield vs Davies (ℓ ≤ {lmax}, n ≤ {nmax}, Γc = {Gc}, bơm V.C.1) ===")
    print(f"điều kiện secular Γc|u_ij|²/|ω_ij|: trung vị {np.median(ratio):.3f}, 90% {np.percentile(ratio,90):.2f}, "
          f"max {ratio.max():.1f}; tỉ lệ cặp > 0.1: {np.mean(ratio > 0.1):.1%}")
    print(f"{'phương trình':>32} {'⟨ℓ⟩':>7} {'⟨H⟩':>7} {'Tr ρ²':>7} {'min eig ρ':>10} {'max Re λ':>9} {'khe':>7} "
          f"{'ℓ→ℓ−1':>7} {'ΔE/bước':>8}")
    for r_ in rows:
        print(f"{r_[0]:>32} {r_[1]:7.3f} {r_[2]:7.3f} {r_[3]:7.4f} {r_[4]:10.1e} {r_[5]:9.1e} {r_[6]:7.4f} {r_[7]:7.1%} {r_[8]:8.3f}")
    print(f"khoảng cách vết giữa hai ρ dừng: {0.5*np.sum(np.abs(np.linalg.eigvalsh(rows[0][9]-rows[1][9]))):.4f}\n")


# ======================================================================= A3
def scaled_model(s, lmax, nmax):
    keep = (UH.I_, UH.DE, UH.G, UH.KAP)
    UH.I_, UH.DE, UH.G, UH.KAP = 1.0 / s, 1.0 / s, 0.1 / math.sqrt(s), (1 / 3) / math.sqrt(s)
    try:
        m = rotor_osc(lmax, nmax)
    finally:
        UH.I_, UH.DE, UH.G, UH.KAP = keep
    return m


def _a3_traj(args):
    s, lmax, nmax, seed, Tmax, dt = args
    m = scaled_model(s, lmax, nmax)
    H = m["H"]
    Ls = davies(H, m["u"], 0.5)[0] + [pump(m, s=s)]
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    rng = random.Random(seed)
    psi = np.zeros(H.shape[0], complex); psi[int(np.where(m["ls"] == 0)[0][0]) * (nmax + 1)] = 1
    t, r = 0.0, rng.random()
    jumps, lsum, lc = [], 0.0, 0
    step = 0
    while t < Tmax:
        psi = U @ psi; t += dt; step += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = [np.linalg.norm(L @ psi)**2 for L in Ls]
            k = 0 if rng.random() < ws[0] / sum(ws) else 1
            Eb = np.real(np.vdot(psi, H @ psi))
            psi = Ls[k] @ psi; psi /= np.linalg.norm(psi)
            if k == 0 and t > Tmax / 4:
                jumps.append((t, s * (Eb - np.real(np.vdot(psi, H @ psi)))))   # năng lượng THỰC
            r = rng.random()
        if t > Tmax / 4 and step % 50 == 0:
            p = psi / np.linalg.norm(psi)
            lsum += s * np.real(np.vdot(p, m["Lop"] @ p)); lc += 1
    return jumps, lsum / lc


def part_a3(N=12, Tmax=600.0, dt=0.01):
    print("=== A3. Giới hạn cổ điển ħ_eff = s → 0 (động lực học cổ điển giữ nguyên; năng lượng thực theo ħω0) ===")
    print(f"{'s':>6} {'cơ sở':>12} {'⟨ℓ⟩ thực = s⟨ℓ̃⟩':>16} {'ΔE/bước TB':>11} {'trung vị':>9} {'bước/đợt':>9} "
          f"{'ΔE/đợt':>8} {'nhảy/T0':>8} {'công suất (ħω0/T0)':>19}")
    cfg = {1.0: (8, 6), 0.5: (12, 7), 0.25: (20, 8), 0.15: (32, 9)}
    with Pool() as pool:
        for s, (lmax, nmax) in cfg.items():
            res = pool.map(_a3_traj, [(s, lmax, nmax, 4000 + i, Tmax, dt) for i in range(N)])
            dE = np.array([j[1] for rr in res for j in rr[0]])
            sizes, ens = [], []
            for rr in res:
                cur, ce, last = 0, 0.0, None
                for tj, e in rr[0]:
                    if last is None or tj - last < T0:
                        cur += 1; ce += e
                    else:
                        sizes.append(cur); ens.append(ce); cur, ce = 1, e
                    last = tj
            span = N * 0.75 * Tmax / T0
            print(f"{s:6.2f} {f'ℓ≤{lmax},n≤{nmax}':>12} {np.mean([rr[1] for rr in res]):16.3f} {dE.mean():11.3f} "
                  f"{np.median(dE):9.3f} {np.mean(sizes):9.2f} {np.mean(ens):8.3f} {len(dE)/span:8.3f} "
                  f"{dE.sum()/span:19.3f}", flush=True)
    print("(cổ điển Part I: 5.12 va chạm/đợt × 0.40–0.46 ħω0 = 2.0–2.4 ħω0/đợt)\n")


# ======================================================================= A4
def part_a4(lmax=5, nmax=4):
    m = rotor_osc(lmax, nmax)
    print(f"=== A4. Độ nhạy theo (Γc, Γ0, ℓ0) — Liouvillian chính xác, ℓ ≤ {lmax}, n ≤ {nmax} ===")
    print(f"{'Γc':>5} {'Γ0':>5} {'ℓ0':>4} {'⟨ℓ⟩':>7} {'Tr ρ²':>7} {'ℓ→ℓ−1 (bước)':>14} {'(năng lượng)':>13} {'n→n−1':>7} "
          f"{'ΔE/bước':>8} {'nhảy Carnot/T0':>15}")
    grid = [(0.5, 0.2, 4)] + [(g, 0.2, 4) for g in (0.25, 1.0, 2.0)] + [(0.5, g, 4) for g in (0.1, 0.4, 0.8)] + \
           [(0.5, 0.2, l) for l in (3, 5)] + [(0.25, 0.1, 3), (2.0, 0.8, 5)]
    for Gc, G0, l0 in grid:
        Ld = davies(m["H"], m["u"], Gc)[0]
        rho = steady(m["H"], Ld + [pump(m, G0, l0)])
        lm, eh, pu, mn = summary(m, rho)
        fr, eJ, rt = flux(m, rho, m["H"], m["u"], Gc)
        a, b_ = fr.get("ℓ→ℓ−1", (0, 0)); c_ = fr.get("n→n−1", (0, 0))[0]
        print(f"{Gc:5.2f} {G0:5.2f} {l0:4d} {lm:7.3f} {pu:7.4f} {a:14.1%} {b_:13.1%} {c_:7.1%} {eJ:8.3f} {rt:15.3f}", flush=True)
    print("(dòng tính theo gần đúng secular trên quần thể trạng thái riêng; hàng đầu = V.C.1)\n")


# ======================================================================= A5
def vdp_ops(m, gam=0.1, Om=1.05 * 3.205):
    lsf = m["ls"].astype(float)
    gup = np.sqrt(np.maximum(lsf, 0) + 0.5)
    gdn = np.where(lsf > 0, lsf**1.5, 0.0)
    Lup = math.sqrt(gam) * np.kron(m["up"] @ np.diag(gup), np.eye(m["nmax"] + 1))
    Ldn = math.sqrt(gam) / Om * np.kron(m["up"].conj().T @ np.diag(gdn), np.eye(m["nmax"] + 1))
    return [Lup, Ldn]


def part_a5(lmax=6, nmax=4):
    m = rotor_osc(lmax, nmax)
    print(f"=== A5. Bơm pha không kết hợp (V.C.1) vs bơm Van der Pol lượng tử (Rayleigh, γ = 0.1, Ω0 = 3.365), ℓ ≤ {lmax} ===")
    print(f"{'bơm':>34} {'Γc':>5} {'⟨ℓ⟩':>7} {'Tr ρ²':>7} {'ℓ→ℓ−1':>7} {'n→n−1':>7} {'ΔE/bước':>8} {'nhảy Carnot/T0':>15}")
    for name, P in (("bão hòa e^{iθ}f(ℓ̂), Γ0 = 0.2", [pump(m)]), ("Van der Pol lượng tử", vdp_ops(m))):
        for Gc in (0.0, 0.05, 0.5):
            Ld = davies(m["H"], m["u"], Gc)[0] if Gc > 0 else []
            rho = steady(m["H"], Ld + P)
            lm, eh, pu, mn = summary(m, rho)
            if Gc > 0:
                fr, eJ, rt = flux(m, rho, m["H"], m["u"], Gc)
                s1 = f"{fr.get('ℓ→ℓ−1',(0,0))[0]:7.1%} {fr.get('n→n−1',(0,0))[0]:7.1%} {eJ:8.3f} {rt:15.3f}"
            else:
                s1 = f"{'—':>7} {'—':>7} {'—':>8} {'—':>15}"
            print(f"{name:>34} {Gc:5.2f} {lm:7.3f} {pu:7.4f} {s1}", flush=True)
    print()


# ======================================================================= A6
def part_a6(lmax=5, nmax=4):
    m = rotor_osc(lmax, nmax)
    E, V = np.linalg.eigh(m["H"])
    print(f"=== A6. Davies ở nhiệt độ hữu hạn (cân bằng chi tiết), bơm V.C.1, ℓ ≤ {lmax}, n ≤ {nmax} ===")
    print(f"{'k_BT':>8} {'⟨ℓ⟩':>7} {'⟨H⟩':>7} {'⟨n⟩':>7} {'Tr ρ²':>7} {'nhảy xuống/T0':>14} {'nhảy lên/T0':>12} "
          f"{'p(+3)/p(−3) cặp':>16}")
    il3 = [j for j in range(len(E)) if abs(np.real(V[:, j].conj() @ m["Lop"] @ m["Lop"] @ V[:, j]) - 9) < 0.6][:2]
    for kT in (0.0, 0.1, 0.5):
        Ls = davies(m["H"], m["u"], 0.5, kT)[0]
        rho = steady(m["H"], Ls + [pump(m)])
        ev = lambda O: float(np.real(np.trace(O @ rho)))
        dn = ev(Ls[0].conj().T @ Ls[0]) * T0
        upr = ev(Ls[1].conj().T @ Ls[1]) * T0 if len(Ls) > 1 else 0.0
        p = np.real(np.diag(V.conj().T @ rho @ V))
        pr = p[il3[1]] / p[il3[0]] if len(il3) == 2 and p[il3[0]] > 0 else float('nan')
        print(f"{kT*10:6.1f}MeV {ev(m['Lop']):7.3f} {ev(m['H']):7.3f} {ev(m['Nn']):7.4f} {np.real(np.trace(rho@rho)):7.4f} "
              f"{dn:14.3f} {upr:12.4f} {pr:16.3f}", flush=True)
    print("(cặp ±3: hai trạng thái riêng gần suy biến có ⟨ℓ²⟩ ≈ 9; tỉ số quần thể của chúng)\n")


# ======================================================================= A7
def heom_steady(Hs, Q, Lp, lam, gD, kT, K, depth):
    """HEOM Drude–Lorentz (+ bơm Lindblad): dựng ma trận sinh đầy đủ, giải trực tiếp trạng thái dừng."""
    beta = 1 / kT
    nus = [gD] + [2 * math.pi * k / beta for k in range(1, K + 1)]
    cs = [lam * gD * (1 / math.tan(beta * gD / 2) - 1j)] + \
         [4 * lam * gD / beta * nk / (nk**2 - gD**2) for nk in nus[1:]]
    tail = 4 * lam * gD / beta * sum(1 / (nk**2 - gD**2) for nk in [2 * math.pi * k / beta for k in range(K + 1, 20000)])
    idx = [n for n in itertools.product(range(depth + 1), repeat=K + 1) if sum(n) <= depth]
    pos = {n: i for i, n in enumerate(idx)}
    NA, D = len(idx), Hs.shape[0]
    D2 = D * D
    I = np.eye(D)
    L0 = liouvillian(Hs, [Lp])
    QQ = Q @ Q
    dbl = np.kron(QQ, I) - 2 * np.kron(Q, Q.T) + np.kron(I, QQ.T)          # [Q,[Q,·]]
    comm = np.kron(Q, I) - np.kron(I, Q.T)                                  # [Q,·]
    G = np.zeros((NA * D2, NA * D2), complex)
    for i, n in enumerate(idx):
        sl = slice(i * D2, (i + 1) * D2)
        G[sl, sl] = L0 - sum(n[k] * nus[k] for k in range(K + 1)) * np.eye(D2) - tail * dbl
        for k in range(K + 1):
            n2 = list(n); n2[k] += 1
            j = pos.get(tuple(n2))
            if j is not None:
                G[sl, j * D2:(j + 1) * D2] = -1j * comm
            if n[k] > 0:
                n3 = list(n); n3[k] -= 1
                j = pos[tuple(n3)]
                G[sl, j * D2:(j + 1) * D2] = -1j * n[k] * (cs[k] * np.kron(Q, I) - np.conj(cs[k]) * np.kron(I, Q.T))
    b = np.zeros(NA * D2, complex)
    G[0, :] = 0; G[0, :D2] = I.reshape(-1); b[0] = 1
    x = np.linalg.solve(G, b)
    r0 = x[:D2].reshape(D, D)
    r0 = (r0 + r0.conj().T) / 2
    res = float(np.linalg.norm(G[1:] @ x - b[1:]))
    return r0 / np.trace(r0).real, res, NA


def part_a7(lmax=4, nmax=1, gD=2.0, kT=0.1):
    m = rotor_osc(lmax, nmax)
    H, u = m["H"], m["u"]
    Lp = pump(m)
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    print(f"=== A7. HEOM vs Davies cùng bể Drude (γ_D = {gD}, k_BT = {kT}), ghép qua û, + bơm V.C.1; ℓ ≤ {lmax}, n ≤ {nmax} ===")
    print(f"{'λ':>6} {'phương pháp':>26} {'⟨ℓ⟩':>7} {'⟨H⟩':>7} {'Tr ρ²':>7} {'min eig':>9} {'dư số':>11}")
    for lam in (0.05, 0.25, 1.0):
        J = lambda w: 2 * lam * gD * w / (w * w + gD * gD)
        w = E[:, None] - E[None, :]
        nb = lambda x: 1 / np.expm1(np.clip(x, 1e-12, None) / kT)
        Adn = np.where(w < 0, np.sqrt(np.abs(2 * J(-w) * (nb(-w) + 1))) * ue, 0)       # j → i với E_i < E_j
        Aup = np.where(w > 0, np.sqrt(np.abs(2 * J(w) * nb(w))) * ue, 0)
        Ld = [V @ Adn @ V.conj().T, V @ Aup @ V.conj().T]
        rD = steady(H, Ld + [Lp])
        ev = lambda O, r: float(np.real(np.trace(O @ r)))
        print(f"{lam:6.2f} {'Davies (Markov, secular)':>26} {ev(m['Lop'],rD):7.3f} {ev(H,rD):7.3f} "
              f"{np.real(np.trace(rD@rD)):7.4f} {np.linalg.eigvalsh(rD).min():9.1e} {'—':>11}", flush=True)
        for K, depth in ((3, 2), (6, 2), (2, 3), (3, 3)):
            rH, drift, NA = heom_steady(H, u, Lp, lam, gD, kT, K, depth)   # drift = dư số của nghiệm
            print(f"{lam:6.2f} {f'HEOM K={K}, độ sâu {depth} ({NA} ADO)':>26} {ev(m['Lop'],rH):7.3f} {ev(H,rH):7.3f} "
                  f"{np.real(np.trace(rH@rH)):7.4f} {np.linalg.eigvalsh(rH).min():9.1e} {drift:11.1e}", flush=True)
    print("(HEOM: phân rã Drude + K số hạng Matsubara + hiệu chỉnh đuôi Ishizaki–Tanimura; trạng thái dừng giải trực tiếp)\n")


# ======================================================================= A8
def part_a8():
    import giai_so_wdw_pw_e7 as W
    W.LMAX, W.NMAX = 3, 2
    NT = 10
    ls, H, Ls = W.channels()
    E, V = np.linalg.eigh(H)
    n = np.arange(W.NMAX + 1)

    def packet(th0, l0):
        cr = np.exp(-(ls - l0)**2 / (4 * 0.7**2) - 1j * ls * th0); cr /= np.linalg.norm(cr)
        al = ((1 / 3) * math.sin(th0) + 1j * (1 / 3) * l0 * math.cos(th0)) / math.sqrt(2)
        cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in n]) * math.exp(-abs(al)**2 / 2)
        return np.kron(cr, cx / np.linalg.norm(cx))
    Lop = np.kron(np.diag(ls.astype(float)), np.eye(W.NMAX + 1))
    il = [j for j in range(len(E)) if abs(np.real(V[:, j].conj() @ Lop @ Lop @ V[:, j]) - 4) < 0.8]
    pairs = {"C: bó sóng ràng buộc θ0 = 0 / π (ℓ0 = 2)": (packet(0, 2), packet(math.pi, 2)),
             "E: cơ bản / mức |ℓ| = 2": (V[:, 0], V[:, il[0]])}

    def ent(M):
        w = np.linalg.eigvalsh(M); w = w[w > 1e-13]
        return float(-np.sum(w * np.log2(w)))
    print(f"=== A8. Darwinism với {NT} bản ghi (S {len(ls)*(W.NMAX+1)} chiều) ===")
    rng = random.Random(1)
    for tdiv in (5, 10):
        W.TAU = T0 / tdiv
        U, d = W.collision_unitary(H, Ls)
        print(f"-- τ = T0/{tdiv} ({NT} bản ghi = {NT/tdiv:.1f} T0)")
        for name, (A_, B_) in pairs.items():
            psi = np.zeros((d,) + (3,) * NT + (2,), complex)
            psi[(slice(None),) + (0,) * NT + (0,)] = A_ / math.sqrt(2)
            psi[(slice(None),) + (0,) * NT + (1,)] = B_ / math.sqrt(2)
            for k in range(NT):
                psi = W.apply_tick(psi, U, d, k)

            def red(keep):
                other = [ax for ax in range(NT + 2) if ax not in keep]
                Mx = np.transpose(psi, keep + other).reshape(int(np.prod([psi.shape[k] for k in keep])), -1)
                return ent(Mx @ Mx.conj().T) if Mx.shape[0] <= Mx.shape[1] else ent(Mx.conj().T @ Mx)
            SR = red([NT + 1])
            Iv = []
            for m_ in range(0, NT + 1):
                frs = list(itertools.combinations(range(1, NT + 1), m_)); rng.shuffle(frs)
                vals = []
                for F in frs[:3]:
                    F = list(F)
                    vals.append(SR + (red(F) if F else 0.0) - red([NT + 1] + F))
                Iv.append(np.mean(vals))
            md = next((k for k, v in enumerate(Iv) if v >= 0.9), None)
            print(f"   {name}: I(R:F) (bit) |F| = 0…{NT}: " + " ".join(f"{v:.2f}" for v in Iv) +
                  f";  R_0.1 = {NT/md if md else 0:.2f}", flush=True)
    print()


# ======================================================================= A9
def part_a9(N=24, Tmax=1500.0, dt=0.01):
    from giai_so_part2_master import _traj
    with Pool() as pool:
        res = pool.map(_traj, [(5000 + i, Tmax, dt, False) for i in range(N)])
    first, burstE, ks = [], [], []
    for rr in res:
        cj = [(j[0], j[2]) for j in rr[0] if j[1] == "Carnot"]
        cur, last = [], None
        for tj, e in cj:
            if last is None or tj - last < T0:
                cur.append(e)
            else:
                first.append(cur[0]); burstE.append(sum(cur)); ks.append(len(cur)); cur = [e]
            last = tj
    first, burstE, ks = np.array(first) * HW, np.array(burstE) * HW, np.array(ks)
    Sn = 8.0
    print(f"=== A9. Phổ nucleon phát ra: đợt lượng tử V.C.1 ({len(first)} đợt) vs Part I vs bay hơi Weisskopf ===")
    print(f"năng lượng bước đầu của đợt: TB {first.mean():.1f} MeV, trung vị {np.median(first):.1f}; vượt S_n = 8 MeV ngay "
          f"bước đầu: {np.mean(first > Sn):.1%}")
    print(f"năng lượng cả đợt: TB {burstE.mean():.1f} MeV; đợt vượt S_n: {np.mean(burstE > Sn):.1%}")
    # exciton: sau k bước n = 2k + 1, ⟨ε⟩ = 2(E_exc − S_n)/(n + 1) (Part I, Eq. 22), tính tại bước đầu vượt S_n
    eps = []
    for rr in res:
        cj = [(j[0], j[2] * HW) for j in rr[0] if j[1] == "Carnot"]
        cur, last = [], None
        for tj, e in cj + [(1e18, 0)]:
            if last is not None and tj - last >= T0:
                acc = np.cumsum(cur)
                hit = np.nonzero(acc > Sn)[0]
                if len(hit):
                    k = hit[0] + 1; n_ = 2 * k + 1
                    eps.append(2 * (acc[hit[0]] - Sn) / (n_ + 1))
                cur = []
            cur.append(e); last = tj
    eps = np.array(eps)
    a = 70 / 8.0
    print(f"{'mô hình':>40} {'bước tới S_n':>13} {'⟨ε⟩ nucleon (MeV)':>18}")
    print(f"{'lượng tử (đợt V.C.1, exciton n = 2k+1)':>40} {'1 (trung vị)':>13} {eps.mean():18.2f}")
    print(f"{'cổ điển Part I (Eq. 22)':>40} {'3 (trung vị)':>13} {'0.9–1.6':>18}")
    for Ex in (12.0, 20.0, burstE.mean()):
        T = math.sqrt(Ex / a)
        print(f"{f'bay hơi Weisskopf E* = {Ex:.1f} MeV (T = {T:.2f})':>40} {'—':>13} {2*T:18.2f}")
    print("(Weisskopf: phổ ε e^{−ε/T}, ⟨ε⟩ = 2T, T = √(E*/a), a = A/8 MeV⁻¹ — chuẩn thống kê, không phải dữ liệu đo)\n")


if __name__ == "__main__":
    parts = {f"A{i}": f for i, f in enumerate([None, part_a1, part_a2, part_a3, part_a4, part_a5, part_a6, part_a7,
                                                part_a8, part_a9]) if f}
    for name in sys.argv[1:] or list(parts):
        parts[name.upper()]()
