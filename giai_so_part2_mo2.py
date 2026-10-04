"""Ba câu hỏi còn mở sau A1–A9 của Part II (E_dual, HocThuyetE (7).pdf).

Q1  Gãy kiểu Part I trong phương trình chủ: giếng ràng buộc lực cực đại R (như A1) + kênh Carnot CHỈ tác dụng lên trạng
    thái riêng không liên kết (P_ngoài > p*, hạt đã gãy khỏi ràng buộc) — Part I: tiêu tán chỉ xảy ra ở va chạm tái khóa,
    pha khóa bảo toàn. Đo: ngưỡng ⟨|ℓ|⟩ lúc gãy, va chạm/đợt (so 5.12), năng lượng/va chạm, thời gian lưu (CV).
Q2  Darwinism với bản ghi góc: xung Carnot truyền dọc ∇f = (−κ cosθ, 1) nên mỗi nucleon nhận xung ghi cosθ. Mỗi tích tắc
    τ = T0/10: tiến hóa e^{−iHτ} rồi ancilla qubit quay góc ε cosθ. I(R:F) theo |F| cho cặp bó sóng ràng buộc θ0 = 0/π
    và cặp sóng dừng ±ℓ (trải đều theo góc).
Q3  HEOM ở liên kết mạnh λ = 1 (A7 chưa hội tụ): tích phân thời gian tới trạng thái dừng cho K = 6–8, độ sâu 3–4.
Chạy: python3 giai_so_part2_mo2.py [Q1|Q2|Q3]
"""
import itertools
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

import giai_so_part2_mo as P

T0 = 2 * math.pi


# ======================================================================= Q1
def q1_model(pstar):
    m, H, Pout = P.a1_model()
    E, V = np.linalg.eigh(H)
    po = np.array([float(np.real(V[:, j].conj() @ Pout @ V[:, j])) for j in range(len(E))])
    free = po > pstar
    ue = V.conj().T @ m["u"] @ V
    A = math.sqrt(0.5) * np.triu(ue, 1) * free[None, :]            # chỉ từ trạng thái không liên kết j
    Lc = V @ A @ V.conj().T
    return m, H, Pout, [Lc, P.pump(m)], E, V, po, free


def _q1_traj(args):
    pstar, seed, Tmax, dt = args
    m, H, Pout, Ls, E, V, po, free = q1_model(pstar)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * dt)) @ np.linalg.inv(Rm)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = V[:, 0].astype(complex)
    t, r = 0.0, rng.random()
    jumps = []
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
                jumps.append((t, Eb - np.real(np.vdot(psi, H @ psi)), lb))
            r = rng.random()
    return jumps


def part_q1(N=24, Tmax=1200.0, dt=0.01):
    print("=== Q1. Gãy kiểu Part I: Carnot chỉ tác dụng lên trạng thái không liên kết (giếng R = 3, ℓ ≤ 6, n ≤ 12) ===")
    m, H, Pout, Ls, E, V, po, free = q1_model(0.1)
    ll, nn = P.labels(m, V)
    print("phân bố P_ngoài trên trạng thái riêng: " +
          ", ".join(f"{a}–{b}: {np.sum((po>=a)&(po<b))}" for a, b in ((0, .01), (.01, .1), (.1, .3), (.3, .6), (.6, 1.01))))
    for pstar in (0.1, 0.3):
        fr = po > pstar
        lowE = np.argsort(np.where(fr, E, np.inf))[:5]
        print(f"p* = {pstar}: {fr.sum()} trạng thái không liên kết; thấp nhất E = " +
              ", ".join(f"{E[j]:.2f} (|ℓ| {ll[j]:.1f}, n {nn[j]:.1f})" for j in lowE))
    with Pool() as pool:
        for pstar in (0.1, 0.3):
            res = pool.map(_q1_traj, [(pstar, 6000 + i, Tmax, dt) for i in range(N)])
            sizes, dws, eJ, lon, ens = [], [], [], [], []
            for jl in res:
                jl = [j for j in jl if j[0] > Tmax / 6]
                cur, ce, last = [], 0.0, None
                for tj, e, lb in jl:
                    eJ.append(e)
                    if last is None or tj - last >= T0:
                        if cur:
                            sizes.append(len(cur)); ens.append(ce)
                        if last is not None:
                            dws.append((tj - last) / T0)
                        cur, ce = [], 0.0
                        lon.append(lb)
                    cur.append(tj); ce += e; last = tj
            s, d_, e_ = np.array(sizes), np.array(dws), np.array(eJ)
            mean = s.mean() if len(s) else float('nan')
            print(f"-- p* = {pstar}: {len(s)} đợt; va chạm/đợt TB {mean:.2f} (hình học cùng TB: P(1) = {1/mean:.0%}; đo "
                  f"P(1) = {np.mean(s==1):.0%}); ΔE/va chạm TB {e_.mean():.3f} ħω0 (trung vị {np.median(e_):.3f}); "
                  f"ΔE/đợt {np.mean(ens):.3f}")
            print(f"   ⟨|ℓ|⟩ lúc khởi phát đợt TB {np.mean(lon):.2f} (trung vị {np.median(lon):.2f}); thời gian lưu giữa đợt "
                  f"TB {d_.mean():.2f} T0, CV {d_.std()/d_.mean():.2f}", flush=True)
    print("(Part I: 5.12 va chạm/đợt, 0.40–0.46 ħω0/va chạm, θ̇_br = 3.205, lưu 4.5–6 T0 với CV 1.0–1.1)\n")


# ======================================================================= Q2
def part_q2(NT=12, eps=0.5):
    m = P.rotor_osc(3, 2)
    H = m["H"]
    E, V = np.linalg.eigh(H)
    ls = m["ls"]
    d = H.shape[0]
    n = np.arange(m["nmax"] + 1)

    def packet(th0, l0):
        cr = np.exp(-(ls - l0)**2 / (4 * 0.7**2) - 1j * ls * th0); cr /= np.linalg.norm(cr)
        al = ((1 / 3) * math.sin(th0) + 1j * (1 / 3) * l0 * math.cos(th0)) / math.sqrt(2)
        cx = np.array([al**k / math.sqrt(math.factorial(k)) for k in n]) * math.exp(-abs(al)**2 / 2)
        return np.kron(cr, cx / np.linalg.norm(cx))
    L2 = m["Lop"] @ m["Lop"]
    dub = [j for j in range(len(E)) if abs(np.real(V[:, j].conj() @ L2 @ V[:, j]) - 4) < 0.6][:2]
    pairs = {"C: bó sóng ràng buộc θ0 = 0 / π (ℓ0 = 2)": (packet(0, 2), packet(math.pi, 2)),
             "E: cặp sóng dừng ±2 (cos 2θ / sin 2θ)": (V[:, dub[0]], V[:, dub[1]])}
    cosS = np.kron(m["cos"], np.eye(m["nmax"] + 1))
    cw, cv = np.linalg.eigh(cosS)
    tau = T0 / 10
    Uh = (V * np.exp(-1j * E * tau)) @ V.conj().T
    # ancilla |0⟩ → cos(ε c)|0⟩ + sin(ε c)|1⟩ với c = trị riêng cosθ
    K0 = (cv * np.cos(eps * cw)) @ cv.conj().T
    K1 = (cv * np.sin(eps * cw)) @ cv.conj().T

    def ent(M):
        w = np.linalg.eigvalsh(M); w = w[w > 1e-13]
        return float(-np.sum(w * np.log2(w)))
    print(f"=== Q2. Darwinism với bản ghi góc cosθ (ε = {eps}, {NT} bản ghi, τ = T0/10; S ℓ ≤ 3, n ≤ 2) ===")
    for name, (A_, B_) in pairs.items():
        psi = np.zeros((d,) + (1,) + (2,), complex)        # (S, bản ghi…, R)
        psi[:, 0, 0] = A_ / math.sqrt(2); psi[:, 0, 1] = B_ / math.sqrt(2)
        psi = psi[:, 0, :]                                  # (S, R)
        recs = []
        cur = psi                                           # trục: S, q1..qk, R
        for k in range(NT):
            cur = np.tensordot(Uh, cur, axes=(1, 0))
            a0 = np.tensordot(K0, cur, axes=(1, 0)); a1 = np.tensordot(K1, cur, axes=(1, 0))
            cur = np.stack([a0, a1], axis=-2)              # chèn qubit mới trước trục R
        # cur có trục (S, q1, …, qNT, R)
        Rax = NT + 1

        def red(keep):
            other = [ax for ax in range(NT + 2) if ax not in keep]
            Mx = np.transpose(cur, keep + other).reshape(int(np.prod([cur.shape[k] for k in keep])), -1)
            return ent(Mx @ Mx.conj().T) if Mx.shape[0] <= Mx.shape[1] else ent(Mx.conj().T @ Mx)
        SR = red([Rax]); rng = random.Random(3); Iv = []
        for m_ in range(NT + 1):
            frs = list(itertools.combinations(range(1, NT + 1), m_)); rng.shuffle(frs)
            Iv.append(np.mean([SR + (red(list(F)) if F else 0.0) - red([Rax] + list(F)) for F in frs[:4]]))
        md = next((k for k, v in enumerate(Iv) if v >= 0.9), None)
        rho_RS = None
        ov = abs(np.vdot(A_, B_))
        print(f"   {name}: ⟨A|B⟩ = {ov:.3f}; I(R:F) (bit) |F| = 0…{NT}: " + " ".join(f"{v:.2f}" for v in Iv) +
              f";  R_0.1 = {NT/md if md else 0:.1f}", flush=True)
    print("(bình nguyên I ≈ 1 bit ở mảnh nhỏ = thông tin cổ điển về nhánh được nhân bản; 2 bit ở toàn bộ = phần lượng tử)\n")


# ======================================================================= Q3
def heom_time(Hs, Q, Lp, lam, gD, kT, K, depth, rho0, dt=0.02, tmax=400.0, tol=1e-7):
    beta = 1 / kT
    nus = [gD] + [2 * math.pi * k / beta for k in range(1, K + 1)]
    cs = [lam * gD * (1 / math.tan(beta * gD / 2) - 1j)] + [4 * lam * gD / beta * nk / (nk**2 - gD**2) for nk in nus[1:]]
    tail = 4 * lam * gD / beta * sum(1 / (nk**2 - gD**2) for nk in [2 * math.pi * k / beta for k in range(K + 1, 20000)])
    idx = [n for n in itertools.product(range(depth + 1), repeat=K + 1) if sum(n) <= depth]
    pos = {n: i for i, n in enumerate(idx)}
    NA = len(idx)
    up = np.array([[pos.get(tuple(n[:k] + (n[k] + 1,) + n[k + 1:]), -1) for k in range(K + 1)] for n in idx])
    dn = np.array([[pos[tuple(n[:k] + (n[k] - 1,) + n[k + 1:])] if n[k] > 0 else -1 for k in range(K + 1)] for n in idx])
    gsum = np.array([sum(n[k] * nus[k] for k in range(K + 1)) for n in idx])
    nk = np.array(idx, float)
    LdL = Lp.conj().T @ Lp; Q2 = Q @ Q

    def rhs(R):
        out = -1j * (Hs @ R - R @ Hs) + Lp @ R @ Lp.conj().T - 0.5 * (LdL @ R + R @ LdL) - gsum[:, None, None] * R
        out -= tail * (Q2 @ R - 2 * Q @ R @ Q + R @ Q2)
        for k in range(K + 1):
            mu = up[:, k] >= 0
            Rp = R[up[mu, k]]; out[mu] += -1j * (Q @ Rp - Rp @ Q)
            md = dn[:, k] >= 0
            Rm = R[dn[md, k]]
            out[md] += -1j * nk[md, k][:, None, None] * (cs[k] * (Q @ Rm) - np.conj(cs[k]) * (Rm @ Q))
        return out
    D = Hs.shape[0]
    R = np.zeros((NA, D, D), complex); R[0] = rho0
    t, last = 0.0, R[0].copy()
    while t < tmax:
        for _ in range(int(5 / dt)):
            k1 = rhs(R); k2 = rhs(R + dt / 2 * k1); k3 = rhs(R + dt / 2 * k2); k4 = rhs(R + dt * k3)
            R = R + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += 5
        dr = float(np.linalg.norm(R[0] - last)) / 5; last = R[0].copy()
        if dr < tol:
            break
    r0 = (R[0] + R[0].conj().T) / 2
    return r0 / np.trace(r0).real, dr, t, NA


def part_q3(lam=1.0, gD=2.0, kT=0.1):
    m = P.rotor_osc(4, 1)
    H, u = m["H"], m["u"]
    Lp = P.pump(m)
    E, V = np.linalg.eigh(H)
    ue = V.conj().T @ u @ V
    J = lambda w: 2 * lam * gD * w / (w * w + gD * gD)
    w = E[:, None] - E[None, :]
    nb = lambda x: 1 / np.expm1(np.clip(x, 1e-12, None) / kT)
    Adn = np.where(w < 0, np.sqrt(np.abs(2 * J(-w) * (nb(-w) + 1))) * ue, 0)
    Aup = np.where(w > 0, np.sqrt(np.abs(2 * J(w) * nb(w))) * ue, 0)
    rD = P.steady(H, [V @ Adn @ V.conj().T, V @ Aup @ V.conj().T, Lp])
    ev = lambda O, r: float(np.real(np.trace(O @ r)))
    print(f"=== Q3. HEOM liên kết mạnh λ = {lam} (γ_D = {gD}, k_BT = {kT}), tích phân thời gian tới dừng; ℓ ≤ 4, n ≤ 1 ===")
    print(f"{'phương pháp':>28} {'⟨ℓ⟩':>7} {'⟨H⟩':>7} {'Tr ρ²':>7} {'min eig':>9} {'‖dρ/dt‖':>9} {'t':>5}")
    print(f"{'Davies':>28} {ev(m['Lop'],rD):7.3f} {ev(H,rD):7.3f} {np.real(np.trace(rD@rD)):7.4f} "
          f"{np.linalg.eigvalsh(rD).min():9.1e} {'—':>9} {'—':>5}", flush=True)
    for K, depth in ((6, 3), (8, 3), (6, 4)):
        rH, dr, t, NA = heom_time(H, u, Lp, lam, gD, kT, K, depth, rD)
        print(f"{f'HEOM K={K}, độ sâu {depth} ({NA} ADO)':>28} {ev(m['Lop'],rH):7.3f} {ev(H,rH):7.3f} "
              f"{np.real(np.trace(rH@rH)):7.4f} {np.linalg.eigvalsh(rH).min():9.1e} {dr:9.1e} {t:5.0f}", flush=True)
    print("(khởi động từ trạng thái dừng Davies; dừng khi ‖dρ/dt‖ < 1e-7)\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["Q1", "Q2", "Q3"]:
        {"Q1": part_q1, "Q2": part_q2, "Q3": part_q3}[name.upper()]()
