"""Hai vũ trụ, mỗi vũ trụ một động cơ E_dual, chạm nhau ở tầng chân không qua kênh trao đổi.

Mỗi động cơ: rô-to ℓ ≤ 4 ⊗ hạt n ≤ 1 (18 chiều, đủ mức gãy ℓ = 4), H Bảng I (giai_so_u_hermite_drude.basis),
Carnot L_C = √0.5 V triu(V†uV) V†, bơm √(0.2 a) · (bậc thang ℓ → ℓ+1 có bão hòa F).
Ghép (cường độ κ): L_{A→B} = √κ · L̂_C^A ⊗ Û^B (A tiêu tán đúng lúc B được nâng), L_{B→A} tương tự.
Hệ ghép 324 chiều, giải bằng MCWF (mỗi quỹ đạo = một nhánh của trạng thái lịch sử chung).

Đo: năng lượng vượt chân không ΔE_A, ΔE_B; tốc độ bản ghi Carnot riêng r_A, r_B (đồng hồ nội tại);
tốc độ trao đổi (thứ duy nhất mỗi vũ trụ "thấy" của vũ trụ kia); tỉ số tốc độ thời gian r_A/r_B.
Chạy: python3 giai_so_hai_vu_tru_e_dual.py   (≈ 5–10 phút, 8 lõi)
"""
import math
import random
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
DT, TMAX, WARM = 0.01, 800.0, 100.0


def engine():
    import giai_so_u_hermite_drude as UH
    ls, cos_t, sin_t, pth, H, u = UH.basis(lmax=4, nmax=1)
    E, V = np.linalg.eigh(H)
    Lc = V @ (math.sqrt(0.5) * np.triu(V.conj().T @ u @ V, 1)) @ V.conj().T
    up = np.diag(np.ones(len(ls) - 1), -1).astype(complex)
    F = np.diag(np.sqrt(np.maximum(0.0, 1 - ls.astype(float)**2 / 16.0)))
    Up = np.kron(up @ F, np.eye(2))
    return H, Lc, Up, E[0]


def build(aA, aB, kap):
    H1, Lc, Up, E0 = engine()
    d = H1.shape[0]; I = np.eye(d)
    H = np.kron(H1, I) + np.kron(I, H1)
    Ls, tags = [], []
    def add(L, tag):
        if np.abs(L).max() > 0:
            Ls.append(L); tags.append(tag)
    add(np.kron(Lc, I), 'cA'); add(np.kron(I, Lc), 'cB')
    add(math.sqrt(0.2 * aA) * np.kron(Up, I), 'pA'); add(math.sqrt(0.2 * aB) * np.kron(I, Up), 'pB')
    if kap > 0:
        add(math.sqrt(kap) * np.kron(Lc, Up) / math.sqrt(0.5), 'xAB')     # A tiêu tán → B được nâng
        add(math.sqrt(kap) * np.kron(Up, Lc) / math.sqrt(0.5), 'xBA')
    return H, Ls, tags, np.kron(H1, I), np.kron(I, H1), E0, d


def _traj(args):
    aA, aB, kap, seed = args
    H, Ls, tags, HA, HB, E0, d = build(aA, aB, kap)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, R = np.linalg.eig(Heff)
    U = R @ np.diag(np.exp(-1j * w * DT)) @ np.linalg.inv(R)
    rng = random.Random(seed)
    psi = np.zeros(d * d, complex); psi[0] = 1.0                 # cả hai ở trạng thái cơ sở thấp nhất
    t, r = 0.0, rng.random()
    cnt = {tg: 0 for tg in tags}; eA = eB = 0.0; ns = 0
    while t < TMAX:
        psi = U @ psi; t += DT
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            j = min(int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random())), len(Ls) - 1)
            psi = Ls[j] @ psi; psi /= np.linalg.norm(psi)
            if t > WARM:
                cnt[tags[j]] += 1
            r = rng.random()
        if t > WARM and int(t / DT) % 20 == 0:
            p = psi / np.linalg.norm(psi)
            eA += np.vdot(p, HA @ p).real; eB += np.vdot(p, HB @ p).real; ns += 1
    span = (TMAX - WARM) / T0
    return dict(case=(aA, aB, kap), cnt={k: v / span for k, v in cnt.items()},
                dEA=eA / ns - E0, dEB=eB / ns - E0)


CASES = [
    # (a_A, a_B, κ)
    (1.0, 1.0, 0.0), (1.0, 0.3, 0.0), (1.0, 3.0, 0.0),              # tách rời: năng lượng khác ⟹ thời gian khác
    (1.0, 0.3, 0.05), (1.0, 0.3, 0.2),                               # chạm nhau ở tầng chân không
    (0.0, 1.0, 0.0), (0.0, 1.0, 0.05), (0.0, 1.0, 0.2),              # A đóng băng (không bơm); B có làm A tan băng?
    (0.0, 0.0, 0.2),                                                 # cả hai đóng băng
]


def main(N=16):
    with Pool(8) as p:
        res = p.map(_traj, [(a, b, k, 9000 + 37 * i + s) for i, (a, b, k) in enumerate(CASES) for s in range(N)])
    print("Hai vũ trụ E_dual (mỗi bên 18 chiều, hệ ghép 324 chiều, MCWF; "
          f"{N} nhánh × {(TMAX - WARM) / T0:.0f} T0 mỗi trường hợp)")
    print(" a_A  a_B   κ    | ΔE_A   ΔE_B (ħω0) | r_A    r_B  (bản ghi Carnot riêng /T0) | r_A/r_B | "
          "trao đổi A→B  B→A (/T0) = điều mỗi bên 'thấy' của bên kia")
    for c in CASES:
        rs = [x for x in res if x['case'] == c]
        g = lambda k: np.mean([x['cnt'].get(k, 0.0) for x in rs])
        se = lambda k: np.std([x['cnt'].get(k, 0.0) for x in rs]) / math.sqrt(len(rs))
        rA = g('cA') + g('xAB'); rB = g('cB') + g('xBA')          # mỗi lần tiêu tán = một bản ghi của chính vũ trụ đó
        dA = np.mean([x['dEA'] for x in rs]); dB = np.mean([x['dEB'] for x in rs])
        ratio = rA / rB if rB > 0 else float('inf')
        print(f" {c[0]:3.1f}  {c[1]:3.1f}  {c[2]:4.2f} | {dA:5.2f}  {dB:5.2f}       | {rA:5.3f} {rB:5.3f} "
              f"(±{se('cA'):.3f}, ±{se('cB'):.3f})          | {ratio:6.2f}  | {g('xAB'):6.3f}  {g('xBA'):6.3f}")


if __name__ == '__main__':
    main()
