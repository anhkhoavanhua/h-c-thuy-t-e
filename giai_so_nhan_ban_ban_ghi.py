"""Tăng độ nhân bản bản ghi Carnot bằng dòng thác exciton (ánh xạ Ericson của Part I, HocThuyetE (7).pdf).

Cơ chế: va chạm Carnot tạo một cặp hạt–lỗ (n = 2 nucleon đổi trạng thái); va chạm nucleon–nucleon bên trong tiếp tục tạo
cặp mới với tốc độ Williams λ₊(n, E) = (2π/ħ)|M|² g³E²/(2(n+1)), |M|² = K/(A³E), K = 400 MeV³, tới n_eq = √(1.39 g E).
Mỗi nucleon bị kích thích là một bản sao của sự kiện "đã có va chạm". Môi trường đọc được: N_act nucleon gần mặt Fermi
(N_act = 30); nhiễu nền: mỗi nucleon không tham gia dòng thác vẫn bị kích thích nhiệt với xác suất q.
Đo: phân bố n(t); I(R:F) theo kích thước mảnh m (R = qubit tham chiếu: có/không va chạm, tiên nghiệm 1/2); R_δ = N_act/m_δ
(I ≥ 0.9 bit); xác suất hai người quan sát với hai mảnh rời nhau cùng kết luận đúng.
So với cơ chế hiện tại (một bản ghi/va chạm, mục A8: R_δ ≈ 1.4–1.7).
Chạy: python3 giai_so_nhan_ban_ban_ghi.py   (≈ 1 phút)
"""
import math

import numpy as np

HBAR = 6.582119569e-22
T0 = 4.1357e-22
A = 70
g = A / 13.0
NACT = 30


def lam_plus(n, E, K=400.0):
    return 2 * math.pi / HBAR * (K / (A**3 * E)) * g**3 * E**2 / (2 * (n + 1)) * T0   # mỗi T0


def n_distribution(E, t, nsamp=20000, rng=None):
    """Gillespie: n(t) từ n = 2 lúc va chạm, mỗi bước n → n + 2, chặn ở n_eq."""
    neq = int(round(math.sqrt(1.39 * g * E)))
    neq -= neq % 2
    out = np.empty(nsamp, int)
    for s in range(nsamp):
        n, tt = 2, 0.0
        while n < neq:
            tt += rng.exponential(1 / lam_plus(n, E))
            if tt > t:
                break
            n += 2
        out[s] = min(n, NACT)
    return out, neq


def hyper(j, N, K, m):
    if j < max(0, m - (N - K)) or j > min(K, m):
        return 0.0
    return math.comb(K, j) * math.comb(N - K, m - j) / math.comb(N, m)


def pk_given(m, nvals, q):
    """P(K_F = k) cho mảnh m nucleon: nhánh A (có va chạm, n phân bố nvals) và nhánh B (chỉ nhiễu)."""
    PA = np.zeros(m + 1)
    ns, cnt = np.unique(nvals, return_counts=True)
    w = cnt / cnt.sum()
    for n, wn in zip(ns, w):
        for j in range(0, min(n, m) + 1):
            hj = hyper(j, NACT, n, m)
            if hj == 0:
                continue
            for b in range(0, m - j + 1):
                PA[j + b] += wn * hj * math.comb(m - j, b) * q**b * (1 - q)**(m - j - b)
    PB = np.array([math.comb(m, k) * q**k * (1 - q)**(m - k) for k in range(m + 1)])
    return PA, PB


def mutual_info(PA, PB):
    P = 0.5 * (PA + PB)
    I = 0.0
    for pa, pb, p in zip(PA, PB, P):
        for px in (pa, pb):
            if px > 0:
                I += 0.5 * px * math.log2(px / p)
    return I


def agreement(m, nvals, q, rng, trials=20000):
    """Hai người quan sát, mỗi người m nucleon RỜI NHAU, kết luận MAP; xác suất cùng kết luận và đúng."""
    PA, PB = pk_given(m, nvals, q)
    dec = PA > PB                                             # k → kết luận "có va chạm"
    same = correct = 0
    for _ in range(trials):
        ev = rng.random() < 0.5
        exc = np.zeros(NACT, bool)
        if ev:
            n = nvals[rng.integers(len(nvals))]
            exc[rng.choice(NACT, n, replace=False)] = True
        exc |= rng.random(NACT) < q
        perm = rng.permutation(NACT)
        ka, kb = exc[perm[:m]].sum(), exc[perm[m:2 * m]].sum()
        da, db = dec[ka], dec[kb]
        same += da == db
        correct += (da == ev) and (db == ev)
    return same / trials, correct / trials


def main():
    rng = np.random.default_rng(9)
    print("=== Dòng thác exciton: số nucleon mang bản sao sự kiện theo thời gian ===")
    print(f"A = {A}, g = {g:.2f} MeV⁻¹, N_act = {NACT}; λ₊(n=2, E=4.1) = {lam_plus(2, 4.1):.2f}/T0, "
          f"λ₊(n=2, E=20) = {lam_plus(2, 20):.2f}/T0")
    cases = {"một va chạm (E = 4.1 MeV)": 4.1, "cả đợt (E = 20 MeV)": 20.0}
    dist = {}
    for name, E in cases.items():
        row = []
        for t in (0.1, 0.5, 1.0, 2.0, 5.0):
            nv, neq = n_distribution(E, t, 4000, rng)
            dist[(E, t)] = nv
            row.append(f"t={t:g}T0: ⟨n⟩={nv.mean():.1f}")
        print(f"   {name} (n_eq = {neq}): " + "; ".join(row))
    print()

    print("=== Độ nhân bản R_δ = N_act/m_δ (I(R:F) ≥ 0.9 bit) ===")
    print(f"{'trường hợp':>28} {'q nền':>6} {'t (T0)':>7} {'I(1 nucleon)':>13} {'m_δ':>5} {'R_δ':>6} {'I(m=5)':>7}")
    for name, E in cases.items():
        for q in (0.02, 0.1, 0.2):
            for t in (0.5, 2.0, 5.0):
                nv = dist[(E, t)]
                Is = [mutual_info(*pk_given(m, nv, q)) for m in range(1, NACT + 1)]
                md = next((i + 1 for i, v in enumerate(Is) if v >= 0.9), None)
                R = NACT / md if md else 0.0
                print(f"{name:>28} {q:6.2f} {t:7.1f} {Is[0]:13.3f} {md if md else '—':>5} {R:6.1f} {Is[4]:7.3f}")
    print("(cơ chế hiện tại — một bản ghi mỗi va chạm, mục A8: R_δ ≈ 1.4–1.7)\n")

    print("=== Hai người quan sát, mỗi người đọc m nucleon riêng: cùng kết luận / cùng đúng ===")
    print(f"{'trường hợp':>28} {'q':>5} {'t':>5} " + " ".join(f"{'m='+str(m):>14}" for m in (1, 3, 5, 10)))
    for name, E in cases.items():
        for q in (0.02, 0.2):
            t = 2.0
            nv = dist[(E, t)]
            cells = []
            for m in (1, 3, 5, 10):
                s, c = agreement(m, nv, q, rng, 6000)
                cells.append(f"{s:6.1%}/{c:6.1%}")
            print(f"{name:>28} {q:5.2f} {t:5.1f} " + " ".join(f"{x:>14}" for x in cells))
    print("(ô: P(cùng kết luận) / P(cả hai đúng); so B2 bài trước: bất đồng tới 13% thời gian)\n")


if __name__ == "__main__":
    main()
