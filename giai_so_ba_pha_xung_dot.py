"""Ghép E_quá khứ + E_hiện tại + E_tương lai: các xung đột có thể xảy ra, và người quan sát (bể Carnot) làm đồng hồ
Page–Wootters (HocThuyetE (7).pdf; mô hình lượng tử hiệu chỉnh tốt nhất của Part II: giếng ràng buộc R, Carnot chỉ trên
trạng thái không liên kết p* = 0.3, bơm Γ0 = 2 (ℓ0 = 4), bể rô-to Γt = 1 ở k_BT = 1 MeV — khớp lưu 5.0 T0, CV 0.92,
ngưỡng 3.59 của Part I).

Mỗi quỹ đạo MCWF = một nhánh (một bản ghi môi trường) của trạng thái lịch sử. Bản ghi: Carnot (c), bơm (p), bể (b).
Pha theo bản ghi Carnot: Hiện tại = trong đợt (các bản ghi Carnot cách nhau < T0), Tương lai = T0 sau đợt, Quá khứ = còn lại.
(A) Dòng thời gian ba pha: tỉ lệ thời gian, độ dài từng pha, chuyển tiếp.
(B) Xung đột: B1 chồng chập (tỉ lệ nhánh ở mỗi pha tại một thời điểm — trước khi đọc bản ghi); B2 hai người quan sát
    mỗi người nhận một bản ghi Carnot với xác suất p (dư thừa không hoàn hảo, R_δ ≈ 1.5) — tỉ lệ thời gian họ gán pha
    khác nhau, theo loại; B3 đồng hồ không phân biệt: hai mốc pha liên tiếp rơi vào cùng số đọc.
(C) Người quan sát làm đồng hồ PW: đồng hồ = bộ đếm bản ghi (chỉ Carnot; mọi bản ghi). Tiêu chí: đơn điệu, độ phân giải
    σ_t(k), tính độc lập với hệ (tương quan tốc độ tích tắc ↔ trạng thái hệ), động lực học thấy được giữa hai tích tắc.
Chạy: python3 giai_so_ba_pha_xung_dot.py   (≈ 10 phút, 8 lõi)
"""
import math
import random
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
TMAX, DT, WARM = 2400.0, 0.01, 400.0


def _traj(seed):
    import giai_so_part2_mo4 as R
    m, H, Ls, V = R.r1_model(0.1, Gt=1.0, G0=2.0, l0=4.0)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * DT)) @ np.linalg.inv(Rm)
    Labs = np.abs(m["Lop"])
    rng = random.Random(seed)
    psi = V[:, 0].astype(complex)
    t, r, k = 0.0, rng.random(), 0
    rec, ts, lv = [], [], []
    names = "cpbb"
    while t < TMAX:
        psi = U @ psi; t += DT; k += 1
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            j = min(int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random())), len(Ls) - 1)
            psi = Ls[j] @ psi; psi /= np.linalg.norm(psi)
            rec.append((t, names[j], np.real(np.vdot(psi, Labs @ psi))))
            r = rng.random()
        if k % 10 == 0:
            p = psi / np.linalg.norm(psi)
            ts.append(t); lv.append(np.real(np.vdot(p, Labs @ p)))
    return rec, np.array(ts), np.array(lv)


def phases(ctimes, grid):
    """Nhãn pha trên lưới thời gian từ danh sách thời điểm bản ghi Carnot: 0 Quá khứ, 1 Hiện tại, 2 Tương lai."""
    lab = np.zeros(len(grid), int)
    if len(ctimes) == 0:
        return lab
    bursts, cur = [], [ctimes[0]]
    for a, b in zip(ctimes, ctimes[1:]):
        if b - a < T0:
            cur.append(b)
        else:
            bursts.append((cur[0], cur[-1])); cur = [b]
    bursts.append((cur[0], cur[-1]))
    for s, e in bursts:
        lab[(grid >= s) & (grid <= e)] = 1
        lab[(grid > e) & (grid <= e + T0) & (lab == 0)] = 2
    return lab


def main(N=16):
    with Pool() as pool:
        res = pool.map(_traj, range(23000, 23000 + N))
    grid = np.arange(WARM, TMAX, 0.05)
    labs, ctl = [], []
    for rec, ts, lv in res:
        c = np.array([x[0] for x in rec if x[1] == "c"])
        labs.append(phases(c, grid)); ctl.append(c)
    L = np.array(labs)
    names = ["Quá khứ", "Hiện tại", "Tương lai"]
    print(f"=== A. Dòng thời gian ba pha ghép ({N} nhánh × {(TMAX-WARM)/T0:.0f} T0) ===")
    for i in range(3):
        runs = []
        for lab in L:
            d = np.diff(np.concatenate([[0], (lab == i).astype(int), [0]]))
            st, en = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
            runs += list((en - st) * 0.05 / T0)
        print(f"   {names[i]:>10}: {np.mean(L == i):6.1%} thời gian; độ dài một lần TB {np.mean(runs):.2f} T0 "
              f"(trung vị {np.median(runs):.2f})")
    trans = np.zeros((3, 3), int)
    for lab in L:
        ch = np.nonzero(np.diff(lab))[0]
        for i in ch:
            trans[lab[i], lab[i + 1]] += 1
    print("   ma trận chuyển pha (hàng → cột: QK, HT, TL):")
    for i in range(3):
        print(f"      {names[i]:>10} → " + "  ".join(f"{trans[i, j]:5d}" for j in range(3)))
    print("   (chu trình đúng: QK → HT → TL → QK; TL → HT = gãy lại trước khi kịp khóa ổn định)\n")

    print("=== B. Các xung đột ===")
    occ = np.array([np.mean(L[:, j] == i) for i in range(3) for j in range(L.shape[1])]).reshape(3, -1)
    print(f"B1 chồng chập (trước khi đọc bản ghi): tại một thời điểm, tỉ lệ nhánh ở QK/HT/TL = "
          f"{occ[0].mean():.2f}/{occ[1].mean():.2f}/{occ[2].mean():.2f}; thời điểm có ≥ 2 pha cùng hiện diện trong các nhánh: "
          f"{np.mean((occ > 0).sum(axis=0) >= 2):.0%}, đủ cả 3 pha: {np.mean((occ > 0).sum(axis=0) == 3):.0%}")
    rng = np.random.default_rng(5)
    print("B2 hai người quan sát, mỗi người nhận một bản ghi Carnot với xác suất p (độc lập):")
    print(f"   {'p':>5} {'thời gian bất đồng':>19} {'A:QK–B:HT':>10} {'A:HT–B:TL':>10} {'A:QK–B:TL':>10} {'số đợt A/B so thật':>19}")
    for pr in (1.0, 0.9, 0.7, 0.5, 0.3):
        dis, t01, t12, t02, nb = [], [], [], [], []
        for c in ctl:
            a = c[rng.random(len(c)) < pr]; b = c[rng.random(len(c)) < pr]
            la, lb = phases(a, grid), phases(b, grid)
            dis.append(np.mean(la != lb))
            t01.append(np.mean(((la == 0) & (lb == 1)) | ((la == 1) & (lb == 0))))
            t12.append(np.mean(((la == 1) & (lb == 2)) | ((la == 2) & (lb == 1))))
            t02.append(np.mean(((la == 0) & (lb == 2)) | ((la == 2) & (lb == 0))))
            nt = 1 + np.sum(np.diff(c) >= T0) if len(c) else 0
            na = 1 + np.sum(np.diff(a) >= T0) if len(a) else 0
            nb.append(na / max(nt, 1))
        print(f"   {pr:5.1f} {np.mean(dis):19.1%} {np.mean(t01):10.1%} {np.mean(t12):10.1%} {np.mean(t02):10.1%} "
              f"{np.mean(nb):19.2f}")
    ties, tot = 0, 0
    for (rec, ts, lv), c in zip(res, ctl):
        allt = np.array([x[0] for x in rec])
        lab = phases(c, grid)
        ch = grid[np.nonzero(np.diff(lab))[0]]
        for a, b in zip(ch, ch[1:]):
            tot += 1
            ties += np.sum((allt > a) & (allt <= b)) == 0
    print(f"B3 đồng hồ = MỌI bản ghi: {ties}/{tot} cặp mốc pha liên tiếp không có tích tắc nào ở giữa ({ties/max(tot,1):.1%}) — "
          "đồng hồ không phân biệt được\n")

    print("=== C. Người quan sát (bể Carnot) làm đồng hồ Page–Wootters ===")
    for label, keep in (("chỉ bản ghi Carnot", "c"), ("mọi bản ghi (c+p+b)", "cpb")):
        tk, lafter, rates, lvs = [], [], [], []
        for rec, ts, lv in res:
            rr = [x for x in rec if x[1] in keep and x[0] > WARM]
            tk.append(np.array([x[0] for x in rr]) - WARM)
            lafter += [x[2] for x in rr]
            win = np.arange(WARM, TMAX - T0, T0)
            cnt = np.histogram([x[0] for x in rr], bins=np.append(win, win[-1] + T0))[0]
            mean_l = np.array([lv[(ts >= a) & (ts < a + T0)].mean() for a in win])
            rates.append(cnt); lvs.append(mean_l)
        rate = np.mean([len(x) for x in tk]) / ((TMAX - WARM) / T0)
        sig = []
        for kk in (5, 20, 80):
            v = [x[kk] / T0 for x in tk if len(x) > kk]
            sig.append((kk, np.mean(v), np.std(v)))
        mono = all(np.all(np.diff(x) > 0) for x in tk)
        corr = np.corrcoef(np.concatenate(rates), np.concatenate(lvs))[0, 1]
        gaps = np.concatenate([np.diff(x) for x in tk]) / T0
        print(f"-- đồng hồ = {label}: tốc độ {rate:.2f} tích tắc/T0; đơn điệu: {'có' if mono else 'KHÔNG'}; "
              f"khoảng giữa tích tắc TB {gaps.mean():.2f} T0, CV {gaps.std()/gaps.mean():.2f}")
        print("   độ phân giải: " + "; ".join(f"tích tắc {kk}: t = {m:.1f} ± {s:.1f} T0 (tương đối {s/m:.0%})"
                                            for kk, m, s in sig))
        print(f"   phụ thuộc vào hệ: tương quan (số tích tắc trong 1 T0) ↔ ⟨|ℓ|⟩ = {corr:+.2f}; "
              f"⟨|ℓ|⟩ ngay sau tích tắc: TB {np.mean(lafter):.2f} ± {np.std(lafter):.2f}")
    print()


if __name__ == "__main__":
    main()
