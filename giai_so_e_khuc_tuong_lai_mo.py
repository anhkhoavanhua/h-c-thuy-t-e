"""Xung đột cả ba pha cùng lúc và "E khúc": tương lai trong Page–Wootters có bị định đoạt không?

"Định đoạt" được đo bằng ĐỘ MỞ của tương lai: entropy có điều kiện của bản ghi kế tiếp khi đã biết TOÀN BỘ bản ghi quá
khứ trên một nhánh, h = ⟨−Σ_j p_j ln p_j⟩ với p_j = ‖K_j ψ_nhánh‖² (nat/T0). h = 0 ⟺ tương lai đã quyết định.
"E khúc" = hiện tại có bề dày: mép sau (bản ghi đã viết, cố định) và mép trước (bản ghi sắp tới, còn mở);
bề dày ε* = 1/h là thời gian để tương lai mở thêm 1 nat.

  A  Mô hình Part III (S 36 chiều, Kraus mỗi nhịp τ = T0/25, ancilla {·, Carnot, bơm}):
     h theo chế độ: PW đóng (không môi trường) / chỉ Carnot / chỉ bơm / cân bằng nhiệt / E_dual (cả ba pha xung đột).
  B  Mô hình MCWF hiệu chỉnh của Part II (giếng gãy, Carnot chỉ trên trạng thái không liên kết, bơm, bể rô-to):
     tách nhánh từ một trạng thái chung — xác suất pha (Quá khứ/Hiện tại/Tương lai) sau Δ; entropy pha H(Δ)
     → chân trời dự đoán; mép sau: pha đã qua do bản ghi cố định.
Chạy: python3 giai_so_e_khuc_tuong_lai_mo.py [A|B]   (A ≈ 1 phút; B ≈ 3–6 phút, 8 lõi)
"""
import math
import random
import sys
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi


# ======================================================================= A
def openness(Ks, psi0, nk, M, rng, late=500):
    """Trung bình theo nhánh của entropy kết quả nhịp kế tiếp (nat/nhịp), trên đoạn cuối `late` nhịp."""
    psis = np.tile(psi0.astype(complex), (M, 1))
    acc = np.zeros(M); cnt = 0
    for k in range(nk):
        outs = np.stack([psis @ K.T for K in Ks], 0)             # (nK, M, d)
        p = np.real(np.einsum('jmd,jmd->jm', outs.conj(), outs))
        p = np.clip(p, 0, None); p /= p.sum(0, keepdims=True)
        if k >= nk - late:
            acc += -(p * np.log(np.where(p > 0, p, 1))).sum(0); cnt += 1
        u = rng.random(M); c = np.cumsum(p, 0)
        j = (u[None, :] > c).sum(0).clip(0, len(Ks) - 1)
        psis = outs[j, np.arange(M)] / np.sqrt(p[j, np.arange(M)])[:, None]
    return acc / cnt                                             # nat/nhịp, mỗi nhánh


def part_a(M=300, nk=1500):
    from giai_so_wdw_pw_e7 import channels, TAU
    from giai_so_kim_chi_gio import kraus, psi0
    from giai_so_part2_mo import davies
    import giai_so_u_hermite_drude as UH
    ls, H, (Lc, Lb) = channels()
    _, _, _, _, _, u = UH.basis(lmax=4, nmax=3)
    rng = np.random.default_rng(5)
    per = T0 / TAU
    regimes = [
        ("PW đóng (chỉ H, không môi trường)", [0 * Lc]),
        ("chỉ Carnot (Quá khứ thắng)", [Lc]),
        ("chỉ bơm (Hiện tại thắng)", [Lb]),
        ("cân bằng nhiệt k_BT = 1 (có tích tắc, không chiều)", davies(H, u, 0.5, kT=1.0)[0]),
        ("E_dual: bơm + Carnot (ba pha xung đột)", [Lc, Lb]),
    ]
    print("=== A. Độ mở của tương lai khi đã biết toàn bộ bản ghi quá khứ (nhánh, Part III S 36 chiều) ===")
    print(f"  {M} nhánh × {nk} nhịp ({nk / per:.0f} T0), đo trên {500 / per:.0f} T0 cuối")
    print("  chế độ                                              h (nat/T0)   ± sai số   bề dày E khúc ε* = 1/h")
    out = {}
    for name, Ls in regimes:
        Ks = kraus(H, Ls, TAU)
        hb = openness(Ks, psi0(), nk, M, rng) * per
        h, se = hb.mean(), hb.std() / math.sqrt(M)
        eps = f"{1 / h:.3f} T0" if h > 1e-6 else "∞ (đóng băng)"
        print(f"  {name:50s}  {h:9.4f}   {se:.4f}     {eps}")
        out[name] = h
    print()
    return out


# ======================================================================= B
DT = 0.01


def _model():
    import giai_so_part2_mo4 as R
    m, H, Ls, V = R.r1_model(0.1, Gt=1.0, G0=2.0, l0=4.0)
    Heff = H - 0.5j * sum(L.conj().T @ L for L in Ls)
    w, Rm = np.linalg.eig(Heff)
    U = Rm @ np.diag(np.exp(-1j * w * DT)) @ np.linalg.inv(Rm)
    return Ls, U, V


def _evolve(psi, t, T, U, Ls, rng):
    """MCWF từ (psi, t) tới t+T. Trả thời điểm các bản ghi Carnot (kênh 0) và trạng thái cuối."""
    r = rng.random(); rec = []
    end = t + T
    while t < end:
        psi = U @ psi; t += DT
        nrm = np.vdot(psi, psi).real
        if nrm < r:
            psi /= math.sqrt(nrm)
            ws = np.array([np.linalg.norm(L @ psi)**2 for L in Ls])
            j = min(int(np.searchsorted(np.cumsum(ws / ws.sum()), rng.random())), len(Ls) - 1)
            psi = Ls[j] @ psi; psi /= np.linalg.norm(psi)
            if j == 0:
                rec.append(t)
            r = rng.random()
    return rec, psi / np.linalg.norm(psi), t


def phase_at(ctimes, t):
    """Nhãn pha tại t khi biết MỌI bản ghi Carnot (cả sau t), đúng quy tắc của ba_pha_xung_dot.phases:
    1 Hiện tại nếu t nằm trong một đợt [s, e] (bản ghi cách nhau < T0); 2 Tương lai nếu t ∈ (e, e + T0]; 0 Quá khứ."""
    if not ctimes:
        return 0
    ct = sorted(ctimes); lab = 0
    cur = [ct[0]]
    bursts = []
    for a, b in zip(ct, ct[1:]):
        if b - a < T0:
            cur.append(b)
        else:
            bursts.append((cur[0], cur[-1])); cur = [b]
    bursts.append((cur[0], cur[-1]))
    for s_, e_ in bursts:
        if s_ <= t <= e_:
            return 1
        if e_ < t <= e_ + T0:
            lab = 2
    return lab


def _fork(args):
    seed, n_start, n_child, horizon = args
    rng = random.Random(seed)
    Ls, U, V = _model()
    psi = V[:, 0].astype(complex)
    hist, psi, t = _evolve(psi, 0.0, 300.0, U, Ls, rng)        # làm nóng ⟹ một nhánh "hiện tại" ngẫu nhiên
    res = []
    for s in range(n_start):
        rec, psi, t = _evolve(psi, t, 40.0 + 60.0 * rng.random(), U, Ls, rng)   # trôi tới thời điểm bắt đầu mới
        hist += rec
        now = t
        ph_now = 0 if (not hist or now - hist[-1] >= T0) else 1     # chỉ dùng quá khứ: 'xa đợt' / 'vừa có bản ghi'
        kids = []
        for c in range(n_child):
            r2 = random.Random(seed * 100003 + s * 1009 + c)
            crec, _, _ = _evolve(psi.copy(), now, horizon, U, Ls, r2)
            full = hist + crec
            kids.append([phase_at(full, now + d) for d in DELTAS])
        res.append((ph_now, kids))
    return res


DELTAS = [0.04 * T0, 0.1 * T0, 0.25 * T0, 0.5 * T0, 1 * T0, 2 * T0, 4 * T0, 8 * T0]


def part_b(nproc=8, n_start=6, n_child=60):
    horizon = DELTAS[-1] + 1.5 * T0
    with Pool(nproc) as p:
        out = p.map(_fork, [(4000 + i, n_start, n_child, horizon) for i in range(nproc)])
    rows = [r for o in out for r in o]
    names = ["Quá khứ (bản ghi gần nhất ≥ T0 trước)", "vừa có bản ghi (< T0): Hiện tại hay Tương lai?"]
    # nhãn thật của chính 'bây giờ' (Δ = 0⁺) do tương lai quyết định: tỉ lệ HT/TL ngay tại Δ nhỏ nhất
    print("=== B. Cả ba pha xung đột: tách nhánh từ cùng một trạng thái 'bây giờ' (MCWF Part II) ===")
    print(f"  {len(rows)} điểm 'bây giờ' × {n_child} tương lai mỗi điểm; pha theo bản ghi Carnot (như ba-pha-xung-dot)")
    for ph_now in range(2):
        sel = [k for p_, k in rows if p_ == ph_now]
        if not sel:
            continue
        print(f"\n  bây giờ ở {names[ph_now]} ({len(sel)} điểm):")
        print("    Δ/T0     P(QK)   P(HT)   P(TL)   entropy pha H(Δ) (nat)   tỉ lệ điểm có tương lai KHÔNG duy nhất")
        for i, d in enumerate(DELTAS):
            P = np.zeros(3); Hs = []; nondet = 0
            for kids in sel:
                lab = np.array([kk[i] for kk in kids])
                q = np.bincount(lab, minlength=3) / len(lab)
                P += q; Hs.append(-(q[q > 0] * np.log(q[q > 0])).sum()); nondet += (q.max() < 1)
            P /= len(sel)
            print(f"    {d / T0:5.2f}    {P[0]:.3f}   {P[1]:.3f}   {P[2]:.3f}   {np.mean(Hs):.3f}"
                  f"                     {nondet / len(sel):.2f}")
    print()


# ======================================================================= C, D
def _setup():
    from giai_so_wdw_pw_e7 import channels, TAU
    from giai_so_kim_chi_gio import kraus, psi0
    ls, H, (Lc, Lb) = channels()
    return kraus(H, [Lc, Lb], TAU), psi0(), T0 / TAU            # Ks = [K·, K_C, K_B]


def part_c(M=300, nk=1500, late=500):
    """Người quan sát chỉ thấy mỗi bản ghi với xác suất p (bản ghi bị bỏ lỡ trộn vào 'không có gì').
    Bộ lọc của người quan sát: ρ ← K·ρK·† + (1−p)Σ_j K_jρK_j† khi không thấy gì; ρ ← K_jρK_j† khi thấy j.
    Độ mở với người quan sát: entropy của kết quả QUAN SÁT kế tiếp theo dự đoán của bộ lọc, nat/T0."""
    Ks, psi0, per = _setup()
    rng = np.random.default_rng(9)
    print("=== C. Người quan sát đọc một phần bản ghi (mỗi bản ghi đến tay với xác suất p) ===")
    print("  p      h_quan_sát (nat/T0)   ε* = 1/h (T0)   nhầm: −ln q(kết quả thật) (nat/T0)   ⟨S(ρ_người quan sát)⟩ (nat)")
    for pv in (1.0, 0.9, 0.7, 0.5, 0.3, 0.1, 0.0):
        hs = np.zeros(M); xs = np.zeros(M); Ss = np.zeros(M)
        for mtraj in range(M):
            psi = psi0.astype(complex); rho = np.outer(psi, psi.conj())
            for k in range(nk):
                outs = [K @ psi for K in Ks]
                pt = np.array([np.vdot(o, o).real for o in outs]); pt /= pt.sum()
                j = int(rng.choice(3, p=pt))
                seen = j if (j > 0 and rng.random() < pv) else 0
                Kr = [K @ rho @ K.conj().T for K in Ks]
                q = np.array([np.trace(x).real for x in Kr])
                qobs = np.array([q[0] + (1 - pv) * (q[1] + q[2]), pv * q[1], pv * q[2]])
                qobs = np.clip(qobs, 1e-300, None); qobs /= qobs.sum()
                if k >= nk - late:
                    hs[mtraj] += -(qobs * np.log(qobs)).sum()
                    xs[mtraj] += -math.log(qobs[seen])
                    ev = np.linalg.eigvalsh(rho); ev = ev[ev > 1e-12]
                    Ss[mtraj] += -(ev * np.log(ev)).sum()
                psi = outs[j] / math.sqrt(np.vdot(outs[j], outs[j]).real)
                rho = (Kr[0] + (1 - pv) * (Kr[1] + Kr[2])) if seen == 0 else Kr[seen]
                rho = rho / np.trace(rho).real
            if mtraj == 59 and M > 60:
                break
        n = mtraj + 1
        h = hs[:n].mean() / late * per; x = xs[:n].mean() / late * per; S = Ss[:n].mean() / late
        print(f"  {pv:.1f}    {h:8.3f}              {1 / h if h > 0 else float('inf'):7.3f}         {x:8.3f}"
              f"                          {S:.3f}", flush=True)
    print()


def part_d(M=300, nk=1500, late=500):
    """Tách h = ⟨−Σ p_j ln p_j⟩ theo kết quả j ∈ {· không có gì, C Carnot, B bơm}."""
    Ks, psi0, per = _setup()
    rng = np.random.default_rng(13)
    psis = np.tile(psi0.astype(complex), (M, 1))
    acc = np.zeros(3); rates = np.zeros(3); cnt = 0
    for k in range(nk):
        outs = np.stack([psis @ K.T for K in Ks], 0)
        p = np.real(np.einsum('jmd,jmd->jm', outs.conj(), outs)); p = np.clip(p, 0, None); p /= p.sum(0, keepdims=True)
        if k >= nk - late:
            acc += (-(p * np.log(np.where(p > 0, p, 1)))).mean(1); rates += p.mean(1); cnt += 1
        u = rng.random(M); c = np.cumsum(p, 0); j = (u[None, :] > c).sum(0).clip(0, 2)
        psis = outs[j, np.arange(M)] / np.sqrt(p[j, np.arange(M)])[:, None]
    acc = acc / cnt * per; rates = rates / cnt * per
    tot = acc.sum()
    print("=== D. Tách độ mở tương lai h theo loại bản ghi (E_dual, Part III) ===")
    names = ["· (không có gì: 'chưa xảy ra' cũng là thông tin)", "C Carnot (Hiện tại ↔ Tương lai)", "B bơm (Quá khứ ↔ Hiện tại)"]
    for nm, a, r in zip(names, acc, rates):
        print(f"  {nm:52s}: {a:.3f} nat/T0 ({100 * a / tot:4.1f}%), tần suất {r:.3f} /T0")
    print(f"  tổng h = {tot:.3f} nat/T0;  nat mỗi bản ghi Carnot ≈ {acc[1] / rates[1]:.2f}, mỗi bản ghi bơm ≈ {acc[2] / rates[2]:.2f}\n")


if __name__ == '__main__':
    for p in (sys.argv[1:] or ['A', 'B']):
        {'A': part_a, 'B': part_b, 'C': part_c, 'D': part_d}[p.upper()]()
