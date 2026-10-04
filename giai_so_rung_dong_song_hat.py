"""Rung động của khối lượng m sinh ra các mức năng lượng khác nhau ⟹ lúc thấy sóng, lúc thấy hạt, lúc đồng thời.
Học thuyết E_dual, HocThuyetE (7).pdf, Bảng I (ħω0 = 10 MeV, κ = ℓ0/3, k_BT = 1 MeV, γx = 0.70).

Thước đo lưỡng tính (Englert, D² + V² ≤ 1) cho dao động tử của m:
  V (sóng) = |⟨ψ(t)|ψ(t + T0/2)⟩|  — trạng thái sau nửa chu kỳ có giống chính nó (sóng dừng) không;
  D (hạt)  = khoảng cách vết giữa ψ(t) và ψ(t + T0/2) — phân biệt được hạt đang ở bên này hay bên kia không.
(A) Giải số trong cơ sở Fock (n ≤ 60): trạng thái kết hợp |α⟩ (năng lượng rung E_vib = |α|² ħω0, trải trên nhiều mức)
    và trạng thái riêng |n⟩ (một mức duy nhất) → V, D theo năng lượng; ngưỡng E_vib cho V = 0.9, D = V, D = 0.9.
(B) Quỹ đạo Langevin Bảng I (dùng lại run() của giai_so_dong_co_langevin.py qua một bản sao trong bộ nhớ ghi thêm
    (x, v, khóa?) — file gốc không bị sửa). Hạt bị lái theo quỹ đạo cổ điển ⟹ trạng thái lượng tử là trạng thái kết hợp
    α = (x + i v)/√2 (định lý dao động tử cưỡng bức); n_nhiệt = 1/(e^{10} − 1) ≈ 5·10⁻⁵ bỏ qua.
    Chia thời gian: Quá khứ (khóa), Ngưỡng (±0.05 T0 quanh lần gãy), Hiện tại (tự do), Tương lai (T0 đầu sau tái khóa).
Chạy: python3 giai_so_rung_dong_song_hat.py [A|B]   (A tức thì; B ≈ 1 phút, 8 lõi)
"""
import inspect
import math
import sys
from multiprocessing import Pool

import numpy as np

T0 = 2 * math.pi
NF = 60


def fock_ops():
    a = np.diag(np.sqrt(np.arange(1, NF + 1)), 1).astype(complex)
    return a, np.diag(np.arange(NF + 1) + 0.5)


def coherent(alpha):
    n = np.arange(NF + 1)
    lf = np.array([math.lgamma(k + 1) for k in n])
    c = np.exp(-abs(alpha)**2 / 2 + n * np.log(abs(alpha) + 1e-300) - lf / 2) * np.exp(1j * n * np.angle(alpha))
    if alpha == 0:
        c = np.zeros(NF + 1, complex); c[0] = 1
    return c / np.linalg.norm(c)


def vd(psi, H):
    """V, D giữa ψ và ψ sau nửa chu kỳ (tiến hóa số bằng e^{−iHπ})."""
    phi = np.exp(-1j * np.diag(H) * math.pi) * psi
    V = abs(np.vdot(psi, phi))
    rho = np.outer(psi, psi.conj()) - np.outer(phi, phi.conj())
    D = 0.5 * np.sum(np.abs(np.linalg.eigvalsh(rho)))
    return V, D


def part_a():
    a, H = fock_ops()
    print("=== A. Lưỡng tính theo năng lượng rung của m (cơ sở Fock, giải số) ===")
    print(f"{'E_vib (ħω0)':>11} {'E_vib (MeV)':>11} {'biên độ/ℓ0':>10} {'V kết hợp':>10} {'D kết hợp':>10} {'D²+V²':>7} "
          f"{'ΔE (MeV)':>9} {'V |n⟩':>6} {'D |n⟩':>6} {'đọc':>10}")
    for E in (0.0, 0.053, 0.173, 0.22, 0.38, 0.415, 1.0, 2.05, 3.0):
        al = math.sqrt(E)
        psi = coherent(al)
        V, D = vd(psi, H)
        dE = math.sqrt(max(0.0, np.real(np.vdot(psi, H @ H @ psi)) - np.real(np.vdot(psi, H @ psi))**2))
        n = int(round(E))
        e_n = np.zeros(NF + 1, complex); e_n[n] = 1
        Vn, Dn = vd(e_n, H)
        tag = "sóng" if V >= 0.9 else "hạt" if D >= 0.9 else "đồng thời"
        print(f"{E:11.3f} {E*10:11.2f} {2*al/math.sqrt(2):10.3f} {V:10.4f} {D:10.4f} {V*V+D*D:7.4f} {dE*10:9.2f} "
              f"{Vn:6.3f} {Dn:6.3f} {tag:>10}")
    print("giải tích: V = e^{−2E_vib}, D = √(1 − e^{−4E_vib}); V = 0.9 ⟺ E_vib = 0.053 ħω0 (0.53 MeV); "
          "D = V ⟺ ln2/4 = 0.173 (1.73 MeV); D = 0.9 ⟺ 0.415 (4.15 MeV)")
    print("trạng thái riêng |n⟩ (năng lượng xác định): V = 1, D = 0 ở MỌI mức ⟹ sóng dừng thuần\n")


def _run_full_factory():
    import giai_so_dong_co_langevin as L
    src = inspect.getsource(L.run).replace("def run(", "def run_full(").replace(
        "samples.append(s[0])", "samples.append((t, s[0], s[1], s[2], s[3], locked))")
    exec(src, L.__dict__)
    return L


def _traj(a):
    from giai_so_kramers_A import setup
    cfg, seed, T, every = a
    L = _run_full_factory()
    setup(cfg)
    wc = cfg[2]
    r = L.run_full(0.1, 1.05 * wc, 0.1, T, 2e-3, seed, gx=0.7, w_start=0.82 * wc, sample_every=every)
    return np.array(r['samples'], float), r['events']


def part_b(N=8):
    from giai_so_burst_g import cfg_for
    cfg = cfg_for(0.1)
    every = 0.01
    with Pool() as pool:
        res = pool.map(_traj, [(cfg, 95000 + i, 60 * T0, every) for i in range(N)])
    groups = {k: [] for k in ("Quá khứ", "Ngưỡng", "Hiện tại", "Tương lai")}
    for S, ev in res:
        t, x, v, lk = S[:, 0], S[:, 3], S[:, 4], S[:, 5] > 0.5
        keep = t > 5 * T0
        brk = np.array([e[1] for e in ev if e[0] == 'gãy'])
        lock = np.array([e[1] for e in ev if e[0] == 'khóa'])
        near_b = np.zeros_like(keep)
        for tb in brk:
            near_b |= np.abs(t - tb) < 0.05 * T0
        after_l = np.zeros_like(keep)
        for tl in lock:
            after_l |= (t > tl) & (t < tl + T0)
        E = 0.5 * (x**2 + v**2)
        for name, msk in (("Ngưỡng", near_b), ("Hiện tại", ~lk & ~near_b), ("Tương lai", lk & after_l & ~near_b),
                          ("Quá khứ", lk & ~after_l & ~near_b)):
            groups[name].append(E[keep & msk])
    print(f"=== B. Quỹ đạo Bảng I (1 MeV), {N} × 60 T0: năng lượng rung của m và lưỡng tính theo pha ===")
    print(f"{'pha':>10} {'tỉ lệ t/g':>9} {'E_vib TB (MeV)':>15} {'max':>6} {'V TB':>6} {'D TB':>6} "
          f"{'% thời gian sóng (V≥0.9)':>25} {'đồng thời':>10} {'hạt (D≥0.9)':>12}")
    tot = sum(sum(len(e) for e in g) for g in groups.values())
    for name in ("Quá khứ", "Ngưỡng", "Hiện tại", "Tương lai"):
        E = np.concatenate(groups[name])
        V = np.exp(-2 * E); D = np.sqrt(1 - np.exp(-4 * E))
        w, p = np.mean(V >= 0.9), np.mean(D >= 0.9)
        print(f"{name:>10} {len(E)/tot:9.1%} {E.mean()*10:15.2f} {E.max()*10:6.1f} {V.mean():6.3f} {D.mean():6.3f} "
              f"{w:25.1%} {1-w-p:10.1%} {p:12.1%}")
    print("(E_vib = ½(x² + v²) ħω0 = năng lượng rung trên điểm không; V = e^{−2E_vib}, D = √(1 − V²))\n")


if __name__ == "__main__":
    for name in sys.argv[1:] or ["A", "B"]:
        {"A": part_a, "B": part_b}[name.upper()]()
