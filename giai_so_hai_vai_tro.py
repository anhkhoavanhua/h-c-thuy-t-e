"""Hai vai trò trong cơ chế Page–Wootters của E_dual (HocThuyetE (7).pdf, Part III; nối tiếp giai_so_kim_chi_gio.py).

Vai trò 1: kim lý tưởng K̂ = Σ k|k⟩⟨k| (chuỗi xung τ = T0/25) tạo ra thời gian t = kτ và phương trình Lindblad.
  1a  Lindblad nổi lên O(τ) với MỌI mô hình bể (Carnot + bơm; chỉ Carnot) — không phụ thuộc bể.
  1b  Không cần đọc bản ghi: trung bình các nhánh bản ghi (quỹ đạo lượng tử) = Λ^k ρ0; từng nhánh thì thuần và lệch xa.
Vai trò 2: bể Carnot tạo mũi tên qua entropy bản ghi thô, không qua S(E) toàn cục.
  2a  Từng nhánh: surprisal s_k = −ln p(r_1…r_k) của chuỗi bản ghi (cơ sở con trỏ {·, C, B}) không bao giờ giảm;
      trung bình ⟨s_k⟩ = H(r_1…r_k) (Shannon) tăng tuyến tính với tốc độ h; so S(E) = S(ρ_S) bão hòa, Σ_i S(E_i),
      và H(N_C) của riêng bộ đếm Carnot (chỉ tăng ~ log k).
  2b  Mũi tên tắt khi tiêu tán tắt: chỉ Carnot ở T = 0, hệ về trạng thái cơ bản ⟹ h → 0.
Phép thử chéo: mũi tên của bể có cùng chiều với k? Hai nhánh độc lập: P(s_{k+Δk} > s_k) theo Δk; Δk cần cho 99%.
Chạy: python3 giai_so_hai_vai_tro.py
"""
import math
import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

from giai_so_wdw_pw_e7 import TAU, T0, channels, vn
from giai_so_kim_chi_gio import kraus, psi0, lind_at, apply, tdist, marginal_entropy, sandwich, _cmul

M_TRAJ, NK = 4000, 1500                      # số nhánh bản ghi, số tích tắc (60 T0)


def branches(Ks, p0, nk, M, rng, keep=(25, 125)):
    """Lấy mẫu M nhánh bản ghi: ψ → K_j ψ/√p_j, s += −ln p_j. Trả s[k] (nk+1, M), số Carnot, ρ trung bình tại keep."""
    d = len(p0)
    Psi = np.tile(p0, (M, 1)).astype(complex)
    s = np.zeros((nk + 1, M)); nc = np.zeros(M, int); rho_avg = {}
    KT = [np.ascontiguousarray(K.T) for K in Ks]
    for k in range(1, nk + 1):
        Phi = [_cmul(Psi, KT[j]) for j in range(len(Ks))]
        P = np.stack([np.sum(np.abs(F)**2, axis=1) for F in Phi], axis=1)
        P /= P.sum(axis=1, keepdims=True)
        u = rng.random(M)[:, None]
        j = np.minimum((u > np.cumsum(P, axis=1)).sum(axis=1), len(Ks) - 1)
        pj = P[np.arange(M), j]
        Psi = np.choose(j[:, None], Phi) / np.sqrt(pj)[:, None]
        s[k] = s[k - 1] - np.log(pj)
        nc += (j == 1)
        if k in keep:
            rho_avg[k] = _cmul(Psi.T, Psi.conj()) / M
    return s, nc, rho_avg, Psi


def counting_entropy(Ks, rho0, nk, ks, nmax=60):
    """H(N_C) chính xác tại các k: đệ quy ρ_n(k)."""
    d = rho0.shape[0]
    rn = np.zeros((nmax + 1, d, d), complex); rn[0] = rho0
    out = {}
    for k in range(1, nk + 1):
        new = sandwich(Ks[0], rn) + sandwich(Ks[2], rn)
        new[1:] += sandwich(Ks[1], rn[:-1])
        new[np.abs(new) < 1e-200] = 0
        rn = new
        if k in ks:
            p = np.real(np.einsum('naa->n', rn)); p = p[p > 1e-300]
            out[k] = float(-np.sum(p * np.log(p)))
    return out


def role1(H, Ls, rho0, rng):
    print("=== Vai trò 1. Kim lý tưởng K̂: thời gian t = kτ và phương trình Lindblad ===")
    print("1a. Lindblad nổi lên với mọi mô hình bể (sai số vết tại t = 1 T0; tỉ số khi chia đôi τ):")
    for lab, LL in (("Carnot + bơm (E_dual)", Ls), ("chỉ Carnot (T = 0)", [Ls[0]]), ("chỉ bơm", [Ls[1]])):
        ref = lind_at(H, LL, rho0, [T0])[0]
        errs = []
        for m in range(4):
            t = TAU / 2**m; Kx = kraus(H, LL, t); r = rho0.copy()
            for _ in range(int(round(T0 / t))):
                r = apply(Kx, r)
            errs.append(tdist(r, ref))
        print(f"   {lab:24s}: " + ", ".join(f"{e:.2e}" for e in errs)
              + "   tỉ số " + ", ".join(f"{errs[i]/errs[i+1]:.2f}" for i in range(3)))
    print("   ⟹ K̂ tái lập Lindblad O(τ) bất kể bể là gì: thời gian và phương trình chuyển động đến từ K̂, không từ bể\n")

    print(f"1b. Không cần đọc bản ghi ({M_TRAJ} nhánh, τ = T0/25):")
    Ks = kraus(H, Ls, TAU)
    s, nc, ravg, _ = branches(Ks, psi0(), 125, M_TRAJ, rng)
    for k in (25, 125):
        r = rho0.copy()
        for _ in range(k):
            r = apply(Ks, r)
        print(f"   k = {k:3d} ({k*TAU/T0:.0f} T0): ‖trung bình các nhánh − Λ^k ρ0‖ = {tdist(ravg[k], r):.3f} "
              f"(nhiễu thống kê ~1/√M = {1/math.sqrt(M_TRAJ):.3f}); mỗi nhánh thuần (Tr ρ² = 1), "
              f"Tr(Λ^k ρ0)² = {np.real(np.trace(r @ r)):.3f}")
    print("   ⟹ ρ_S(t) do K̂ quy định là trung bình trên MỌI bản ghi: người đọc kim không cần đọc bể\n")


def role2(H, Ls, rho0, rng):
    print("=== Vai trò 2. Bể Carnot: mũi tên qua entropy bản ghi thô ===")
    Ks = kraus(H, Ls, TAU)
    s, nc, _, _ = branches(Ks, psi0(), NK, M_TRAJ, rng)
    ds = np.diff(s, axis=0)
    print(f"2a. {M_TRAJ} nhánh × {NK} tích tắc: số bước s_k GIẢM trên một nhánh = {int(np.sum(ds < 0))} / {ds.size}; "
          f"bước tăng nhỏ nhất {ds.min():.2e} nat")
    ks = (1, 7, 25, 125, 250, 500, 1000, 1500)
    HN = counting_entropy(Ks, rho0, NK, set(ks))
    r = rho0.copy(); SE = {0: 0.0}; Sm = {0: 0.0}; acc = 0.0
    for k in range(1, NK + 1):
        acc += marginal_entropy(Ks, r); r = apply(Ks, r)
        if k in ks:
            SE[k] = vn(r); Sm[k] = acc
    print(f"{'k':>5} {'t/T0':>5} {'S(E)=S(ρ_S)':>12} {'Σ S(E_i)':>9} {'H(chuỗi) = ⟨s_k⟩':>18} {'sd(s_k)':>8} {'H(N_C)':>7}")
    for k in ks:
        print(f"{k:5d} {k*TAU/T0:5.2f} {SE[k]:12.3f} {Sm[k]:9.2f} {s[k].mean():12.2f} ± {s[k].std()/math.sqrt(M_TRAJ):.2f}"
              f" {s[k].std():8.2f} {HN[k]:7.2f}")
    h = (s[NK].mean() - s[NK - 500].mean()) / (500 * TAU / T0)
    D = (s[NK].var() - s[NK - 500].var()) / (500 * TAU / T0)
    print(f"   tốc độ entropy bản ghi h = {h:.2f} nat/T0 (Σ S(E_i): {(Sm[1500]-Sm[1000])/20:.2f}); khuếch tán D = {D:.2f} nat²/T0")
    print(f"   H(N_C) 1000→1500: +{HN[1500]-HN[1000]:.2f} nat (~½ ln, gần bão hòa) — bộ đếm Carnot riêng KHÔNG mang tốc độ entropy"
          " tuyến tính; mũi tên nằm ở toàn chuỗi bản ghi (thời điểm + loại)\n")

    print("2b. Mũi tên tắt khi tiêu tán tắt (chỉ Carnot, T = 0, xuất phát trạng thái riêng thứ 6):")
    Hd, Vd = np.linalg.eigh(H)
    Kc = kraus(H, [Ls[0]], TAU)
    sc, _, _, _ = branches(Kc, Vd[:, 5], 2500, 2000, rng, keep=())
    r = np.outer(Vd[:, 5], Vd[:, 5].conj())
    for k in (25, 125, 500, 1000, 2500):
        rr = r.copy()
        for _ in range(k):
            rr = apply(Kc, rr)
        rate = (sc[k].mean() - sc[k - 25].mean())
        print(f"   k = {k:4d} ({k*TAU/T0:5.1f} T0): S(E) = {vn(rr):.3f}, H(chuỗi) = {sc[k].mean():.3f}, "
              f"tăng trong 1 T0 vừa qua = {rate:.4f} nat")
    print("   ⟹ hệ về trạng thái cơ bản ⟹ không còn bản ghi mới, h → 0: mũi tên chỉ có khi động cơ còn tiêu tán (bơm + Carnot)\n")
    return s, h, D


def cross(s, h, D):
    print("=== Phép thử chéo: mũi tên của bể có cùng chiều với kim K̂? ===")
    k0 = 500
    print(f"   cùng một nhánh: s_(k+Δk) ≥ s_k với mọi Δk ≥ 0 (vai trò 2a: 0 lần giảm) ⟹ không bao giờ ngược chiều K̂")
    print(f"   hai nhánh độc lập (hai người quan sát, mỗi người một bản ghi), từ k0 = {k0} (20 T0):")
    rng = np.random.default_rng(3)
    a = s[k0]
    need = None
    for dk in (1, 5, 25, 50, 125, 250, 300, 400, 600):
        b = rng.permutation(s[k0 + dk])
        p = float(np.mean(b > a))
        if need is None and p >= 0.99:
            need = dk
        print(f"      Δk = {dk:4d} ({dk*TAU/T0:5.2f} T0): P(s_b(k0+Δk) > s_a(k0)) = {p:.3f}")
    print("   độ phân giải đồng hồ entropy (một nhánh, đọc t̂ = s_k/h): σ_t = sd(s_k)/h = "
          + ", ".join(f"{s[k].std()/h:.2f} T0 (t = {k*TAU/T0:g} T0)" for k in (25, 125, 500, 1500))
          + "; cần 0.071 T0 cho 99%")
    print(f"   ⟹ mũi tên của bể luôn cùng chiều K̂ trên một nhánh; giữa hai nhánh cần Δk ≳ {need} tích tắc "
          f"({need*TAU/T0 if need else float('nan'):.1f} T0) để xếp đúng 99%: bể cho CHIỀU, K̂ cho SỐ ĐO\n")


def main():
    ls, H, Ls = channels()
    rho0 = np.outer(psi0(), psi0().conj())
    rng = np.random.default_rng(11)
    print(f"S {H.shape[0]} chiều; τ = T0/25; ancilla {{·, C (Carnot), B (bơm)}}\n")
    role1(H, Ls, rho0, rng)
    s, h, D = role2(H, Ls, rho0, rng)
    cross(s, h, D)


if __name__ == "__main__":
    main()
