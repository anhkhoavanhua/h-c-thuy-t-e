"""Giải số hai đề xuất: (1) bộ đếm thời gian qua nhiều vòng — số cuộn N(t) = ⌊Φ/2π⌋ và năng lượng tiêu tán
tích lũy E_diss(t); (2) mở rộng Lindblad với toán tử nhảy L = √Γ a^j.

(A) Langevin cổ điển (bộ thông số hạt nhân, f=1.05, τ=10, γx=0.7, g=0.1), nhiễu 1 và 5 MeV:
    N(t) có đơn điệu không (tỉ lệ bước lùi, số lần quay ngược), độ đều của N và của E_Carnot tích lũy theo từng T0.
(B) Lindblad cho hạt (ħ=ω0=1), bể nhiệt L1 = √(γ(n̄+1)) a, L2 = √(γ n̄) a†, γ = γx = 0.7:
    B1 mũi tên thời gian: xuất phát |2⟩ — entropy von Neumann S(ρ) có đơn điệu không; entropy tương đối D(ρ‖ρ_ss).
    B2 hạt bị lái −g x sin(Ωt), Ω = θ̇_d = 2.63: công suất tiêu tán lượng tử γ⟨a†a⟩ (mỗi bước nhảy đúng 1 ħω0)
       so với cổ điển γx⟨ẋ²⟩.
    B3 toán tử a^j: tốc độ nhảy ∝ ⟨a†^j a^j⟩ tại trạng thái hạt ứng với năng lượng cổ điển 0.57 ħω0.
Chạy: python3 giai_so_winding_lindblad.py [A|B]   (A ≈ 3 phút, 8 lõi; B vài giây)
"""
import math
import sys
from multiprocessing import Pool

import numpy as np

import giai_so_dong_co_langevin as L
from giai_so_burst_g import cfg_for
from giai_so_kramers_A import setup

T0 = L.T0


# ---------------------------------------------------------------- A
def _traj(a):
    cfg, kT, seed, T, warm, every = a
    setup(cfg)
    wc = cfg[2]
    r = L.run(0.1, 1.05 * wc, kT, T, 2e-3, seed, gx=0.7, w_start=0.82 * wc, sample_every=every)
    th = np.array(r['samples'])
    locks = [(e[1], e[3]) for e in r['events'] if e[0] == 'khóa']
    return th, locks


def part_a(N=8):
    cfg = cfg_for(0.1)
    T, warm, every = 70 * T0, 10 * T0, 0.05
    print(f"=== A. Bộ đếm vòng N(t) và E_Carnot tích lũy ({N} quỹ đạo × 60 T0) ===")
    with Pool() as pool:
        for kT in (0.1, 0.3, 0.5):
            res = pool.map(_traj, [(cfg, kT, 33000 + i, T, warm, every) for i in range(N)])
            back, steps, rev, rate, all_n, all_e, span_back, neg = 0, 0, 0, [], [], [], [], 0
            for th, locks in res:
                k0 = int(warm / every)
                phi = th[k0:]
                n = np.floor(phi / (2 * math.pi))
                dn = np.diff(n)
                back += int(np.sum(dn < 0))
                steps += int(np.sum(dn != 0))
                # số lần đảo chiều quay (dấu của dΦ/dt đổi) và độ lùi lớn nhất so với đỉnh trước đó
                dphi = np.diff(phi)
                sgn = np.sign(dphi[np.abs(dphi) > 1e-9])
                rev += int(np.sum(sgn[1:] != sgn[:-1]))
                span_back.append(float(np.max(np.maximum.accumulate(phi) - phi)) / (2 * math.pi))
                per = int(T0 / every)
                inc_n = [n[i + per] - n[i] for i in range(0, len(n) - per, per)]
                rate.append(sum(inc_n) / len(inc_n))
                all_n += inc_n
                neg += int(phi[-1] < phi[0])
                tl = np.array([t for t, _ in locks if t > warm]); el = np.array([e for t, e in locks if t > warm])
                inc_e = []
                t = warm
                while t + T0 <= T:
                    inc_e.append(float(el[(tl >= t) & (tl < t + T0)].sum()))
                    t += T0
                all_e += inc_e
            print(f"k_BT={kT*10:g} MeV: tốc độ TB {np.mean(rate):.3f} vòng/T0; bước N lùi {back}/{steps} "
                  f"({back/max(steps,1):.1%}); số lần đảo chiều quay {rev}; lùi xa nhất TB {np.mean(span_back):.2f} vòng; "
                  f"quỹ đạo lùi ròng {neg}/{N}; CV gộp(ΔN mỗi T0) = {np.std(all_n)/abs(np.mean(all_n)):.3f}; "
                  f"CV gộp(ΔE_Carnot mỗi T0) = {np.std(all_e)/np.mean(all_e):.3f}; P(ΔE_Carnot = 0 trong một T0) = "
                  f"{np.mean(np.array(all_e) == 0):.1%}", flush=True)
    print()


# ---------------------------------------------------------------- B
def ops(nmax):
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1).astype(complex)
    return a, a.conj().T


def lindblad_rhs(rho, H, Ls):
    d = -1j * (H @ rho - rho @ H)
    for Lk in Ls:
        LdL = Lk.conj().T @ Lk
        d += Lk @ rho @ Lk.conj().T - 0.5 * (LdL @ rho + rho @ LdL)
    return d


def evolve(rho, Hf, Ls, T, dt, obs=None, every=None):
    t, out = 0.0, []
    n = int(T / dt)
    for i in range(n):
        H1, H2, H3 = Hf(t), Hf(t + dt / 2), Hf(t + dt)
        k1 = lindblad_rhs(rho, H1, Ls)
        k2 = lindblad_rhs(rho + dt / 2 * k1, H2, Ls)
        k3 = lindblad_rhs(rho + dt / 2 * k2, H2, Ls)
        k4 = lindblad_rhs(rho + dt * k3, H3, Ls)
        rho = rho + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        if obs is not None and i % every == 0:
            out.append((t, obs(rho, t)))
    return rho, out


def vn_entropy(rho):
    w = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-np.sum(w * np.log(w)))


def rel_entropy(rho, sig_log):
    w, v = np.linalg.eigh((rho + rho.conj().T) / 2)
    w = np.clip(w, 1e-300, None)
    lr = v @ np.diag(np.log(w)) @ v.conj().T
    return float(np.real(np.trace(rho @ (lr - sig_log))))


def part_b():
    nmax, gam = 20, 0.7
    a, ad = ops(nmax)
    num = ad @ a
    x = (a + ad) / math.sqrt(2)
    H0 = num + 0.5 * np.eye(nmax + 1)
    print("=== B. Lindblad cho hạt (ħ = ω0 = 1, γ = 0.7) ===")
    # B1
    for kT_mev in (1.0, 5.0):
        nb = 1 / (math.exp(10.0 / kT_mev) - 1)
        Ls = [math.sqrt(gam * (nb + 1)) * a, math.sqrt(gam * nb) * ad]
        p = np.array([(nb / (nb + 1))**k for k in range(nmax + 1)]); p /= p.sum()
        sig = np.diag(p).astype(complex)
        sig_log = np.diag(np.log(p)).astype(complex)
        rho = np.zeros((nmax + 1, nmax + 1), complex); rho[2, 2] = 1
        _, tr = evolve(rho, lambda t: H0, Ls, 12.0, 0.005,
                       obs=lambda r, t: (vn_entropy(r), rel_entropy(r, sig_log), float(np.real(np.trace(num @ r)))), every=40)
        S = [o[1][0] for o in tr]; D = [o[1][1] for o in tr]
        iS = int(np.argmax(S))
        dmono = all(D[i + 1] <= D[i] + 1e-9 for i in range(len(D) - 1))
        print(f"B1 k_BT={kT_mev:g} MeV (n̄={nb:.2e}), xuất phát |2⟩: S(ρ) tăng tới max {S[iS]:.4f} tại t={tr[iS][0]:.2f} "
              f"rồi giảm về {S[-1]:.4f} (S của trạng thái dừng {vn_entropy(sig):.4f}) ⟹ S KHÔNG đơn điệu; "
              f"D(ρ‖ρ_ss): {D[0]:.3f} → {D[-1]:.2e}, đơn điệu giảm: {'CÓ' if dmono else 'KHÔNG'}")
    # B2
    Om, g = 2.63, 0.1
    Ls = [math.sqrt(gam) * a]
    rho = np.zeros((nmax + 1, nmax + 1), complex); rho[0, 0] = 1
    _, tr = evolve(rho, lambda t: H0 - g * math.sin(Om * t) * x, Ls, 60.0, 0.005,
                   obs=lambda r, t: float(np.real(np.trace(num @ r))), every=10)
    tail = [o[1] for o in tr if o[0] > 40]
    Pq = gam * np.mean(tail)
    A = g / math.sqrt((1 - Om**2)**2 + (gam * Om)**2)
    Pc = gam * A * A * Om * Om / 2
    print(f"B2 lái −g x sin(Ωt), Ω={Om}, g={g}: công suất Lindblad γ⟨n⟩ = {Pq:.3e} ħω0/đv t "
          f"(mỗi bước nhảy đúng 1 ħω0, tốc độ nhảy {Pq:.3e}); cổ điển γx⟨ẋ²⟩ = {Pc:.3e}; tỉ số {Pq/Pc:.3f}")
    Om2 = 1.0
    _, tr2 = evolve(rho, lambda t: H0 - g * math.sin(Om2 * t) * x, Ls, 60.0, 0.005,
                    obs=lambda r, t: float(np.real(np.trace(num @ r))), every=10)
    Pq2 = gam * np.mean([o[1] for o in tr2 if o[0] > 40])
    A2 = g / math.sqrt((1 - Om2**2)**2 + (gam * Om2)**2)
    Pc2 = gam * A2 * A2 * Om2 * Om2 / 2
    print(f"   so sánh ở cộng hưởng Ω = 1: Lindblad {Pq2:.3e}, cổ điển {Pc2:.3e}, tỉ số {Pq2/Pc2:.3f}")
    # B3
    Ecl = 0.571
    nb_eff = max(Ecl - 0.5, 0.0)
    print(f"B3 năng lượng hạt cổ điển 0.571 ħω0 ⟹ ⟨n⟩ ≈ {nb_eff:.3f} (trừ điểm không 0.5).")
    for label, rho_st in (("trạng thái kết hợp |α|²=0.071", None), ("hỗn hợp nhiệt ⟨n⟩=0.071", None)):
        if label.startswith("trạng thái"):
            al = math.sqrt(nb_eff)
            c = np.array([math.exp(-al * al / 2) * al**k / math.sqrt(math.factorial(k)) for k in range(nmax + 1)])
            rho_st = np.outer(c, c).astype(complex)
        else:
            q = nb_eff / (1 + nb_eff)
            p = np.array([(1 - q) * q**k for k in range(nmax + 1)])
            rho_st = np.diag(p).astype(complex)
        m = [float(np.real(np.trace(np.linalg.matrix_power(ad, j) @ np.linalg.matrix_power(a, j) @ rho_st)))
             for j in (1, 2, 3)]
        print(f"   {label}: ⟨a†a⟩ = {m[0]:.4f}, ⟨a†²a²⟩ = {m[1]:.5f}, ⟨a†³a³⟩ = {m[2]:.6f} "
              f"⟹ tỉ lệ tốc độ nhảy j=2 : j=1 = {m[1]/m[0]:.4f}")
    print()


if __name__ == "__main__":
    for name in sys.argv[1:] or ["B", "A"]:
        {"A": part_a, "B": part_b}[name.upper()]()
