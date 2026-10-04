"""Phương trình Schrödinger trên đa tạp ràng buộc (pha Quá khứ, Eq. (6)) của HocThuyetE (7).pdf, Bảng I.

Đơn vị ħ = m = ω0 = ℓ0 = 1; I = 1, k = 1, ΔE = 1, g = 0.1, κ = 1/3, R = 3, θ̇_br = 3.205.
Lagrangian rút gọn trên x = κ sinθ:  L = ½ M(θ) θ̇² − U(θ),
  M(θ) = I + mκ² cos²θ,  U(θ) = −ΔE cosθ + ½kκ² sin²θ − gκ sin²θ.
(S1) Lượng tử hóa: ba cách sắp toán tử
     LB  : Laplace–Beltrami trên metric M(θ)  (−½ M^{-1/2} ∂ M^{-1/2} ∂ trong độ đo √M dθ)
     dC  : LB + thế hình học da Costa V_g = −K²/8 (ràng buộc thực hiện bằng giam cứng, K = độ cong
           của đường cong trong tọa độ trọng-khối (√I θ, √m x))
     naive: −½ ∂ M^{-1} ∂ (độ đo dθ)
     Phổ, số mức dưới năng lượng gãy E_br = ½M(0)θ̇_br² + U(0), so với đếm Bohr–Sommerfeld.
(S2) Khoảng cách mức quanh ngưỡng so với ⟨ΔE_Carnot⟩ = 0.40–0.46 ħω0 (V.B.2) và 2.05 ħω0/bước (V.C.1).
(S3) TDSE: bó sóng rô-to tại điểm làm việc θ̇_d = 0.821 θ̇_br; độ kết hợp pha |⟨e^{iθ}⟩| theo t, so với tập hợp
     cổ điển cùng phân bố ban đầu (Liouville) và với thời gian lưu 4.5–6 T0 của bản thảo.
Chạy: python3 giai_so_schrodinger_rang_buoc.py   (≈ 15 s)
"""
import math

import numpy as np

I_, M_, K_, DE, G_, KAP = 1.0, 1.0, 1.0, 1.0, 0.1, 1 / 3
WBR = 3.2048
T0 = 2 * math.pi
N = 512
TH = 2 * math.pi * np.arange(N) / N
H_ = TH[1] - TH[0]


def Mf(th):
    return I_ + M_ * KAP**2 * np.cos(th)**2


def Uf(th):
    return -DE * np.cos(th) + (0.5 * K_ * KAP**2 - G_ * KAP) * np.sin(th)**2


def hamiltonian(kind):
    """Ma trận H trên lưới tuần hoàn, biểu diễn độ đo phẳng (φ = M^{1/4} ψ cho LB/dC)."""
    Dm = (np.roll(np.eye(N), 1, axis=1) - np.eye(N)) / H_          # sai phân tiến: (φ_{i+1} − φ_i)/h
    mid = TH + H_ / 2
    if kind == "naive":
        T = 0.5 * Dm.T @ np.diag(1 / Mf(mid)) @ Dm
        V = Uf(TH)
    else:
        W = np.diag(Mf(TH)**-0.25)
        T = 0.5 * W @ Dm.T @ np.diag(Mf(mid)**-0.5) @ Dm @ W
        V = Uf(TH).copy()
        if kind == "dC":
            Kc = math.sqrt(I_ * M_) * KAP * np.abs(np.sin(TH)) / Mf(TH)**1.5
            V = V - Kc**2 / 8
    return T + np.diag(V)


def bohr_sommerfeld(E):
    """Số trạng thái N(E) = ∮ p dθ /(2π) (quay) hoặc tương ứng cho dao động (libration)."""
    p = np.sqrt(np.maximum(0.0, 2 * Mf(TH) * (E - Uf(TH))))
    J = np.sum(p) * H_
    if np.all(E > Uf(TH)):
        return 2 * J / (2 * math.pi)        # hai chiều quay ±
    return J * 2 / (2 * math.pi)            # libration: ∮ = 2∫ p dθ trên vùng được phép


def part_s1():
    Ebr = 0.5 * Mf(0.0) * WBR**2 + Uf(0.0)
    print("=== S1. Phổ Schrödinger trên đa tạp ràng buộc (Bảng I) ===")
    print(f"U(0) = {Uf(0.0):.4f}, U(π) = {Uf(math.pi):.4f} (rào 2ΔE = {Uf(math.pi)-Uf(0.0):.3f});  "
          f"E_br = ½M(0)θ̇_br² + U(0) = {Ebr:.4f} ħω0 = {Ebr*10:.2f} MeV")
    res = {k: np.linalg.eigvalsh(hamiltonian(k)) for k in ("LB", "dC", "naive")}
    _, V = np.linalg.eigh(hamiltonian("LB"))
    kk = np.fft.fftfreq(N, d=1 / N)
    Mbar = 1 / np.mean(1 / Mf(TH))
    print(f"{'mức':>4} {'LB':>9} {'LB + da Costa':>14} {'naive':>9} {'√⟨p_θ²⟩':>8} {'ℓ²/2M̄':>7} {'loại':>10}")
    for i in range(12):
        e = res["LB"][i]
        f = np.abs(np.fft.fft(V[:, i]))**2
        prms = math.sqrt(np.sum(f * kk**2) / np.sum(f))
        kind = "dao động" if e < Uf(math.pi) else "quay"
        print(f"{i:4d} {e:9.4f} {res['dC'][i]:14.4f} {res['naive'][i]:9.4f} {prms:8.3f} "
              f"{round(prms)**2/(2*Mbar):7.3f} {kind:>10}")
    print(f"(M̄ = 1/⟨1/M⟩ = {Mbar:.4f}; ℓ = làm tròn √⟨p_θ²⟩)")
    for k, ev in res.items():
        n = int(np.sum(ev < Ebr))
        print(f"{k:>6}: {n} mức dưới E_br (Bohr–Sommerfeld N(E_br) = {bohr_sommerfeld(Ebr):.2f})")
    d = np.max(np.abs(res["dC"][:10] - res["LB"][:10]))
    d2 = np.max(np.abs(res["naive"][:10] - res["LB"][:10]))
    print(f"chênh lệch do cách sắp toán tử (10 mức đầu): da Costa ≤ {d:.4f}, naive ≤ {d2:.4f} ħω0\n")
    return res["LB"], Ebr


def part_s2(ev, Ebr):
    print("=== S2. Khoảng cách mức quanh ngưỡng gãy (LB) ===")
    i = int(np.sum(ev < Ebr))
    print("các mức quanh E_br (cặp gần suy biến ±ℓ của chuyển động quay):")
    for j in range(max(0, i - 6), i + 3):
        mark = "  ← trên ngưỡng" if ev[j] >= Ebr else ""
        print(f"  E_{j} = {ev[j]:8.4f}{mark}")
    # bước nhảy đi xuống nhỏ nhất (giữa các doublet khác nhau) từ doublet cao nhất dưới ngưỡng
    top = ev[i - 1]
    lower = ev[ev < top - 0.05]
    step = top - lower.max()
    print(f"mức quay cao nhất của pha Quá khứ: |ℓ| = 3 (E = {top:.3f}); mức đầu trên ngưỡng: |ℓ| = 4 (E = {ev[i]:.3f})")
    print(f"bước đi xuống nhỏ nhất ℓ 3→2: {step:.3f} ħω0 = {step*10:.1f} MeV;  bước ℓ 4→3: {ev[i]-top:.3f} ħω0")
    print(f"so với ⟨ΔE_Carnot⟩ cổ điển 0.40–0.46 ħω0 (V.B.2): bước lượng tử lớn hơn {step/0.46:.1f}–{step/0.40:.1f} lần")
    print(f"so với MCWF V.C.1: ≈ 2.05 ħω0 mỗi bước nhảy\n")


def classical_traj(th0, w0, T, dt=2e-3):
    """Eq. (6) không cản, không lái: RK4 theo (θ, θ̇), vector hóa cho cả tập hợp."""
    def acc(th, w):
        Mm = Mf(th)
        dM = -2 * M_ * KAP**2 * np.cos(th) * np.sin(th)
        dU = DE * np.sin(th) + (K_ * KAP**2 - 2 * G_ * KAP) * np.sin(th) * np.cos(th)
        return (-0.5 * dM * w * w - dU) / Mm
    th, w = np.array(th0, float), np.array(w0, float)
    out, t = [], 0.0
    nst = int(round(T / dt))
    rec = max(1, nst // 400)
    for s in range(nst + 1):
        if s % rec == 0:
            out.append((t, np.abs(np.mean(np.exp(1j * th))), th.copy()))
        k1t, k1w = w, acc(th, w)
        k2t, k2w = w + dt / 2 * k1w, acc(th + dt / 2 * k1t, w + dt / 2 * k1w)
        k3t, k3w = w + dt / 2 * k2w, acc(th + dt / 2 * k2t, w + dt / 2 * k2w)
        k4t, k4w = w + dt * k3w, acc(th + dt * k3t, w + dt * k3w)
        th = th + dt / 6 * (k1t + 2 * k2t + 2 * k3t + k4t)
        w = w + dt / 6 * (k1w + 2 * k2w + 2 * k3w + k4w)
        t += dt
    return out


def part_s3():
    wd = 0.821 * WBR
    p0 = Mf(0.0) * wd
    Hm = hamiltonian("LB")
    E, V = np.linalg.eigh(Hm)
    Tmax = 6 * T0
    print(f"=== S3. TDSE bó sóng rô-to tại θ̇_d = 0.821 θ̇_br = {wd:.3f} (p_θ = {p0:.3f} ħ) ===")
    print(f"{'Δp_θ (ħ)':>9} {'Δθ (rad)':>9} {'t mất pha LT (T0)':>18} {'t mất pha CĐ (T0)':>18} "
          f"{'|⟨e^iθ⟩| LT @ 4.5T0':>20} {'CĐ @ 4.5T0':>11} {'hồi sinh LT quanh T_rev (t, T0)':>32}")
    rng = np.random.default_rng(1)
    for dp in (0.5, 1.0, 2.0):
        sig = 1 / (2 * dp)
        d = (TH + math.pi) % (2 * math.pi) - math.pi
        phi = np.exp(-d**2 / (4 * sig**2) + 1j * p0 * d)
        phi /= np.linalg.norm(phi)
        c = V.conj().T @ phi
        ts = np.linspace(0, Tmax, 1201)
        coh = []
        for t in ts:
            ph = V @ (np.exp(-1j * E * t) * c)
            coh.append(abs(np.sum(np.abs(ph)**2 * np.exp(1j * TH))))
        coh = np.array(coh)
        # tập hợp cổ điển: cùng phân bố Gauss (θ, p_θ) ban đầu
        th0 = rng.normal(0, sig, 4000)
        pp = rng.normal(p0, dp, 4000)
        out = classical_traj(th0, pp / Mf(th0), Tmax)
        tc = np.array([o[0] for o in out]); cc = np.array([o[1] for o in out])
        c0 = coh[0]
        tq = ts[np.argmax(coh < c0 / math.e)] if np.any(coh < c0 / math.e) else float('nan')
        tcl = tc[np.argmax(cc < cc[0] / math.e)] if np.any(cc < cc[0] / math.e) else float('nan')
        i45 = np.argmin(np.abs(ts - 4.5 * T0)); j45 = np.argmin(np.abs(tc - 4.5 * T0))
        late = (ts > 1.5 * T0) & (ts < 2.6 * T0)               # quanh hồi sinh đầy đủ T_rev = 4πM̄ ≈ 2.1 T0
        ir = np.argmax(np.where(late, coh, 0))
        print(f"{dp:9.2f} {sig:9.3f} {tq/T0:18.3f} {tcl/T0:18.3f} {coh[i45]:20.3f} {cc[j45]:11.3f} "
              f"{coh[ir]:12.3f} ({ts[ir]/T0:.2f})")
    print("(LT = lượng tử, CĐ = tập hợp cổ điển Liouville cùng phân bố đầu; mất pha = |⟨e^iθ⟩| < 1/e giá trị đầu)\n")


if __name__ == "__main__":
    ev, Ebr = part_s1()
    part_s2(ev, Ebr)
    part_s3()
