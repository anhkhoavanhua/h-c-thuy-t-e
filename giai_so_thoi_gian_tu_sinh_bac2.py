"""Giải tiếp sàn 1 − F ≈ 0.05 của đồng hồ rô-to ở J ≈ 3ħ: phần nào là GÁN SAI THỜI GIAN, phần nào KHÔNG KHỬ ĐƯỢC.

Cùng mô hình N3 (giai_so_part3 / giai_so_dong_ho_tu_sinh): H = L²/2I + H_x − g sinθ x, I = J/Ω, g = 1, Ω = 0.3.
Rô-to mở cuộn (δ = 1/4, chu kỳ 4 T_rot) và cửa sổ tw = 5 T_rot để loại phần tái hồi.
Với mỗi góc đọc φ: χ(φ) = ⟨φ|Ψ⟩ (trạng thái hạt điều kiện), so với nghiệm Schrödinger ψ(t) của hạt:
  F_0(φ)    = |⟨ψ(φ/Ω)|χ(φ)⟩|²           thời gian WKB bậc 0: t = φ/Ω (rô-to quay đều)
  F_best(φ) = max_t |⟨ψ_ref(t)|χ(φ)⟩|²   gán lại thời gian tối ưu t*(φ) (giới hạn trên của mọi hàm t(φ))
  F_mf(φ)   = độ trùng với nghiệm tự hợp (Ehrenfest): hạt lái bởi θ_c(t), rô-to cổ điển chịu phản lực −g cosθ ⟨x⟩.
Phần không khử được = 1 − max_ψ∈{mọi trạng thái thuần của hạt có thể nối với χ(0) bằng Schrödinger} ≥ độ hỗn hợp:
  P_mix(φ) = 1 − Tr ρ_x(φ)² với ρ_x lấy trên cửa sổ góc nhỏ quanh φ (vướng víu rô-to–hạt).
Chạy: python3 giai_so_thoi_gian_tu_sinh_bac2.py   (≈ 2–4 phút)
"""
import math

import numpy as np


def build(J, delta=0.25, g=1.0, Om=0.3, nmax=10, half=16):
    I = J / Om
    step = int(round(1 / delta))
    ls = J + delta * np.arange(-half * step, half * step + 1)
    n = ls.size
    sh = np.zeros((n, n)); sh[np.arange(step, n), np.arange(n - step)] = 1.0
    sin_t = (sh - sh.T) / 2j
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    x = (a + a.T) / math.sqrt(2)
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    H = np.kron(np.diag(ls**2 / (2 * I)), np.eye(nmax + 1)) + np.kron(np.eye(n), Hx) - g * np.kron(sin_t, x)
    return dict(I=I, ls=ls, n=n, x=x, Hx=Hx, H=H, g=g, Om=Om, nmax=nmax)


def wdw_state(m, J, tw_rot=5.0):
    E, Uv = np.linalg.eigh(m['H'])
    sl = max(1.0, math.sqrt(J / 3))
    c = np.exp(-(m['ls'] - J)**2 / (4 * sl * sl)).astype(complex); c /= np.linalg.norm(c)
    chi0 = np.zeros(m['nmax'] + 1, complex); chi0[0] = 1
    amp = Uv.conj().T @ np.kron(c, chi0)
    E0 = float(E[np.argmax(np.abs(amp))])
    Trot = 2 * math.pi / m['Om']
    Psi = (Uv @ (np.exp(-((E - E0) * tw_rot * Trot)**2 / 2) * amp)).reshape(m['n'], m['nmax'] + 1)
    return Psi, chi0


def cond(m, Psi, phi):
    chi = np.exp(1j * m['ls'] * phi) @ Psi
    return chi / np.linalg.norm(chi)


def schrod_path(m, theta_of_t, T, nt=3000):
    """ψ(t) của hạt dưới H_x − g sin(θ(t)) x trên lưới t ∈ [0, T]."""
    psi = np.zeros(m['nmax'] + 1, complex); psi[0] = 1
    dt = T / nt; out = [psi.copy()]
    for i in range(nt):
        th = theta_of_t((i + 0.5) * dt)
        w_, v_ = np.linalg.eigh(m['Hx'] - m['g'] * math.sin(th) * m['x'])
        psi = v_ @ (np.exp(-1j * w_ * dt) * (v_.conj().T @ psi))
        out.append(psi.copy())
    return np.array(out), dt


def mean_field(m, J, T, nt=3000):
    """Ehrenfest tự hợp: I θ̈ = g cosθ ⟨x⟩ (phản lực của hạt lên rô-to), hạt lái bởi θ(t)."""
    psi = np.zeros(m['nmax'] + 1, complex); psi[0] = 1
    th, w = 0.0, J / m['I']
    dt = T / nt; ths, psis = [th], [psi.copy()]
    for _ in range(nt):
        xm = float(np.real(np.vdot(psi, m['x'] @ psi)))
        w += 0.5 * dt * m['g'] * math.cos(th) * xm / m['I']
        th += dt * w
        w_, v_ = np.linalg.eigh(m['Hx'] - m['g'] * math.sin(th - 0.5 * dt * w) * m['x'])
        psi = v_ @ (np.exp(-1j * w_ * dt) * (v_.conj().T @ psi))
        xm = float(np.real(np.vdot(psi, m['x'] @ psi)))
        w += 0.5 * dt * m['g'] * math.cos(th) * xm / m['I']
        ths.append(th); psis.append(psi.copy())
    return np.array(ths), np.array(psis)


def purity_window(m, Psi, phi, width=0.15, nq=9):
    """ρ_x khi đọc góc với độ phân giải hữu hạn (cửa sổ ±width rad): đo vướng víu rô-to–hạt còn sót."""
    rho = np.zeros((m['nmax'] + 1,) * 2, complex)
    for p in np.linspace(phi - width, phi + width, nq):
        chi = np.exp(1j * m['ls'] * p) @ Psi
        rho += np.outer(chi, chi.conj())
    rho /= np.trace(rho)
    return float(np.real(np.trace(rho @ rho)))


def run(J, g=1.0):
    m = build(J, g=g)
    Psi, chi0 = wdw_state(m, J)
    Om = m['Om']; Trot = 2 * math.pi / Om
    T = 0.6 * Trot
    ref, dt = schrod_path(m, lambda t: Om * t, T)
    ths, mf = mean_field(m, J, T)
    tgrid = np.arange(len(ref)) * dt
    rows = []
    for f in (0.1, 0.2, 0.3, 0.4):
        phi = f * 2 * math.pi
        chi = cond(m, Psi, phi)
        k0 = int(round(phi / Om / dt))
        F0 = abs(np.vdot(ref[k0], chi))**2
        ov = np.abs(ref.conj() @ chi)**2
        kb = int(np.argmax(ov)); Fb = ov[kb]
        kmf = int(np.argmin(np.abs(ths - phi)))                     # nghiệm tự hợp tại lúc θ_c(t) = φ
        Fmf = abs(np.vdot(mf[kmf], chi))**2
        pur = purity_window(m, Psi, phi)
        rows.append((f, 1 - F0, 1 - Fb, tgrid[kb] / (phi / Om), 1 - Fmf, tgrid[kmf] / (phi / Om), 1 - pur))
    print(f"\ng = {g}, J = {J}ħ (I = {m['I']:.1f}, chiều {m['n'] * (m['nmax'] + 1)}), rô-to mở cuộn δ = 1/4, tw = 5 T_rot")
    print("  φ/2π   1−F₀ (t=φ/Ω)   1−F_best   t*/(φ/Ω)   1−F_mf (Ehrenfest)   t_mf/(φ/Ω)   1−Tr ρ_x² (cửa sổ ±0.15)")
    for r in rows:
        print(f"  {r[0]:.1f}    {r[1]:.4f}         {r[2]:.4f}     {r[3]:.3f}      {r[4]:.4f}              {r[5]:.3f}        {r[6]:.4f}")
    a = np.array(rows)
    print(f"  max: 1−F₀ = {a[:, 1].max():.4f}, 1−F_best = {a[:, 2].max():.4f}, 1−F_mf = {a[:, 4].max():.4f}")
    return a


# ======================================================================= thông số đúng của E_dual (Bảng I)
def edual_check(J=3, g=0.1, dE=1.0, delta=1.0, tw_rot=50.0, nmax=10, half=12, nphi=8, alpha=0.0, I=1.0, sl=0.7, harm=1):
    """I = 1, ω0 = 1, θ̇ ≈ J/I ≈ 3ω0 (nhanh, không đoạn nhiệt), thế −ΔE cosθ, ghép −g x sinθ.
    Tham chiếu: hạt lái bởi rô-to CỔ ĐIỂN θ_c(t) cùng năng lượng (con lắc quay, không phản lực).
    Đọc góc φ = θ_c(t) trên vòng đầu; 1 − F = 1 − |⟨ψ(t)|χ(φ)⟩|²."""
    step = int(round(1 / delta))
    ls = J + delta * np.arange(-half * step, half * step + 1)
    n = ls.size
    sh = np.zeros((n, n)); sh[np.arange(step, n), np.arange(n - step)] = 1.0
    cos_t = (sh + sh.T) / 2
    shh = np.zeros((n, n)); shh[np.arange(harm * step, n), np.arange(n - harm * step)] = 1.0
    sin_t = (shh - shh.T) / 2j                                # sin(harm·θ): dời ℓ đi ±harm
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1); x = (a + a.T) / math.sqrt(2)
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    H = (np.kron(np.diag(ls**2 / (2 * I)) - dE * cos_t, np.eye(nmax + 1)) + np.kron(np.eye(n), Hx)
         - g * np.kron(sin_t, x))
    E, Uv = np.linalg.eigh(H)
    c = np.exp(-(ls - J)**2 / (4 * sl * sl)).astype(complex); c /= np.linalg.norm(c)
    chi0 = np.array([alpha**k / math.sqrt(math.factorial(k)) for k in range(nmax + 1)], complex)
    chi0 /= np.linalg.norm(chi0)                              # trạng thái kết hợp |α⟩ (α = 0: cơ bản)
    amp = Uv.conj().T @ np.kron(c, chi0)
    E0 = float(E[np.argmax(np.abs(amp))])                     # tâm cửa sổ = trị riêng có biên độ lớn nhất (như N3)
    Psi = (Uv @ (np.exp(-((E - E0) * tw_rot * 2 * math.pi / (J / I))**2 / 2) * amp)).reshape(n, nmax + 1)
    # rô-to cổ điển: ½Iθ̇² − ΔE cosθ = ½I(J/I)² − ΔE  (bắt đầu θ = 0, ℓ = J)
    Erot = 0.5 * J * J / I - dE
    T1 = 0.0; th = 0.0; dt = 2e-4; ths = [0.0]; ts = [0.0]
    while th < 2 * math.pi * 0.95:
        w = math.sqrt(2 * (Erot + dE * math.cos(th)) / I); th += w * dt; T1 += dt
        ths.append(th); ts.append(T1)
    ths = np.array(ths); ts = np.array(ts)
    psi = chi0.copy(); out = []; targets = [2 * math.pi * (k + 1) / (nphi + 1) for k in range(nphi)]
    k = 0
    for i in range(1, len(ts)):
        thm = 0.5 * (ths[i] + ths[i - 1])
        w_, v_ = np.linalg.eigh(Hx - g * math.sin(harm * thm) * x)
        psi = v_ @ (np.exp(-1j * w_ * (ts[i] - ts[i - 1])) * (v_.conj().T @ psi))
        while k < nphi and ths[i] >= targets[k]:
            chi = np.exp(1j * ls * targets[k]) @ Psi; chi /= np.linalg.norm(chi)
            out.append((1 - abs(np.vdot(psi, chi))**2, 1 - abs(np.vdot(psi, chi0))**2,
                        1 - abs(np.vdot(psi, _free(Hx, chi0, ts[i])))**2)); k += 1
    o = np.array(out)
    return o[:, 0].max(), float(o[:, 0].mean()), n * (nmax + 1), o[:, 1].max(), o[:, 2].max()


def _free(Hx, chi0, t):
    return np.exp(-1j * np.diag(Hx) * t) * chi0


def clock_condition():
    """Điều kiện phổ của đồng hồ: bậc năng lượng rô-to ħθ̇ phải ≪ lượng tử của hạt ħω0.
    Quét θ̇ = J/I (ΔE = 0 để rô-to quay đều), hạt chuyển động thật (α = 0.707i), g = 0.1, S¹, WdW tw = 50."""
    print("\n=== Điều kiện phổ đồng hồ: quét θ̇ = J/I (ΔE = 0, α = 0.707i, g = 0.1, S¹, tw = 50 T_rot) ===")
    print("   J    I      θ̇/ω0   max 1−F   TB 1−F   | nền 1−|⟨χ0|ψ⟩|²")
    for J in (3, 10, 30):
        for thd in (3.0, 1.0, 0.3, 0.1):
            I = J / thd
            half = 12 if J == 3 else 14
            mx, mean, dim, b0, bf = edual_check(J=J, g=0.1, dE=0.0, delta=1.0, tw_rot=50.0, alpha=0.707j, I=I,
                                                half=half)
            print(f"  {J:3d}  {I:6.1f}  {thd:5.2f}    {mx:.4f}    {mean:.4f}   |  {b0:.4f}", flush=True)


def clock_width():
    """Điều kiện Page–Wootters: độ trải năng lượng đồng hồ θ̇·σ_ℓ phải phủ độ trải năng lượng của hạt (≈ |α|ħω0)."""
    print("\n=== Độ rộng bó sóng rô-to σ_ℓ (J = 30, ΔE = 0, α = 0.707i, g = 0.1, S¹, tw = 50 T_rot) ===")
    print("  θ̇/ω0   σ_ℓ    θ̇·σ_ℓ (ħω0)   max 1−F   TB 1−F")
    for thd in (3.0, 1.0, 0.3):
        for sl in (0.7, 2.0, 5.0, 10.0):
            mx, mean, dim, b0, bf = edual_check(J=30, g=0.1, dE=0.0, delta=1.0, tw_rot=50.0, alpha=0.707j,
                                                I=30 / thd, half=max(14, int(4 * sl) + 4), sl=sl)
            print(f"  {thd:5.2f}  {sl:5.1f}   {thd * sl:7.2f}        {mx:.4f}    {mean:.4f}", flush=True)


def commensurate():
    """Luật cộng hưởng: đồng hồ rô-to S¹ tái lập Schrödinger của hạt khi k·θ̇ ≈ ω0 (k nguyên), vì trong WdW
    bước nhảy ħω0 của hạt phải được trả bằng k bậc ħθ̇ của rô-to."""
    print("\n=== Cộng hưởng phổ đồng hồ–hạt (J = 30, σ_ℓ = 5, ΔE = 0, α = 0.707i, g = 0.1, S¹) ===")
    print("  θ̇/ω0    ω0/θ̇    lệch khỏi ω0/k gần nhất    max 1−F (tw=50)   (tw=5)")
    for thd in (3.0, 2.0, 1.5, 1.0, 0.7, 0.5, 1 / 3, 0.3, 0.25):
        r = 1 / thd; mis = abs(r - round(r)) if r >= 1 else abs(thd - round(thd))
        out = []
        for tw in (50.0, 5.0):
            mx, mean, dim, b0, bf = edual_check(J=30, g=0.1, dE=0.0, delta=1.0, tw_rot=tw, alpha=0.707j,
                                                I=30 / thd, half=24, sl=5.0)
            out.append(mx)
        print(f"  {thd:6.3f}   {r:6.3f}       {mis:.3f}                   {out[0]:.4f}            {out[1]:.4f}", flush=True)


def resonance_mechanism():
    """Cơ chế: trong WdW (một mặt năng lượng, phổ rời rạc) hạt chỉ 'chuyển động' theo góc khi các cặp (ℓ, n) và
    (ℓ − h, n + 1) gần suy biến VÀ được tương tác nối (sin hθ dời ℓ đi ±h). Dự đoán:
      (1) độ rộng cộng hưởng quanh hθ̇ = ω0 tỉ lệ g;   (2) với ghép sin2θ, θ̇ = ω0/2 thành công."""
    print("\n=== Cơ chế cộng hưởng (J = 30, σ_ℓ = 5, ΔE = 0, α = 0.707i, S¹, tw = 50 T_rot) ===")
    print("  (1) quét θ̇ quanh ω0, ghép sinθ:")
    print("   θ̇/ω0    g=0.1     g=0.3")
    for thd in (0.9, 0.95, 0.98, 1.0, 1.02, 1.05, 1.1):
        r = [edual_check(J=30, g=g, dE=0.0, tw_rot=50.0, alpha=0.707j, I=30 / thd, half=24, sl=5.0)[0]
             for g in (0.1, 0.3)]
        print(f"   {thd:5.2f}    {r[0]:.4f}    {r[1]:.4f}", flush=True)
    print("  (2) đổi ghép thành sin(hθ): θ̇ thành công phải là ω0/h")
    print("   h   θ̇/ω0   max 1−F (g = 0.1)")
    for h, thd in ((1, 1.0), (1, 0.5), (2, 0.5), (2, 1.0), (3, 1 / 3), (1, 1 / 3)):
        r = edual_check(J=30, g=0.1, dE=0.0, tw_rot=50.0, alpha=0.707j, I=30 / thd, half=24, sl=5.0, harm=h)[0]
        print(f"   {h}   {thd:5.3f}  {r:.4f}", flush=True)


def edual_table():
    print("\n=== Thông số E_dual: I = 1, θ̇ ≈ 3ω0, ΔE = 1, J = 3ħ, σ_ℓ = 0.7 ===")
    print("  α      g     δ (rô-to)        tw (T_rot)   max 1−F    TB 1−F   | nền: 1−|⟨χ0|ψ⟩|²  1−|⟨ψ_tự do|ψ⟩|² (ảnh hưởng của ghép)")
    for alpha in (0.0, 0.707j):
        for g in (0.1, 0.3, 1.0):
            for delta, tw in ((1.0, 0.35), (1.0, 50.0), (0.25, 5.0)):
                mx, mean, dim, b0, bf = edual_check(g=g, delta=delta, tw_rot=tw, alpha=alpha)
                lab = "S¹ (học thuyết)" if delta == 1 else "mở cuộn ×4"
                print(f"  {abs(alpha):.3f}  {g:4.1f}  {lab:16s} {tw:6.2f}       {mx:.4f}     {mean:.4f}   |  {b0:.4f}            {bf:.4f}", flush=True)


if __name__ == '__main__':
    import sys
    if 'edual' in sys.argv:
        edual_table(); sys.exit()
    if 'clock' in sys.argv:
        clock_condition(); sys.exit()
    if 'width' in sys.argv:
        clock_width(); sys.exit()
    if 'comm' in sys.argv:
        commensurate(); sys.exit()
    if 'mech' in sys.argv:
        resonance_mechanism(); sys.exit()
    summary = []
    for g in (1.0, 0.3, 0.1):                       # g = 1: mô hình N3; g = 0.1: hệ số ghép của E_dual (Bảng I)
        for J in (3, 6, 12, 30):
            a = run(J, g)
            summary.append((g, J, a[:, 1].max(), a[:, 2].max(), a[:, 4].max()))
    print("\nTóm tắt max 1 − F:   g     J    t=φ/Ω     t*(φ) tối ưu   Ehrenfest   J·(1−F_best)")
    for g, J, f0, fb, fm in summary:
        print(f"                     {g:4.1f}  {J:3d}   {f0:.4f}    {fb:.4f}         {fm:.4f}      {J * fb:.3f}")
