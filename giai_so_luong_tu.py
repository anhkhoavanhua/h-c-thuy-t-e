"""Kiểm tra mô hình E_kép bằng cơ học lượng tử: Schrödinger, Dirac, Page–Wootters, Wheeler–DeWitt.

Đơn vị ħ = m = ω0 = ℓ0 = 1 (năng lượng theo ħω0 = 10 MeV), trừ phần Dirac dùng MeV–fm.
  S1  Phổ lượng tử của L (1)–(2) (không liên kết) so với mode cổ điển ω± = √(1±g).
  S2  Liên kết x = κ sinθ dưới dạng giếng "bond" có lực cực đại đúng bằng R = 3: có giữ nổi hạt lượng tử
      trong phạm vi κ = 1/3 không? (trạng thái liên kết, độ rộng, so với κ).
  S3  TDSE 2 chiều (θ, x) với giếng bond: xác suất "gãy" theo mô-men động lượng rô-to ℓ, so với cổ điển.
  D   Dirac 1D (thế vô hướng dao động tử) cho proton: hiệu chỉnh tương đối tính so với Schrödinger.
  PW  Page–Wootters / Wheeler–DeWitt: dùng rô-to làm đồng hồ trong trạng thái toàn cục vượt thời gian;
      so trạng thái điều kiện của hạt tại số đọc θ với TDSE lái bởi θ = Ωt. Hai đồng hồ cùng Ω = 3:
      (I=1, ℓ0=3 — như mô hình, ~3ħ) và (I=10, ℓ0=30 — bán cổ điển, 30ħ).
Chạy: python3 giai_so_luong_tu.py [S1|S2|S3|D|PW|PW2]   (PW2: đồng hồ chậm Ω=0.3, g=1)   (cần numpy; S3 ≈ 1–2 phút)
"""
import math
import sys

import numpy as np

I_, M_, K_, DE, G, KAP, R = 1.0, 1.0, 1.0, 1.0, 0.1, 1 / 3, 3.0
T0 = 2 * math.pi


# ---------------------------------------------------------------- tiện ích cơ sở (ℓ, n)
def rotor_ops(ls):
    """Ma trận cosθ, sinθ trên cơ sở e^{iℓθ} với danh sách ℓ liên tiếp."""
    d = len(ls)
    up = np.diag(np.ones(d - 1), -1).astype(complex)   # e^{iθ}: |ℓ⟩ → |ℓ+1⟩
    dn = up.conj().T
    return (up + dn) / 2, (up - dn) / (2j)


def ho_x(nmax):
    a = np.diag(np.sqrt(np.arange(1, nmax + 1)), 1)
    return (a + a.T) / math.sqrt(2)


# ---------------------------------------------------------------- S1
def part_s1():
    ls = np.arange(-15, 16)
    nmax = 30
    cos_t, sin_t = rotor_ops(ls)
    x = ho_x(nmax)
    Hr = np.diag(ls**2 / (2 * I_)) - DE * cos_t
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    H = np.kron(Hr, np.eye(nmax + 1)) + np.kron(np.eye(len(ls)), Hx) - G * np.kron(sin_t, x)
    E = np.linalg.eigvalsh(H)[:6]
    wp, wm = math.sqrt(1 + G), math.sqrt(1 - G)
    E_rot = np.linalg.eigvalsh(Hr)[:3]
    print("=== S1. Phổ lượng tử của L (1)–(2) ===")
    print("6 mức thấp nhất:", ", ".join(f"{e:.4f}" for e in E))
    print(f"Kích thích E1−E0 = {E[1]-E[0]:.4f}, E2−E0 = {E[2]-E[0]:.4f}")
    print(f"Cổ điển tuyến tính: ω− = {wm:.4f}, ω+ = {wp:.4f}")
    print(f"Con lắc lượng tử riêng (g=0): mức kích thích đầu = {E_rot[1]-E_rot[0]:.4f} (điều hòa: 1.000) "
          f"⟹ phi điều hòa của thế −ΔE cosθ ở thang ħ\n")


# ---------------------------------------------------------------- S2
def bound_states(w, L=40.0, N=2000):
    y = np.linspace(-L, L, N)
    h = y[1] - y[0]
    D = R * w * math.sqrt(math.e)                      # max|dV/dy| = D/w·e^{-1/2} = R
    V = -D * np.exp(-y**2 / (2 * w * w))
    T = (np.diag(np.full(N, 1.0)) - 0.5 * np.diag(np.ones(N - 1), 1) - 0.5 * np.diag(np.ones(N - 1), -1)) / h**2
    E, U = np.linalg.eigh(T + np.diag(V))
    nb = int(np.sum(E < 0))
    rms = math.sqrt(np.sum(y**2 * U[:, 0]**2) / np.sum(U[:, 0]**2)) if nb else float('nan')
    return D, nb, E[0], rms


def part_s2():
    print("=== S2. Giếng bond có lực cực đại R = 3 (độ sâu D = R·w·√e) ===")
    print(f"{'w':>6} {'D (ħω0)':>9} {'số TT liên kết':>15} {'E0 (ħω0)':>10} {'Δy rms':>8} {'Δy/κ':>7}")
    for w in (0.05, 0.1, 0.2, 1 / 3, 0.5, 1.0, 2.0):
        D, nb, E0, rms = bound_states(w)
        print(f"{w:6.3f} {D:9.3f} {nb:15d} {E0:10.4f} {rms:8.3f} {rms/KAP:7.2f}")
    print(f"Muốn Δy ≲ κ/2 = {KAP/2:.3f} cần động năng điểm không ~ ħ²/(2mΔy²) = {1/(2*(KAP/2)**2):.1f} ħω0 "
          f"= {10/(2*(KAP/2)**2):.0f} MeV;  dao động tử tự do có Δx0 = 1/√2 = {1/math.sqrt(2):.3f} > κ\n")


# ---------------------------------------------------------------- S3
def part_s3(w=1.0, T=10 * T0, dt=0.005):
    D = R * w * math.sqrt(math.e)
    Nt, Nx, Lx = 64, 256, 12.0
    th = np.linspace(0, 2 * math.pi, Nt, endpoint=False)
    xs = np.linspace(-Lx, Lx, Nx, endpoint=False)
    TH, X = np.meshgrid(th, xs, indexing='ij')
    Y = X - KAP * np.sin(TH)
    V = -DE * np.cos(TH) + 0.5 * K_ * X**2 - G * X * np.sin(TH) - D * np.exp(-Y**2 / (2 * w * w))
    ll = np.fft.fftfreq(Nt, d=1 / Nt)
    kx = 2 * math.pi * np.fft.fftfreq(Nx, d=xs[1] - xs[0])
    KL, KX = np.meshgrid(ll, kx, indexing='ij')
    Kin = KL**2 / (2 * I_) + KX**2 / (2 * M_)
    eV = np.exp(-0.5j * dt * V)
    eK = np.exp(-1j * dt * Kin)
    wb = math.sqrt(K_ + D / (w * w))                   # tần số bond tại đáy
    sig = 1 / math.sqrt(2 * wb)
    core = np.abs(Y) > 2 * w
    print(f"=== S3. TDSE 2D (θ,x), giếng bond w={w}, D={D:.3f} ħω0, {T/T0:.0f} T0 ===")
    print(f"{'ℓ':>3} {'P_ngoài lõi TB (nửa sau)':>25} {'P_ngoài max':>12} {'cổ điển: gãy?':>14} {'t gãy cổ điển (T0)':>19}")
    steps = int(T / dt)
    for l in range(0, 7):
        psi = np.exp(1j * l * TH) * np.exp(-Y**2 / (4 * sig**2)) * np.exp(1j * KAP * l * np.cos(TH) * Y)
        psi /= np.sqrt(np.sum(np.abs(psi)**2))
        pout = []
        for s in range(steps):
            psi = eV * psi
            psi = np.fft.ifft2(eK * np.fft.fft2(psi))
            psi = eV * psi
            if s % 50 == 0:
                pout.append(float(np.sum(np.abs(psi[core])**2)))
        half = pout[len(pout) // 2:]
        tb, tc = classical_break(l, w, D, T)
        print(f"{l:3d} {sum(half)/len(half):25.4f} {max(pout):12.4f} {('CÓ' if tb else 'không'):>14} "
              f"{(tc/T0 if tb else float('nan')):19.3f}", flush=True)
    print("(P_ngoài = xác suất |x − κ sinθ| > 2w; cổ điển: cùng thế, θ̇0 = ℓ, x0 = 0, ẋ0 = κℓ)\n")


def classical_break(l, w, D, T, dt=1e-3):
    def acc(s):
        th, pt, x, px = s
        y = x - KAP * math.sin(th)
        fb = -D * y / (w * w) * math.exp(-y * y / (2 * w * w))   # −dV_b/dy
        dVdth = DE * math.sin(th) - G * x * math.cos(th) + fb * KAP * math.cos(th)
        dVdx = K_ * x - G * math.sin(th) - fb
        return (pt / I_, -dVdth, px / M_, -dVdx)
    s = [0.0, float(l), 0.0, KAP * l]
    t = 0.0
    while t < T:
        k1 = acc(s); k2 = acc([a + dt / 2 * b for a, b in zip(s, k1)])
        k3 = acc([a + dt / 2 * b for a, b in zip(s, k2)]); k4 = acc([a + dt * b for a, b in zip(s, k3)])
        s = [a + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4) for a, b1, b2, b3, b4 in zip(s, k1, k2, k3, k4)]
        t += dt
        if abs(s[2] - KAP * math.sin(s[0])) > 2 * w:
            return True, t
    return False, T


# ---------------------------------------------------------------- D
def part_d():
    hbarc, mc2 = 197.3269804, 938.27208816
    print("=== D. Dirac 1D (thế vô hướng S = ½kx²) so với Schrödinger, proton ===")
    for hw in (6.920, 10.0, 16.271):
        k = mc2 * hw**2 / hbarc**2                      # MeV/fm²
        l0 = hbarc / math.sqrt(mc2 * hw)
        L, N = 12 * l0, 256
        x = np.linspace(-L, L, N, endpoint=False)
        kk = 2 * math.pi * np.fft.fftfreq(N, d=x[1] - x[0])
        F = np.fft.fft(np.eye(N), axis=0)
        P = np.real_if_close(np.fft.ifft(kk[:, None] * F, axis=0))  # toán tử −i∂ (phổ)
        P = (P + P.conj().T) / 2
        S = 0.5 * k * x**2
        H = np.block([[np.diag(mc2 + S), hbarc * P], [hbarc * P, -np.diag(mc2 + S)]])
        Ed = np.linalg.eigvalsh(H)
        Ed = np.sort(Ed[Ed > 0])[:3] - mc2
        Hs = hbarc**2 * (P @ P) / (2 * mc2) + np.diag(S)
        Es = np.linalg.eigvalsh(Hs)[:3]
        rel = [(a - b) / b for a, b in zip(Ed, Es)]
        print(f"ħω0 = {hw:6.3f} MeV: Dirac E−mc² = " + ", ".join(f"{e:.4f}" for e in Ed)
              + " | Schrödinger = " + ", ".join(f"{e:.4f}" for e in Es)
              + " | lệch " + ", ".join(f"{r:+.3%}" for r in rel) + f" | ħω0/mc² = {hw/mc2:.3%}")
    print()


# ---------------------------------------------------------------- PW / WdW
def part_pw(g=0.3, clocks=((10.0, 30, 3.0), (1.0, 3, 1.0)), fs=(0.05, 0.1, 0.2)):
    Om_show = clocks[0][1] / clocks[0][0]
    print(f"=== PW/WdW. Rô-to làm đồng hồ Page–Wootters (ΔE = 0 cho đồng hồ sạch), g = {g}, Ω = {Om_show:g} ===")
    nmax = 15
    x = ho_x(nmax)
    Hx = np.diag(np.arange(nmax + 1) + 0.5)
    for I, l0, sl in clocks:
        Om = l0 / I
        Trot = 2 * math.pi / Om
        ls = np.arange(l0 - 24, l0 + 25)
        cos_t, sin_t = rotor_ops(ls)
        H = (np.kron(np.diag(ls**2 / (2 * I)), np.eye(nmax + 1)) + np.kron(np.eye(len(ls)), Hx)
             - g * np.kron(sin_t, x))
        E, U = np.linalg.eigh(H)
        c = np.exp(-(ls - l0)**2 / (4 * sl * sl)).astype(complex)
        c /= np.linalg.norm(c)
        chi0 = np.zeros(nmax + 1, complex); chi0[0] = 1
        phi0 = np.kron(c, chi0)
        amp = U.conj().T @ phi0
        E0 = float(E[np.argmax(np.abs(amp))])          # trị riêng có trọng số lớn nhất (nghiệm WdW HΨ = E0Ψ)
        # tham chiếu: TDSE của hạt lái bởi θ = Ωt
        def ref(t, n_sub=2000):
            psi = chi0.copy(); dt = t / n_sub
            for i in range(n_sub):
                tm = (i + 0.5) * dt
                Ht = Hx - g * math.sin(Om * tm) * x
                w_, v_ = np.linalg.eigh(Ht)
                psi = v_ @ (np.exp(-1j * w_ * dt) * (v_.conj().T @ psi))
            return psi
        print(f"-- đồng hồ I={I}, ℓ0={l0} (mô-men {l0}ħ), σ_ℓ={sl}, T_rot={Trot:.3f}")
        for label, tw in (("cửa sổ 0.35 T_rot (PaW, 1 vòng)", 0.35 * Trot), ("cửa sổ 4 T_rot", 4 * Trot),
                          ("WdW: cửa sổ 50 T_rot (≈ HΨ = EΨ)", 50 * Trot)):
            wts = np.exp(-((E - E0) * tw)**2 / 2) * amp
            Psi = (U @ wts).reshape(len(ls), nmax + 1)
            out = []
            for f in fs:
                t = f * Trot
                th = Om * t
                ph = np.exp(1j * ls * th)
                chi = ph @ Psi
                pc = np.linalg.norm(chi)**2
                chi /= np.linalg.norm(chi)
                r = ref(t)
                F = abs(np.vdot(r, chi))**2
                Fstat = abs(np.vdot(chi0, chi))**2
                out.append(f"t={f:.2f}T_rot: F={F:.4f} (tĩnh {Fstat:.4f})")
            # tính tuần hoàn: số đọc θ và θ + 2π cho cùng trạng thái
            th = Om * 0.1 * Trot
            a1 = np.exp(1j * ls * th) @ Psi; a2 = np.exp(1j * ls * (th + 2 * math.pi)) @ Psi
            same = abs(np.vdot(a1 / np.linalg.norm(a1), a2 / np.linalg.norm(a2)))**2
            print(f"   {label}: " + "; ".join(out) + f";  |⟨χ(θ)|χ(θ+2π)⟩|² = {same:.6f}")
        # quá khứ/tương lai trộn: trọng số của các vòng t, t±T_rot trong trạng thái WdW
        print(f"   ⟹ trọng số Gauss của các vòng lân cận trong cửa sổ 4 T_rot: e^(−T_rot²/2τ_w²) = {math.exp(-1/(2*16)):.3f}")
    print()


if __name__ == "__main__":
    parts = {"S1": part_s1, "S2": part_s2, "S3": part_s3, "D": part_d, "PW": part_pw,
             "PW2": lambda: part_pw(g=1.0, clocks=((100.0, 30, 3.0), (10.0, 3, 1.0)), fs=(0.1, 0.2, 0.3, 0.4))}
    for name in sys.argv[1:] or ["S1", "S2", "D", "PW", "PW2", "S3"]:
        parts[name.upper()]()
