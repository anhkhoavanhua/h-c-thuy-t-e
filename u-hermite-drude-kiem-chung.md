# Toán tử vận tốc tương đối $\hat u$ và bể Drude–Lorentz: giải số hai câu hỏi mở của Part II

> Tiếp nối [lindblad-joint-kiem-chung.md](lindblad-joint-kiem-chung.md).
> Số liệu tính bằng [giai_so_u_hermite_drude.py](giai_so_u_hermite_drude.py) (cần numpy). Chạy `python3 giai_so_u_hermite_drude.py`, ≈ 5 phút.
> **Đơn vị:** $\hbar=m=\omega_0=1$ (năng lượng theo $\hbar\omega_0=10$ MeV), $I=1$, $\kappa=1/3$, $\Delta E=1$, $g=0.1$, $\gamma=0.7$.
> **Cơ sở tính:** $\ell\in[-8,8]$ cho rô-to, $n\le6$ cho hạt (phần A–C).

## Tóm tắt

Kết luận của đề xuất ("chưa đủ để viết Part II, nên đưa vào mục Outlook của Part I") **được số liệu ủng hộ**. Các số liệu cụ thể:

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\hat u=\tfrac{\kappa}{2I}\{\cos\theta,\hat p_\theta\}-\hat p_x/m$ là Hermite | ✓ ($\lVert\hat u-\hat u^\dagger\rVert=0$) |
| Dạng khai triển $\tfrac{\kappa}{I}\cos\theta\,\hat p_\theta\,\mathbf{-}\,\tfrac{i\hbar\kappa}{2I}\sin\theta$ | **Sai dấu.** Vì $[\hat p_\theta,\cos\theta]=+i\hbar\sin\theta$, dạng đúng là $\tfrac\kappa I\cos\theta\,\hat p_\theta\,\mathbf{+}\,\tfrac{i\hbar\kappa}{2I}\sin\theta$. Kiểm bằng ma trận: dấu "+" khớp tuyệt đối (sai lệch 0); dấu "−" lệch 2.65 và **không Hermite** (sai lệch 5.29) |
| $L$ Hermite chỉ gây mất kết hợp và **nung nóng**, không làm nguội | **✓ Được xác nhận.** Với $L=\sqrt{0.2}\,\hat u$, xuất phát từ $\lvert\ell=3,n=0\rangle$: $\langle H\rangle$ = 5.00 → 5.65 → 6.33 → 7.72 (tại $t$ = 0, 2.5, 5, 10), **tăng đơn điệu**, hướng tới trần 15.5 của cơ sở cắt |
| Cần toán tử hạ mức $L=\sqrt\Gamma\hat P\hat A_-$; kiểm bằng quỹ đạo lượng tử xem có ra phân bố hình học ($p\approx0.2$) và thời gian chờ mũ không | **Đã mô phỏng** với dạng chuẩn **Davies ở $T=0$** ($\hat A_-$ = phần hạ năng lượng của $\hat u$ trong cơ sở riêng của $H$), 300 quỹ đạo MCWF. **Không tái lập được thống kê cổ điển:** số bước nhảy gần như cố định (**2** bước ở 74% quỹ đạo, CV 0.25, trong khi hình học cùng trung bình có CV 0.72); **mỗi bước mất 1.70 $\hbar\omega_0$** (gấp 4 lần Carnot 0.40); thời gian chờ **không mũ** (CV 1.5) |
| Bể Drude–Lorentz với $\hbar\omega_c\approx4.5$ MeV | **Mâu thuẫn định lượng với Part I.** Vì $\omega_c=0.45\,\omega_0<\omega_0$, ma sát hiệu dụng $\mathrm{Re}\,\hat\gamma(\omega)=\gamma\omega_c^2/(\omega_c^2+\omega^2)$ chỉ còn **16.8%** của γ tại $\omega_0$ và **2.8%** tại $\Omega=2.63$. Bể này **không** tạo ra được ma sát $\gamma_x=0.7$ mà Part I dùng |
| Thời gian nhớ $1/\omega_c\approx1.5\times10^{-22}$ s bằng khoảng cách hai xung Carnot | ✓ Đúng số ($1.46\times10^{-22}$ so với $1.47\times10^{-22}$ s). Đây là trùng hợp, vì cả hai là lựa chọn độc lập |
| "Mọi Lindblad Markov đều là gần đúng thô; **bắt buộc** HEOM" | **Đúng một phần.** Với dao động tử **tuyến tính** ghép bể Drude, trạng thái cân bằng giải **chính xác** được mà không cần HEOM (định lý thăng giáng–tiêu tán). Kết quả: $\langle x^2\rangle=$ **0.473**, $\langle p^2\rangle=0.561$, tích 0.265 ≥ 0.25 (vật lý). So với nghiệm chính xác: Caldeira–Leggett 0.100 (sai, ρ âm), Lindblad–Diósi 0.782 (**+65%**), Lindblad nhiệt $a/a^\dagger$ 0.500 (**+6%**, gần nhất). HEOM chỉ thực sự **bắt buộc** khi thêm phần **phi tuyến** (liên kết, va chạm) |
| Giữ Part I ở cấp cổ điển và bán cổ điển; đưa $\hat u$ và bể không Markov vào "Outlook" | **Đồng ý.** Số liệu ở trên cho thấy chưa có ứng viên lượng tử nào tái lập được thống kê của Part I |

## 1. Dấu trong $\hat u$ (A)

$$\tfrac12\{\cos\theta,\hat p_\theta\}=\cos\theta\,\hat p_\theta+\tfrac12[\hat p_\theta,\cos\theta]=\cos\theta\,\hat p_\theta+\tfrac{i\hbar}{2}\sin\theta$$

$$\Longrightarrow\quad\hat u=\frac\kappa I\cos\theta\,\hat p_\theta+\frac{i\hbar\kappa}{2I}\sin\theta-\frac{\hat p_x}{m}$$

| Kiểm tra (ma trận, bỏ hàng/cột biên của cơ sở cắt) | Chuẩn sai lệch |
|---|---|
| $\tfrac12\{\cos,p\}-(\cos p+\tfrac i2\sin)$ | 0 |
| $\tfrac12\{\cos,p\}-(\cos p-\tfrac i2\sin)$ (đề xuất) | 2.65 |
| Tính Hermite của dạng đề xuất | 5.29 (không Hermite) |

## 2. $L$ Hermite nung nóng hệ (B)

$\Gamma=0.2$, xuất phát $\lvert\ell=3,n=0\rangle$, $\langle H\rangle(0)=5.00$:

| $t$ | 0 | 2.5 | 5 | 10 |
|---|---|---|---|---|
| $\langle H\rangle$ | 5.001 | 5.646 | 6.325 | 7.716 |

## 3. Toán tử hạ mức dạng Davies: quỹ đạo lượng tử (C)

**Thiết lập:**
- $L=\sqrt{0.5}\,\hat u_-$, với $\hat u_-=\sum_{E_j<E_i}\langle j\vert\hat u\vert i\rangle\,\lvert j\rangle\langle i\rvert$.
- Xuất phát $\lvert\ell=3,n=0\rangle$ ($E=5.00$). Trạng thái cơ bản có $E=-0.036$.
- Dừng khi năng lượng cách trạng thái cơ bản dưới $0.05$.

| Chỉ số | Lượng tử (Davies) | Cổ điển (Part I) |
|---|---|---|
| Số bước nhảy tới trạng thái dừng | TB 2.08; phân bố 1:28, **2:223**, 3:47, 4:2; CV 0.25 | Hình học, TB 5.12, $p=0.195$, CV ≈ 0.90 |
| Năng lượng mỗi bước | TB **1.70** $\hbar\omega_0$; 10–90%: 0.67–2.43; chỉ 8.2% dưới 0.5 | TB 0.40 $\hbar\omega_0$ |
| Thời gian chờ | CV **1.51** (không mũ) | Mũ (CV ≈ 1) |

**Đọc bảng:**
- Hệ lượng tử đi xuống thang năng lượng bằng **vài bước lớn**, rời rạc theo phổ của $H$. Nó không đi bằng nhiều bước nhỏ 0.4 $\hbar\omega_0$.
- Mô phỏng này **không có nguồn lái và không có nhiễu**, nên chỉ kiểm được một lần hồi phục, không kiểm được chu trình.
- Dù vậy, ngay cả một lần hồi phục đã cho thống kê khác hẳn cổ điển. Muốn so trọn vẹn, cần thêm nguồn lái lượng tử (ví dụ bơm $\hat\ell$ tăng) và bể nhiệt.

## 4. Bể Drude–Lorentz (D)

### 4.1. Ma sát hiệu dụng theo tần số

$$\mathrm{Re}\,\hat\gamma(\omega)=\gamma\,\frac{\omega_c^2}{\omega_c^2+\omega^2}$$

| $\omega$ | $\mathrm{Re}\,\hat\gamma$ | So với $\gamma=0.7$ |
|---|---|---|
| $\omega_0=1$ | 0.118 | 16.8% |
| $\Omega=2.63$ | 0.020 | **2.8%** |

### 4.2. Cân bằng chính xác của dao động tử trong bể Drude ($k_BT=1$ MeV)

Công thức:

$$\langle x^2\rangle=\frac1\pi\int_0^\infty\!d\omega\,\coth\frac{\omega}{2k_BT}\,\mathrm{Im}\,\chi(\omega),\qquad \chi(\omega)=\frac{1}{\omega_0^2-\omega^2-i\omega\hat\gamma(\omega)}$$

Kết quả theo tần số cắt:

| Tần số cắt | $\langle x^2\rangle$ | $\langle p^2\rangle$ | $\langle x^2\rangle\langle p^2\rangle$ |
|---|---|---|---|
| $\omega_c=0.45$ (4.5 MeV) | **0.4726** | 0.5612 | 0.265 |
| $\omega_c=5$ (50 MeV) | 0.4349 | 0.7935 | 0.345 |
| $\omega_c=20$ | 0.4250 | 1.0164 | 0.432 |

So các xấp xỉ Markov với nghiệm chính xác ($\omega_c=0.45$):

| Xấp xỉ Markov | $\langle x^2\rangle$ | So với chính xác 0.473 |
|---|---|---|
| Caldeira–Leggett | 0.100 | −79%, ρ âm |
| Lindblad–Diósi | 0.782 | +65% |
| Lindblad nhiệt $a/a^\dagger$ | 0.500 | **+6%** |

## 5. Đoạn "Outlook" đề xuất cho Part I (có số liệu)

> "A quantum extension of the present model faces three quantified obstacles. (i) The Hermitian relative-velocity operator $\hat u=\tfrac{\kappa}{2I}\{\cos\theta,\hat p_\theta\}-\hat p_x/m$ used as a jump operator only heats the system; a cooling (Davies-type) projection $\hat u_-$ relaxes the rotor in about two large jumps of $\approx1.7\,\hbar\omega_0$ each, rather than in geometric cascades of $\approx0.4\,\hbar\omega_0$ impacts. (ii) A Drude–Lorentz bath with a cutoff set by the giant-resonance width ($\hbar\omega_c\approx4.5$ MeV) retains only 17% (at $\omega_0$) and 3% (at $\Omega\approx2.6\,\omega_0$) of the classical friction assumed here. (iii) For the linear part, exact quantum Brownian motion gives $\langle x^2\rangle=0.47$, compared with 0.10 (Caldeira–Leggett, non-positive) and 0.78 (Diósi–Lindblad). A faithful quantum treatment of the constrained, impacting system therefore requires non-Markovian methods (e.g., HEOM) together with a quantum model of the drive; this is left for future work."

## Câu hỏi còn mở

- **Nguồn lái lượng tử.** Muốn so chu trình (không chỉ một lần hồi phục), cần một cơ chế bơm năng lượng cho rô-to trong mô tả lượng tử, ví dụ một bể nóng thứ hai hoặc một lực lái tuần hoàn. Học thuyết muốn chọn cơ chế nào?
- **Tần số cắt $\omega_c$.** Lấy độ rộng GDR (4–5 MeV) làm $\omega_c$ làm ma sát gần như biến mất ở các tần số quan trọng. Nếu muốn giữ $\gamma_x\approx0.7$, cần $\omega_c\gg\Omega$, tức $\hbar\omega_c\gtrsim50$ MeV. Thang này ứng với quá trình vật lý nào?
