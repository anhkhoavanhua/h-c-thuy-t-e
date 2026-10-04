# Trung bình hấp thụ mặt $\langle W_D\rangle$ trên hàm sóng và giá trị hiệu dụng của $\gamma_x$

> Tiếp nối [phan-bien-part-I-kiem-chung.md](phan-bien-part-I-kiem-chung.md), câu hỏi mở 1.
> Số liệu tính bằng [giai_so_wd_trung_binh.py](giai_so_wd_trung_binh.py) (chạy `python3 giai_so_wd_trung_binh.py`, vài giây).
> **Câu hỏi mở 2 (rà soát `paper_prc.tex`) chưa làm được:** file không có trong `logos_data/`, `Downloads/` hay `Desktop/`.

## Cách tính

**Thừa số dạng mặt** (đạo hàm Woods–Saxon, bằng 1 tại $r=R$):

$$g(r)=\frac{4e^{(r-R_D)/a_D}}{\big(1+e^{(r-R_D)/a_D}\big)^2}$$

**Thừa số dạng thể tích:** $f(r)=\big[1+e^{(r-R_V)/a_V}\big]^{-1}$.

**Hàm sóng:** dao động tử 3D với $\ell_0=2.04$ fm, $R_{nl}\propto r^le^{-r^2/2\ell_0^2}L_n^{l+1/2}(r^2/\ell_0^2)$.

**Hai hình học:**
- (i) $R=1.2A^{1/3}=4.59$ fm, $a=0.54$ fm, như đề xuất.
- (ii) Koning–Delaroche neutron: $R_D=4.90$ fm, $a_D=0.535$ fm; $R_V=4.58$ fm, $a_V=0.669$ fm.

**Độ cản hiệu dụng:**

$$\gamma_{\text{eff}}=\frac{2\big(W_V\langle f\rangle+W_D\langle g\rangle\big)}{\hbar\omega_0}$$

tại $\Delta E=7$ MeV, với $W_V=0.105$ MeV và $W_D=3.448$ MeV (A = 56). Giá trị đỉnh, tức khi $\langle g\rangle=\langle f\rangle=1$, là **0.711**.

## Kết quả (A = 56)

| Quỹ đạo | $N$ | $r_{\text{rms}}$ (fm) | (i) $\langle g\rangle$ | (i) $\gamma_{\text{eff}}$ | (ii) KD $\langle g\rangle$ | (ii) KD $\gamma_{\text{eff}}$ | So với đỉnh (KD) |
|---|---|---|---|---|---|---|---|
| 0s | 0 | 2.50 | 0.149 | 0.123 | 0.102 | 0.090 | −87% |
| 0p | 1 | 3.23 | 0.319 | 0.238 | 0.240 | 0.183 | −74% |
| 1s | 2 | 3.82 | 0.479 | 0.345 | 0.430 | 0.311 | −56% |
| 0d | 2 | 3.82 | 0.484 | 0.350 | 0.398 | 0.290 | −59% |
| **1p** | 3 | 4.33 | 0.460 | 0.329 | 0.475 | **0.340** | **−52%** |
| **0f** | 3 | 4.33 | 0.596 | 0.424 | 0.533 | **0.381** | **−46%** |
| 2s | 4 | 4.78 | 0.299 | 0.216 | 0.350 | 0.251 | −65% |
| 1d | 4 | 4.78 | 0.362 | 0.259 | 0.413 | 0.295 | −59% |
| **0g** | 4 | 4.78 | 0.630 | 0.445 | 0.612 | **0.432** | **−39%** |

Các quỹ đạo in đậm (1p, 0f, 0g) là quỹ đạo neutron ngay trên mức Fermi của ⁵⁶Fe ($N=30$).

## Kết luận

1. **Trung bình trên hàm sóng làm giảm độ cản đáng kể.** Với các quỹ đạo gần mức Fermi, $\gamma_{\text{eff}}=$ **0.22–0.45**, thấp hơn đỉnh 0.71 từ **37% đến 70%** tùy quỹ đạo và hình học. Quỹ đạo sâu (0s) chỉ còn 0.09–0.12.
2. **$\gamma_x=0.70$ cao hơn giá trị hiệu dụng khoảng 1.6–3 lần.** Muốn trung bình đạt 0.70, nucleon phải ở $\Delta E\approx$ **11–13 MeV** trên mặt Fermi, tức $E\approx+1$ tới $+3$ MeV (vùng tán xạ). Con số này cao hơn năng lượng 5–10 MeV của mô hình.

| Quỹ đạo, hình học | $\Delta E$ để $\gamma_{\text{eff}}=0.70$ | $E$ |
|---|---|---|
| 0f, KD | 13.0 MeV | +3.2 MeV |
| 0g, KD | 11.0 MeV | +1.2 MeV |
| 0f, $R=1.2A^{1/3}$ | 11.5 MeV | +1.7 MeV |

3. **Hệ quả cho Part I nếu hạ $\gamma_x$ về khoảng 0.3–0.45.** Tỉ số tốc độ dừng trên tốc độ gãy, $\dot\theta_d/\dot\theta_{\text{gãy}}=1.05\sqrt{1-\gamma_x\kappa^2/2\gamma_\theta}$, tăng từ **0.821** ($\gamma_x=0.70$) lên **0.909** ($0.45$) và **0.959** ($0.30$). Hệ sẽ nằm sát ngưỡng hơn nhiều, nên các thống kê chính (lề 13.7%, $\tau_{\text{lưu}}\approx4.5\,T_0$, $r\approx0.15/T_0$, bội số hình học 5.1) **sẽ thay đổi** và cần chạy lại. Ở tài liệu trước, riêng việc đổi $\gamma_x$ từ 0.7 xuống 0.5 đã làm $r$ tăng từ 0.16 lên 0.23 mỗi $T_0$.

## Đề xuất cho bài báo

Có hai lựa chọn:
- **(a)** Giữ $\gamma_x=0.70$ và ghi rõ đây là tham số hiệu dụng, ứng với nucleon ở khoảng 11–13 MeV trên mặt Fermi (sau khi lấy trung bình hàm sóng). Kèm một bảng độ nhạy các kết quả theo $\gamma_x$.
- **(b)** Hạ $\gamma_x$ về khoảng 0.35–0.45, phù hợp với quỹ đạo 0f và 0g tại $\Delta E\approx7$ MeV, rồi chạy lại toàn bộ thống kê của Part I.

## Câu hỏi còn mở

- **Chọn (a) hay (b)?** Nếu chọn (b), mình có thể chạy lại bộ thống kê chính (lề ứng suất, $\tau_{\text{lưu}}$, bội số, năng lượng xung, cửa sổ $f$) với $\gamma_x=0.4$.
- **File bản thảo.** Đặt `paper_prc.tex` vào `hoc-thuyet-e/` để làm câu hỏi mở 2.
