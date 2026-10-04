# Cản bức xạ $\gamma_\theta$, Langevin đầy đủ 2 bậc tự do và đối chiếu 7.6–7.9 MeV: giải và kiểm chứng số liệu

> Tiếp nối [dong-co-langevin-kiem-chung.md](dong-co-langevin-kiem-chung.md).
> Số liệu tính bằng [giai_so_langevin_day_du.py](giai_so_langevin_day_du.py) (chạy `python3 giai_so_langevin_day_du.py`, ≈ 6 phút trên 8 lõi).
> [giai_so_dong_co_langevin.py](giai_so_dong_co_langevin.py) được thêm tham số `gx`: lực cản $\gamma_x\dot x$ và lực Langevin $\xi_x$ lên hạt. Trong pha liên kết, $\xi_x$ được chiếu lên $\theta$ qua $\partial x/\partial\theta=\kappa\cos\theta$. Khi `gx=0`, kết quả cũ giữ nguyên.
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}$. Một đơn vị năng lượng bằng $\hbar\omega_0=10$ MeV.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\gamma_\theta=\tfrac23\alpha_{\text{eff}}\hbar\omega_0^2/c^2$ | **Sai thứ nguyên.** $\hbar\omega_0^2/c^2$ có đơn vị kg/s, trong khi $\gamma_\theta$ cần kg·m²/s (thiếu m²). Ngoài ra, phản lực Abraham–Lorentz tỉ lệ với **đạo hàm bậc ba** của vị trí; nó chỉ quy được về ma sát tuyến tính cho một dao động tử gần cộng hưởng |
| $\tau=I/\gamma_\theta\approx10$–$50$ | **Sai 2–3 bậc độ lớn.** Tốc độ tắt bức xạ chuẩn là $\Gamma=\tfrac23\alpha\hbar\omega_0^2/(mc^2)$, nên $\Gamma/\omega_0=\tfrac23\alpha\,\hbar\omega_0/(m_pc^2)$. Cho: $\tau=$ **19 287** ($\alpha=1/137$), **14 074** ($\alpha=0.01$), **1 407** ($\alpha=0.1$), tức 224–3 070 $T_0$. Cản bức xạ **rất yếu** |
| "Vùng $\tau\approx10$ là vùng bức xạ thực tế, loại $\tau=1$" | **Không có căn cứ bức xạ:** $\tau=10$ vẫn nhỏ hơn giá trị bức xạ 140–1 900 lần. Thử $\tau=1000$ (gần giá trị bức xạ): 53/55 lần tái khóa kịp (mục 2) |
| Langevin 2 bậc tự do: $\xi_x$, $\gamma_x=m/\tau_{\text{relax}}=0.789$ | **Đúng dạng.** Đã mô phỏng |
| "Cản $\gamma_x$ mạnh giúp bắt dính hạt nhanh hơn" | **Đúng với $\gamma_x=0.789$:** thời gian tái khóa trung bình giảm từ 2.555 xuống **1.03–1.69**. (Ở tài liệu trước, $\gamma_x=0.1$ chỉ nhanh hơn khoảng 10%, và còn làm tệ hơn khi $e=0.2$) |
| "Xác nhận dứt điểm tái khóa 100%" | **Tái khóa: 100%** (63/63, 63/63, 55/55). **Kịp trước khi mất pha: không đạt 100%**, chỉ **62/63 (98%)** với nhiễu 1 MeV, **55/63 (87%)** với nhiễu 5 MeV, và **53/55 (96%)** với $\tau=1000$ |
| 7.6–7.9 MeV/xung "phổ quát" | **Sụp khi bật $\gamma_x$.** Với đúng cấu hình đề xuất ($\gamma_x=0.789$), Carnot mỗi lần khóa chỉ còn **3.72–4.13 MeV**, vì cản đã tiêu tán năng lượng của hạt trước va chạm. Chỉ khi $\gamma_x=0$ mới giữ được 7.56 MeV |
| 7.6–7.9 MeV $=B/A\approx7.5$–$8.5$ MeV | **Trùng hợp số học.** Xung trong mô hình bằng $(0.764$–$0.793)\,\hbar\omega_0$, tức **tỉ lệ với $\hbar\omega_0=41A^{-1/3}$**, trong khi $B/A$ gần như không đổi. Với ¹⁶O: 12.4–12.9 so với $B/A=7.98$. Với ²⁰⁸Pb: 5.3–5.5 so với 7.87. Chỉ khớp khi $A\approx69$, tức đúng giá trị $\hbar\omega_0=10$ MeV đã **chọn tay** |
| ⁵⁶Fe(n,γ) có vạch 7.63 và 7.64 MeV | **Đúng là vạch thực nghiệm.** Nhưng năng lượng đó là năng lượng tách neutron của ⁵⁷Fe, tức năng lượng liên kết của neutron vừa bị bắt. Mô hình (với $\hbar\omega_0^{\text{Fe}}=10.72$ MeV) cho 8.19–8.50 MeV khi không cản, và cỡ 4 MeV khi có $\gamma_x$ |

## 1. Cản bức xạ: $\tau$ thực sự

$$\Gamma=\frac23\,\frac{\alpha_{\text{eff}}\,\hbar\omega_0^2}{m_pc^2},\qquad \frac{\Gamma}{\omega_0}=\frac23\,\alpha_{\text{eff}}\,\frac{\hbar\omega_0}{m_pc^2}$$

Với $\hbar\omega_0=10$ MeV, $m_pc^2=938.27$ MeV:

| $\alpha_{\text{eff}}$ | $\Gamma/\omega_0$ | $\tau=\omega_0/\Gamma$ (đơn vị $1/\omega_0$) | Theo $T_0$ |
|---|---|---|---|
| 1/137.036 | $5.185\times10^{-5}$ | 19 287 | 3 070 |
| 0.01 | $7.105\times10^{-5}$ | 14 074 | 2 240 |
| 0.1 | $7.105\times10^{-4}$ | 1 407 | 224 |
| *Đề xuất* | — | *10–50* | *1.6–8* |

## 2. Monte Carlo Langevin đầy đủ ($N=64$, cửa sổ $10T_0$)

**Cách đo:** $\tau_{\text{deph}}$ được đo lại cho từng cấu hình, là thời điểm $\langle\Delta\theta^2\rangle=1$ khi không cho gãy. "Kịp" nghĩa là tái khóa trong vòng $\tau_{\text{deph}}$ đo được, tính từ lần gãy đầu tiên.

| $\tau$ | $\gamma_x$ | Nhiễu | Gãy ≥1 lần | Tái khóa | $t$ tái khóa min / TB / max | $\tau_{\text{deph}}$ đo | Kịp | Gãy TB mỗi quỹ đạo | **Carnot mỗi lần khóa** |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0 | 5 MeV | 64/64 | 64/64 | 1.431 / 2.555 / 3.254 | 4.00 | 64/64 | 15.1 | **7.56 MeV** |
| 10 | 0.789 | 1 MeV | 63/64 | 63/63 | 0.002 / 1.687 / 5.824 | 5.50 | **62/63** | 10.6 | **3.99 MeV** |
| 10 | 0.789 | 5 MeV | 63/64 | 63/63 | 0.002 / 1.033 / 6.169 | 3.00 | **55/63** | 22.9 | **4.13 MeV** |
| 1000 | 0.789 | 5 MeV | 55/64 | 55/55 | 0.002 / 0.967 / 5.997 | 4.50 | **53/55** | 5.2 | **3.72 MeV** |

**Đọc kết quả:**

1. **Tái khóa vẫn 100%** sau mỗi lần gãy. Cản $\gamma_x$ làm tái khóa nhanh hơn trung bình, và có những lần khóa gần như tức thì ($t=0.002$).
2. **Đuôi phân bố dài hơn:** thời gian tái khóa tối đa tăng từ 3.25 lên 5.8–6.2. Một phần quỹ đạo vì thế **không kịp** trước khi mất pha: 2–13% tùy cấu hình.
3. **Năng lượng xung giảm khoảng một nửa** (7.56 → 3.7–4.1 MeV). Năng lượng lái giờ tiêu tán qua **hai kênh**: ma sát liên tục $\gamma_x\dot x^2$ và va chạm Carnot. Như vậy "xung 7.6–7.9 MeV" không phải là hằng số của mô hình, mà phụ thuộc $\gamma_x$.
4. **Với $\tau=1000$** (gần giá trị bức xạ), mô-men lái rất yếu. Chỉ 55/64 quỹ đạo gãy trong $10T_0$, trung bình 5.2 lần gãy mỗi quỹ đạo. Nhưng những quỹ đạo đã gãy thì vẫn tái khóa kịp 96%.

## 3. Năng lượng xung so với $B/A$

Xung mô hình (khi $\gamma_x=0$) $=(0.764$–$0.793)\times\hbar\omega_0$, với $\hbar\omega_0=41A^{-1/3}$ MeV.

| Hạt nhân | $\hbar\omega_0$ (MeV) | Xung mô hình (MeV) | $B/A$ thực nghiệm (MeV) |
|---|---|---|---|
| ¹⁶O | 16.271 | 12.43–12.90 | 7.976 |
| ⁴⁰Ca | 11.988 | 9.16–9.51 | 8.551 |
| ⁵⁶Fe | 10.717 | 8.19–8.50 | 8.790 |
| ¹²⁰Sn | 8.312 | 6.35–6.59 | 8.505 |
| ²⁰⁸Pb | 6.920 | 5.29–5.49 | 7.867 |

- Mô hình biến thiên **2.4 lần** trên khoảng $A=16$–$208$, trong khi $B/A$ chỉ biến thiên khoảng 10%.
- Hai giá trị bằng nhau chỉ tại $A\approx(41/10)^3=69$.
- Khi bật $\gamma_x=0.789$, mọi giá trị trong cột mô hình giảm thêm khoảng một nửa.

## 4. Khung đã sửa

| Mục | Sau kiểm chứng |
|---|---|
| Nguồn $\gamma_\theta$ | Cản bức xạ cho $\tau\sim10^3$–$10^4$, không phải 10–50. Công thức đề xuất sai thứ nguyên |
| Langevin đầy đủ | Tái khóa 100%. Kịp trước khi mất pha: 87–98% |
| Năng lượng xung | 7.56 MeV khi $\gamma_x=0$; **3.7–4.1 MeV** khi $\gamma_x=m/\tau_{\text{relax}}$. Tỉ lệ với $\hbar\omega_0$, không phổ quát theo $A$ |
| Đối chiếu $B/A$, ⁵⁶Fe(n,γ) | Trùng hợp số tại $A\approx69$. Không phải dự đoán độc lập |

## Câu hỏi còn mở

- **$\gamma_\theta$ chưa có nguồn thứ nguyên đúng.** Nếu giữ cản bức xạ ($\tau\sim10^3$–$10^4$), mô-men lái rất yếu và chu kỳ gãy thưa. Học thuyết chấp nhận nhịp gãy chậm như vậy, hay cần một cơ chế cản khác?
- **Dự đoán quan sát được thật sự của mô hình là gì?** Năng lượng xung phụ thuộc $\hbar\omega_0$, $\gamma_x$ và $R/\kappa$. Nếu muốn so thực nghiệm, cần cố định các tham số này từ số liệu độc lập **trước**, rồi mới so phổ. Không nên chọn tham số sao cho khớp một con số.
