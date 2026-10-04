# Mô-men lái Van der Pol và nhiễu Langevin: giải và kiểm chứng số liệu

> Tiếp nối [giu-pha-rayleigh-kiem-chung.md](giu-pha-rayleigh-kiem-chung.md).
> Số liệu tính bằng [giai_so_dong_co_langevin.py](giai_so_dong_co_langevin.py).
> - `A`: tất định, ≈ 75 s.
> - `B`, `C`: Monte Carlo 64 quỹ đạo, chạy song song 8 lõi, ≈ 2 phút.
> - `Q`: nhiễu điểm không 5 MeV.
>
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0$. Một đơn vị năng lượng bằng $\hbar\omega_0=10$ MeV, $T_0=2\pi$, $\dot\theta_{\text{gãy}}=3.20479$.
>
> **Tích phân:** RK4 cho phần tất định, cộng một cú hích Euler–Maruyama sau mỗi bước. Trong pha liên kết, cú hích chia cho quán tính hiệu dụng $I+m\kappa^2\cos^2\theta$ và giữ đúng $\dot x=\kappa\cos\theta\,\dot\theta$.

## Tóm tắt kết luận

**Mô-men lái Van der Pol giải được bài toán "động cơ tắt":** chu trình gãy → phát xạ → tái khóa lặp lại **tuần hoàn**, và năng lượng cân bằng đúng.

**Nhiễu Langevin:**
- Hạt tái khóa **64/64** ở mọi mức nhiễu đã thử, trung bình sau $\approx0.4\,T_0$.
- Công thức $\tau_{\text{dephasing}}$ của đề xuất sai so với đo trực tiếp (2–40 lần).
- Con số "5–10 $T_0$" không đạt được: đo được 1.1–2.1 $T_0$ với nhiễu nhiệt 1 MeV, và chỉ còn 0.3–0.6 $T_0$ khi tính nhiễu điểm không lượng tử.

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| Mô-men lái ròng $Q_\theta=\gamma_\theta(1-\dot\theta^2/\Omega_0^2)\dot\theta$ | **Đúng dạng** Rayleigh–Van der Pol: bơm năng lượng khi $\lvert\dot\theta\rvert<\Omega_0$, hãm khi $\lvert\dot\theta\rvert>\Omega_0$ |
| Trạng thái dừng $\dot\theta_{\text{ss}}=\Omega_0$ "ổn định tuyệt đối" | **Xấp xỉ.** Thế $-\Delta E\cos\theta$ vẫn làm $\dot\theta$ gợn quanh $\Omega_0$ (ví dụ $\dot\theta=2.95$ khi $\Omega_0=2.88$). Với lái yếu ($\tau=50$) xuất phát chậm ($\dot\theta_0=0.5$), sau $20T_0$ động cơ **vẫn chưa vượt rào** ($\dot\theta=-0.65$) |
| Chu trình gãy → phát xạ → tái khóa tuần hoàn, động cơ không tắt | **Đúng** (mô phỏng tới $50T_0$). $\tau=10$, $f=1.05$: **tuần hoàn chính xác**, chu kỳ **4.282 = 0.681 $T_0$**, 75 lần gãy. $f=1.1$: tuần hoàn bậc 3 (4.184 / 2.965 / 3.366). $\tau=50$: gần tuần hoàn, 2.5–3.0 $T_0$ |
| Carnot kéo $\dot\theta$ xuống dưới ngưỡng sau mỗi lần khóa | **Đúng khi $f=1.05$** ($\dot\theta$ sau khóa tối đa 3.2042 < 3.2048, sát ngưỡng). **Không đúng khi $f=1.1$**, $\tau=10$: có lần $\dot\theta$ sau khóa là 3.365 > ngưỡng, nên gãy lại ngay |
| Nguồn năng lượng | **Cân bằng đúng:** công lái − tổn hao Carnot = biến thiên năng lượng hệ (ví dụ 56.0652 − 56.5645 = −0.4993, so với ΔE hệ −0.5009). **Mỗi lần khóa phát xạ ≈ 0.76–0.79 đơn vị = 7.6–7.9 MeV** ở cả 4 cấu hình. Năng lượng đến từ một nguồn lái **bên ngoài** mô hình Lagrange (lực suy rộng không bảo toàn), không phải chuyển động vĩnh cửu |
| "Hàng triệu chu kỳ" | **Chưa kiểm:** mới mô phỏng tới $50T_0$ (75 chu kỳ) |
| Phương trình Langevin, $\langle\xi\xi\rangle=2D_\theta\delta$, $D_\theta=\gamma_\theta k_BT$ | **Đúng dạng** (định lý thăng giáng–tiêu tán cổ điển) |
| Chuỗi biến đổi ra $\tau_{\text{dephasing}}$ | **Không nhất quán ở bước giữa:** $1/(2D_{\text{diff}})=\gamma_\theta/(2k_BT)$, không phải $I/(2\gamma_\theta k_BT)$. Kết quả cuối $I/(2k_BT\tau)$ đúng đại số, nhưng chỉ áp cho **rô-to tự do không thế, không lái**, và đó không phải hệ đang xét |
| "$I=1$, $k_BT=0.1$ cho $\tau_{\text{deph}}\sim5$–$10\,T_0$" | **Không tái lập được.** Theo công thức, $\tau_{\text{deph}}=5/\tau$, nên muốn có 5–10 $T_0$ cần $\tau=0.08$–$0.16$ (cản cực mạnh), mà đề xuất không cho $\tau$. **Đo trực tiếp:** $\tau=1$ cho **2.07 $T_0$** (công thức: 0.80); $\tau=10$ cho **1.11 $T_0$** (công thức: 0.08) |
| $k_BT_{\text{nuc}}\approx1$ MeV làm nhiễu | **Thiếu nhiễu lượng tử.** Vì $\hbar\omega_0=10$ MeV $\gg k_BT$, nhiễu điểm không cỡ $\hbar\omega_0/2=5$ MeV **lấn át** nhiễu nhiệt; Langevin cổ điển với 1 MeV đánh giá thấp nhiễu khoảng 5 lần. Với 5 MeV: $\tau_{\text{deph}}$ đo được còn **0.32 $T_0$** ($\tau=1$) và **0.64 $T_0$** ($\tau=10$) |
| Tái khóa trước khi mất pha | **Đúng với nhiễu 1 MeV:** 64/64 (cả $\tau=1$ và 10). **Với nhiễu 5 MeV:** 64/64 khi $\tau=10$, nhưng chỉ **16/64** khi $\tau=1$ |

## 1. Mô-men lái Van der Pol, không nhiễu ($50T_0$)

| $\tau=I/\gamma_\theta$ | $f=\Omega_0/\dot\theta_{\text{gãy}}$ | Gãy / khóa | Khoảng giữa hai lần gãy (5 lần cuối) | $\dot\theta$ sau khóa (min–max) | Công lái | Carnot | Chênh | ΔE hệ |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.05 | 75 / 74 | 4.282 ×5 (**tuần hoàn**, 0.681 $T_0$) | 2.3874–3.2042 | 56.0652 | 56.5645 | −0.4993 | −0.5009 |
| 10 | 1.10 | 91 / 90 | 4.184, 2.965, 3.366, 4.184, 2.965 (**tuần hoàn bậc 3**) | 2.4169–**3.3650** | 70.5427 | 71.3668 | −0.8242 | −0.8269 |
| 50.27 | 1.05 | 18 / 18 | 17.587, 18.438, 21.143, 17.587, 18.438 (≈2.97 $T_0$) | 2.3473–3.1863 | 12.7738 | 13.7395 | −0.9657 | −0.9658 |
| 50.27 | 1.10 | 22 / 22 | 16.720, 15.293, 14.014, 16.718, 15.292 (≈2.48 $T_0$) | 2.3542–3.3188 | 16.0226 | 17.2541 | −1.2315 | −1.2316 |

**Kiểm tra chu trình giới hạn** (không gãy, $\Omega_0=0.9\,\dot\theta_{\text{gãy}}$, xuất phát $\dot\theta=0.5$):
- $\tau=10$: $\dot\theta(20T_0)=2.95$, đã về quanh $\Omega_0$.
- $\tau=50$: $\dot\theta(20T_0)=-0.65$, vẫn còn lắc, chưa quay.

**Năng lượng phát xạ mỗi chu kỳ** (Carnot / số lần khóa):

| Cấu hình | Mỗi lần khóa |
|---|---|
| $\tau=10$, $f=1.05$ | 0.764 = 7.64 MeV |
| $\tau=10$, $f=1.1$ | 0.793 = 7.93 MeV |
| $\tau=50$, $f=1.05$ | 0.763 = 7.63 MeV |
| $\tau=50$, $f=1.1$ | 0.784 = 7.84 MeV |

Năng lượng này **gần như không đổi** theo $\tau$ và $f$: cỡ **7.6–7.9 MeV mỗi xung**. Đây là một dự đoán định lượng của mô hình.

## 2. Khuếch tán pha dưới nhiễu Langevin

**Thiết lập:** không cho gãy, $\Omega_0=0.9\,\dot\theta_{\text{gãy}}$, $N=64$. So pha của mỗi quỹ đạo có nhiễu với quỹ đạo không nhiễu cùng điều kiện đầu. $\tau_{\text{deph}}$ là thời điểm $\langle\Delta\theta^2\rangle=1$ (độ phân giải lấy mẫu 0.5).

| $k_BT$ | $\tau$ ($\gamma_\theta$) | $\langle\Delta\theta^2\rangle(10T_0)$ | $\tau_{\text{deph}}$ **đo** | Đề xuất $I/(2k_BT\tau)$ | Rô-to trên chu trình giới hạn, $2I/(k_BT\tau)$ |
|---|---|---|---|---|---|
| 0.1 (1 MeV) | 1 (1.0) | 8.24 | **13.0 = 2.07 $T_0$** | 5.00 = 0.80 $T_0$ | 20.0 = 3.18 $T_0$ |
| 0.1 (1 MeV) | 10 (0.1) | 38.3 | **7.0 = 1.11 $T_0$** | 0.50 = 0.08 $T_0$ | 2.00 = 0.32 $T_0$ |
| 0.5 (5 MeV, điểm không) | 1 (1.0) | 11 466 | **2.0 = 0.32 $T_0$** | 1.00 = 0.16 $T_0$ | 4.00 = 0.64 $T_0$ |
| 0.5 (5 MeV, điểm không) | 10 (0.1) | 3 367 | **4.0 = 0.64 $T_0$** | 0.10 = 0.016 $T_0$ | 0.40 = 0.064 $T_0$ |

**Nhận xét:**
- **Cả hai công thức giải tích đều lệch so với số đo.** Nguyên nhân: rô-to trong thế $-\Delta E\cos\theta$ **không đẳng thời**. Chu kỳ quay phụ thuộc năng lượng, nên nhiễu năng lượng chuyển thành trôi tần số, và $\langle\Delta\theta^2\rangle$ tăng nhanh hơn tuyến tính.
- Ở nhiễu 5 MeV, $\langle\Delta\theta^2\rangle$ lên tới hàng nghìn: một số quỹ đạo bị đẩy **đảo chiều quay**.

## 3. Monte Carlo tái khóa (lái + nhiễu, $\Omega_0=1.05\,\dot\theta_{\text{gãy}}$, $N=64$, cửa sổ $10T_0$)

| $k_BT$ | $\tau$ | Gãy ≥1 lần | Tái khóa sau lần gãy đầu | Thời gian tái khóa: min / TB / max | Số lần gãy TB mỗi quỹ đạo | Kịp trước $\tau_{\text{deph}}$ **đo** | Kịp trước $\tau_{\text{deph}}$ **đề xuất** |
|---|---|---|---|---|---|---|---|
| 1 MeV | 1 | 64/64 | **64/64** | 1.616 / 2.704 / 3.252 (TB 0.430 $T_0$) | 28.0 | **64/64** (≤ 13.0) | 64/64 (≤ 5.0) |
| 1 MeV | 10 | 64/64 | **64/64** | 1.636 / 2.550 / 3.242 (TB 0.406 $T_0$) | 16.7 | **64/64** (≤ 7.0) | 0/64 (≤ 0.5) |
| 5 MeV | 1 | 64/64 | **64/64** | 0.003 / 2.498 / 3.302 (TB 0.398 $T_0$) | 24.6 | **16/64** (≤ 2.0) | 5/64 (≤ 1.0) |
| 5 MeV | 10 | 64/64 | **64/64** | 1.431 / 2.555 / 3.254 (TB 0.407 $T_0$) | 15.1 | **64/64** (≤ 4.0) | 0/64 (≤ 0.1) |

**Đọc kết quả:**
1. **Tái khóa rất bền:** 100% quỹ đạo tái khóa, trung bình sau 0.40–0.43 $T_0$, gần như không phụ thuộc nhiễu.
2. **Kịp trước khi mất pha (dùng $\tau_{\text{deph}}$ đo trực tiếp):** 64/64 ở 3/4 cấu hình. Cấu hình duy nhất thất bại phần lớn là nhiễu điểm không 5 MeV kèm cản mạnh ($\tau=1$): chỉ 16/64 (25%).
3. **Nếu dùng công thức của đề xuất**, kết luận đổi hẳn (0/64 ở ba cấu hình). Vì vậy phải dùng $\tau_{\text{deph}}$ đo, không dùng công thức rô-to tự do.
4. **Lưu ý:** mô phỏng chỉ đặt nhiễu lên pha $\theta$. Va chạm nucleon–nucleon tác động lên hạt $x$ ($\tau_{\text{relax}}=0.2\,T_0$, nhỏ hơn thời gian tái khóa 0.4 $T_0$) **chưa có** trong mô hình.

## 4. Khung đã sửa

| Mục | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Nguồn lái | Van der Pol, $\dot\theta_{\text{ss}}=\Omega_0$ | ✓ Chu trình tuần hoàn (0.68 $T_0$ khi $\tau=10$, $f=1.05$). Nên chọn $f\approx1.05$: $f=1.1$ gây gãy lại ngay sau khóa |
| Phát xạ mỗi chu kỳ | $\tfrac12m_{\text{eff}}u^2$ | ✓ **≈ 7.6–7.9 MeV/xung**, ổn định theo $\tau$ và $f$ |
| $\tau_{\text{dephasing}}$ | $I/(2k_BT\tau)$, "5–10 $T_0$" | ✗ Đo được 1.1–2.1 $T_0$ (nhiễu 1 MeV) và 0.3–0.6 $T_0$ (nhiễu 5 MeV) |
| Tái khóa trước mất pha | "đủ dài" | ✓ 64/64 với nhiễu nhiệt 1 MeV. Với nhiễu điểm không 5 MeV: ✓ khi $\tau=10$, ✗ (25%) khi $\tau=1$ |

## Câu hỏi còn mở

- **$\gamma_\theta$ (tức $\tau$) chưa được xác định.** Kết quả với nhiễu 5 MeV phụ thuộc mạnh vào nó: 25% khi $\tau=1$, 100% khi $\tau=10$. Học thuyết cần một nguồn độc lập cho độ cản của pha.
- **Nhiễu lên hạt $x$ chưa có trong mô hình.** Va chạm nucleon, với $\tau_{\text{relax}}=0.2\,T_0$, cần một lực Langevin thứ hai trong phương trình (2), kèm $\gamma_x$ tương ứng. Đây là bước kiểm chứng tiếp theo.
- **Năng lượng xung 7.6–7.9 MeV** là một dự đoán cụ thể. Có số liệu thực nghiệm nào để so sánh, ví dụ phổ gamma của hạt nhân ở vùng này?
