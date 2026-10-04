# Phát xạ tiền cân bằng, "bão hòa hình học" và bội số hình học: giải và kiểm chứng số liệu

> Tiếp nối [f-cascade-kiem-chung.md](f-cascade-kiem-chung.md).
> Số liệu tính bằng [giai_so_saturation.py](giai_so_saturation.py) (chạy `python3 giai_so_saturation.py`, ≈ 4.5 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \tau=10,\ \gamma_x=0.70$, $\hbar\omega_0=10$ MeV, $\dot\theta_{\text{gãy}}=3.2048$.

## Tóm tắt kết luận

Đề xuất lần này **đúng về cơ chế**: gãy do gợn tất định cộng một nhiễu nhỏ, và bão hòa khi $k_BT\gtrsim1$ MeV. Mô phỏng xác nhận điều này. **Nhưng các con số định lượng đi kèm thì sai.**

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\tau_{\text{step}}=2R/v_F=10\text{ fm}/0.27c\approx1.2\times10^{-22}$ s | **Đúng số** ($1.235\times10^{-22}$ s) |
| "Trùng khớp hoàn hảo" với khoảng cách xung $1.47\times10^{-22}$ s | **Cùng bậc, lệch 19%.** Tỉ số (khoảng cách xung mô hình)/(2R/v_F) = **1.23 với mọi A** (16, 56, 208). Lý do: cả hai cùng tỉ lệ $A^{1/3}$, vì $T_0\propto1/\hbar\omega_0\propto A^{1/3}$ **theo cách chọn tham số**. Mọi đại lượng tỉ lệ với $T_0$ đều có luật co giãn này, nên sự trùng luật tỉ lệ không phải bằng chứng độc lập |
| Mỗi xung Carnot (4.0–4.6 MeV) là một nucleon phát ra tiền cân bằng, "bứt khỏi bề mặt Fermi" | **Không đủ năng lượng.** Muốn thoát khỏi hạt nhân bền, nucleon cần vượt năng lượng tách $S_n\approx6$–$10$ MeV, lớn hơn năng lượng một xung. Trong mô hình exciton, mỗi bước va chạm tạo một cặp hạt–lỗ **bên trong** hạt nhân; chỉ khi một nucleon tích đủ $>S_n$ mới có phát xạ. Đối chiếu hợp lý hơn: **mỗi xung là một bước exciton nội tại**, không phải một lần phát nucleon |
| Gợn tất định: $\lambda_{\max}^{\text{det}}\approx m\kappa\dot\theta_d^2$, $\Delta\lambda/R\approx1-0.82^2=0.327$ | **Số học đúng, xấp xỉ sai.** Tiệm cận $R\approx m\kappa\dot\theta_{\text{gãy}}^2=3.424$, lệch **+14%** so với $R=3$. Đo trực tiếp (không nhiễu): $M_{\det}=\max\lvert\lambda\rvert/R=$ **0.863** tại $f=1.05$, nên khoảng cách thật tới ngưỡng là **13.7%**, không phải 32.7% |
| $\tau_{\text{lưu}}$ do khoảng cách hình học $\Delta\lambda(f)$ quyết định | ✓ về chiều. Nhưng tại **$f=1.2$: $M_{\det}=1.10>1$**, tức gợn tất định **đã vượt ngưỡng**, nên ở đó gãy là tất định chứ không cần nhiễu. Biên $M_{\det}=1$ nằm ở $f\approx1.14$ |
| Bão hòa khi $k_BT\ge1$ MeV | **✓ Được xác nhận.** Giảm nhiễu: 0.01 và 0.03 MeV cho **0** lần gãy; 0.1 MeV cho tốc độ **0.006**/$T_0$; 0.3 MeV cho **0.100**/$T_0$; 1 MeV cho **0.152**/$T_0$. Từ 1 tới 8 MeV (đo trước): 0.14–0.16/$T_0$. Tốc độ **có** phụ thuộc nhiễu dưới ~1 MeV và bão hòa từ ~1 MeV trở lên |
| Hàng rào hiệu dụng | Từ hai điểm 0.1 và 0.3 MeV (ước lượng thô): $\Delta_{\text{eff}}\approx0.4$ MeV. Nhỏ hơn **40 lần** so với $\Delta K=16.75$ MeV của Kramers. Nhiễu cỡ 1 MeV vì vậy đã vượt hẳn hàng rào, giải thích được hiện tượng bão hòa |
| $p_{\text{trap}}\approx0.20$ là hằng số hình học; $\langle N\rangle=1/p=5.13$ "trùng tuyệt đối" với 5.12 | **Vòng tròn.** $p=0.195$ được **tính từ** $\langle N\rangle=5.12$ ($p=1/5.12$), nên $1/p=5.13$ là đồng nhất thức. Phần kiểm chứng thật sự là **dạng** phân bố (hình học, $\chi^2=10.9$), đã có ở tài liệu trước |
| $\tau_{\text{lưu}}\approx4.5\,T_0\approx1.8\times10^{-21}$ s | ✓ ($4.57\,T_0=1.89\times10^{-21}$ s) |

## 1. Thang thời gian

| A | Khoảng cách xung mô hình ($0.355\,T_0$) | $2R/v_F$ ($R=1.2A^{1/3}$ fm, $v_F=0.275c$) | Tỉ số |
|---|---|---|---|
| 16 | $9.02\times10^{-23}$ s | $7.34\times10^{-23}$ s | 1.23 |
| 56 | $1.37\times10^{-22}$ s | $1.11\times10^{-22}$ s | 1.23 |
| 208 | $2.12\times10^{-22}$ s | $1.72\times10^{-22}$ s | 1.23 |

## 2. Lề ứng suất tất định (không nhiễu)

| $f$ | $M_{\det}=\max\lvert\lambda\rvert/R$ (đo) | Đề xuất $(\dot\theta_d/\dot\theta_{\text{gãy}})^2$ | $\dot\theta$ min–max | $\dot\theta$ TB |
|---|---|---|---|---|
| 0.95 | 0.7236 | 0.5515 | 1.870–2.673 | 2.304 |
| 1.00 | 0.7913 | 0.6111 | 2.024–2.785 | 2.446 |
| **1.05** | **0.8628** | **0.6738** | 2.173–2.898 | 2.575 |
| 1.10 | 0.9381 | 0.7394 | 2.319–3.012 | 2.707 |
| 1.20 | **1.1006** (vượt ngưỡng) | 0.8800 | 2.602–3.243 | 2.974 |

**Đọc bảng:**
- Công thức đề xuất đánh giá thấp lề ứng suất khoảng 0.17–0.22 ở mọi $f$. Nguyên nhân: nó bỏ qua gợn vận tốc (θ̇ dao động ±0.35 quanh giá trị trung bình) và dùng dạng tiệm cận của $R$.
- Theo lề ứng suất đo được:

| Vùng | Cơ chế gãy |
|---|---|
| $f\lesssim1.14$ | Cần nhiễu |
| $f\gtrsim1.14$ | Tất định, không cần nhiễu |

## 3. Tốc độ gãy khi giảm nhiễu ($f=1.05$, 8 quỹ đạo × $60T_0$)

| $k_BT$ | Tổng số gãy | Số pha dài | Tốc độ đợt ($/T_0$) | $\tau_{\text{lưu}}$ ($T_0$) |
|---|---|---|---|---|
| 0.01 MeV | 0 | 0 | 0 | — |
| 0.03 MeV | 0 | 0 | 0 | — |
| 0.1 MeV | 11 | 3 | 0.0063 | 14.0 |
| 0.3 MeV | 112 | 48 | 0.100 | 7.59 |
| 1 MeV | 385 | 73 | 0.152 | 4.63 |
| 2–8 MeV (tài liệu trước) | — | — | 0.140–0.158 | — |

Đường cong có hai vùng rõ rệt:
- **Vùng kích hoạt** ($k_BT\lesssim0.3$ MeV): tốc độ tăng nhanh theo nhiễu.
- **Vùng bão hòa** ($k_BT\gtrsim1$ MeV): tốc độ gần như phẳng.

Điều này khớp với hình ảnh "gợn tất định + xung nhiễu nhỏ": hàng rào hiệu dụng chỉ ~0.4 MeV, nên nhiễu cỡ 1 MeV đã vượt qua gần như ở mỗi vòng quay có cơ hội.

## 4. Bảng chuẩn đã sửa

| Thành phần | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Thời gian lưu | Mũ, bão hòa hình học theo $f$ | ✓ Cơ chế đúng. Hàng rào hiệu dụng ~0.4 MeV, bão hòa từ ~1 MeV. Khoảng cách tới ngưỡng tại $f=1.05$ là **13.7%**, không phải 32.7%. Với $f\gtrsim1.14$, gãy là tất định |
| Đợt xả vi chấn | Hình học, $\langle N\rangle\approx5$ | ✓ Dạng hình học đúng. "$1/p=5.13$ khớp 5.12" là đồng nhất thức, không phải kiểm chứng |
| Bản chất phát xạ | Phát nucleon tiền cân bằng | ✗ khi hiểu là phát nucleon ra ngoài, vì 4–4.6 MeV < $S_n$. Có thể hiểu là **bước exciton nội tại**. Thang thời gian cùng bậc với $2R/v_F$ (tỉ số 1.23), nhưng điều này có sẵn do cách chọn $T_0$ |
| Năng lượng xung | $\approx0.4\,\hbar\omega_0$ | ✓ |

## Câu hỏi còn mở

- **Kiểm "bước exciton" bằng số.** Nếu mỗi xung là một bước exciton, thì số bước trước khi phát nucleon phải đủ để tích lũy vượt $S_n\approx8$ MeV, tức cần khoảng 2 xung trở lên cộng lại. Mô hình hiện tại chưa có biến theo dõi năng lượng tích lũy của từng nucleon. Cần thêm một biến như vậy nếu muốn so với phổ nucleon tiền cân bằng thực nghiệm.
- **Vùng $f\gtrsim1.14$** (gãy tất định) có ứng với trạng thái vật lý nào không? Nếu không, bảng tham số nên ghi rõ ràng buộc $f<1.14$.
