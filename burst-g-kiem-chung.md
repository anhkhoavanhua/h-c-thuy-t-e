# Thời gian lưu do ghép cặp $g$ và "biến cố E gãy = đợt vi chấn ≈ GDR": giải và kiểm chứng số liệu

> Tiếp nối [rice-poisson-kiem-chung.md](rice-poisson-kiem-chung.md).
> Số liệu tính bằng [giai_so_burst_g.py](giai_so_burst_g.py) (chạy `python3 giai_so_burst_g.py`, ≈ 6.5 phút trên 8 lõi).
> Sự kiện khóa trong [giai_so_dong_co_langevin.py](giai_so_dong_co_langevin.py) nay ghi thêm năng lượng Carnot của từng lần khóa.
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ e=0,\ \tau=10,\ \gamma_x=0.70,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}(g)$, $\hbar\omega_0=10$ MeV.
> **Định nghĩa:** "Pha đồng bộ dài" là lần khóa kéo dài ≥ $T_0$. "Đợt vi chấn" là chuỗi gãy/khóa nằm giữa hai pha đồng bộ dài.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $V_{\text{eff}}=-\Delta E\cos\theta+\tfrac12(k\kappa^2-g\kappa)\sin^2\theta$ | **Sai hệ số của $g$.** Rút gọn từ $L$ (xem [lien-ket-dong-bo-e-gay-kiem-chung.md](lien-ket-dong-bo-e-gay-kiem-chung.md), mục 2.1) cho $V_{\text{eff}}=-\Delta E\cos\theta+(\tfrac12k\kappa^2-g\kappa)\sin^2\theta$, tức hệ số của $g\kappa$ là 1 chứ không phải ½ |
| $g\kappa\sin2\theta$ là nhiễu loạn tham số làm $\dot\theta$ dao động với chu kỳ π | **Đúng** như một thành phần tất định. Nhưng xem dòng dưới |
| Trùng pha giữa $\ddot\theta$ và $-g\sin\theta$ **kích hoạt** $\lvert\lambda\rvert\ge R$ | **Bị bác bỏ.** Mô phỏng **không nhiễu** với $g=0$, 0.1, 0.2, 0.3: **0 lần gãy** trong $50T_0$ ở mọi $g$. Riêng số hạng $g$ không kích hoạt được gãy |
| "Hỗn loạn yếu, hồi quy Poincaré ⟹ Poisson" | **Không áp dụng được.** Pha đồng bộ là hệ tất định **1 bậc tự do, tự trị**, nên không thể hỗn loạn (cần ít nhất 1.5 bậc tự do). Hồi quy Poincaré cũng **không** kéo theo phân bố Poisson. Dạng mũ đo được đến từ **nhiễu Langevin** |
| $\tau_{\text{lưu}}\propto1/g$; $g\to0$ ⟹ $\tau_{\text{lưu}}\to\infty$ | **Bị bác bỏ.** Với **$g=0$: $\tau_{\text{lưu}}=4.51\,T_0$ hữu hạn** (74 pha dài, 425 lần gãy). Theo $g=0$ / 0.05 / 0.1 / 0.2: $\tau_{\text{lưu}}=4.51$ / 5.34 / 6.09 / 5.06 $T_0$, gần như phẳng (sai số thống kê ~15%). $g\cdot\tau_{\text{lưu}}$ thay đổi 0 → 1.01, không phải hằng số |
| Thời gian lưu dạng mũ | ✓ CV phần dư 0.95–1.25 ở mọi $g$ |
| Đợt vi chấn gồm $N_{\text{sub}}\approx5$–$10$ xung, $\Delta t_{\text{burst}}\approx0.5$–$1.0\,T_0$ | **Phân bố rộng hơn nhiều.** 1 MeV: $N$ TB **4.54**, trung vị 3, khoảng 10–90% là **1–9**; thời lượng TB **1.20 $T_0$**. 5 MeV: $N$ TB **12.0**, trung vị 8, khoảng 10–90% là **2–27**; thời lượng **2.58 $T_0$** |
| $\Delta E_{\text{biến cố}}\approx5.2\times4.59=23.8$ MeV | **Số học đúng, nhưng ghép hai nguồn không cùng điều kiện.** 4.59 MeV là giá trị của ⁵⁶Fe ($\hbar\omega_0=10.72$), còn 5.2 là số đếm ở $\hbar\omega_0=10$. **Đo trực tiếp:** tổng Carnot mỗi đợt ở 1 MeV có TB **17.96**, trung vị **12.72**, khoảng 10–90% **5.3–36.4 MeV**. Ở 5 MeV: TB **55.0**, trung vị 43.2, khoảng 10–90% **8.5–139 MeV** |
| ≈ đỉnh GDR 15–25 MeV | **Không tương thích về bản chất:** (1) năng lượng mỗi đợt **phụ thuộc mạnh vào nhiễu** (×3 từ 1 lên 5 MeV) và **trải rộng 7 lần**, trong khi GDR là một đỉnh cộng hưởng xác định; (2) một đợt gồm trung bình 4.5 xung **riêng rẽ**, mỗi xung ~4 MeV, nên phổ phát xạ là các xung ~4 MeV, không phải một lượng tử 15–25 MeV như GDR; (3) $80A^{-1/3}$ chỉ đúng cho hạt nhân nặng (A=16: 31.7 so với thực nghiệm ~22–24 MeV) |
| Thời gian xả ~$10^{-21}$ s "tương thích GDR" | Thời lượng đợt đo được là $5.0\times10^{-22}$ s (1 MeV) và $1.07\times10^{-21}$ s (5 MeV). Thời gian sống GDR là $\hbar/\Gamma\approx$ **$1.5\times10^{-22}$ s** ($\Gamma\approx4.5$ MeV), ngắn hơn **3–7 lần** |

## 1. Tất định: $g$ không kích hoạt gãy

| $g$ | $\dot\theta_{\text{gãy}}$ | Số lần gãy (không nhiễu, $50T_0$) |
|---|---|---|
| 0.0 | 3.2554 | 0 |
| 0.1 | 3.2048 | 0 |
| 0.2 | 3.1533 | 0 |
| 0.3 | 3.1010 | 0 |

## 2. Thời gian lưu theo $g$ (nhiễu 1 MeV, 8 quỹ đạo × $60T_0$)

| $g$ | $\dot\theta_{\text{gãy}}$ | Tổng số gãy | Số pha dài | $\tau_{\text{lưu}}$ TB ($T_0$) | Phần dư TB ($T_0$) | CV phần dư | $g\cdot\tau_{\text{lưu}}$ |
|---|---|---|---|---|---|---|---|
| **0** | 3.2554 | 425 | 74 | **4.509** | 3.509 | 1.171 | 0 |
| 0.05 | 3.2302 | 404 | 66 | 5.338 | 4.338 | 1.252 | 0.267 |
| 0.1 | 3.2048 | 315 | 64 | 6.088 | 5.088 | 1.187 | 0.609 |
| 0.2 | 3.1533 | 375 | 72 | 5.056 | 4.056 | 0.953 | 1.011 |

Nếu $\tau_{\text{lưu}}\propto1/g$ thì $\tau_{\text{lưu}}(0.05)/\tau_{\text{lưu}}(0.2)$ phải bằng 4. Đo được **1.06**. Thời gian lưu gần như không phụ thuộc $g$, và vẫn hữu hạn khi $g=0$.

## 3. Cấu trúc đợt vi chấn ($g=0.1$, 16 quỹ đạo × $60T_0$)

| Chỉ số | 1 MeV | 5 MeV |
|---|---|---|
| Số đợt | 119 | 113 |
| Số lần khóa mỗi đợt: TB / trung vị / 10–90% | 4.54 / 3 / 1–9 | 12.01 / 8 / 2–27 |
| Thời lượng TB | 1.20 $T_0$ = $4.96\times10^{-22}$ s | 2.58 $T_0$ = $1.07\times10^{-21}$ s |
| **Tổng Carnot mỗi đợt**: TB / trung vị / 10–90% | **17.96 / 12.72 / 5.31–36.42 MeV** | **55.01 / 43.22 / 8.49–139.33 MeV** |
| Carnot mỗi xung trong đợt (TB) | 3.96 MeV | 4.58 MeV |
| Xung lớn nhất mỗi đợt (TB) | 8.07 MeV | 19.70 MeV |

**GDR thực nghiệm để so sánh:**
- Thang $31.2A^{-1/3}+20.6A^{-1/6}$: 18.7 MeV (A=56), 13.7 MeV (A=208).
- Dạng $80A^{-1/3}$: 31.7 / 20.9 / 13.5 MeV (A=16 / 56 / 208).

**Đọc kết quả:**
- Tổng năng lượng TB ở 1 MeV (18 MeV) **tình cờ** gần GDR của ⁵⁶Fe. Nhưng phân bố năng lượng các đợt trải từ 5 tới 36 MeV, và tăng gấp 3 khi nhiễu lên 5 MeV. Một đại lượng phụ thuộc mạnh vào nhiệt độ và trải rộng như vậy không thể đồng nhất với một đỉnh cộng hưởng.
- Ngoài ra, đơn vị phát xạ trong mô hình là **từng xung Carnot** (~4 MeV). Tổng của nhiều xung riêng rẽ không tạo ra một lượng tử năng lượng lớn.

## 4. Bảng kiến trúc đã sửa

| Khía cạnh | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Cơ chế nhịp gãy | "Thời gian lưu Poisson do nhiễu loạn ghép cặp $g$, $\tau_{\text{lưu}}\propto1/g$" | Thời gian lưu dạng mũ ✓, **do nhiễu Langevin**, gần như độc lập $g$ (4.5–6.1 $T_0$ với $g=0$–0.2). $g$ không kích hoạt gãy tất định |
| Biến cố E gãy | Đợt 5–10 xung, 0.5–1 $T_0$ | Đợt có $N$ trung vị 3 (1 MeV) và 8 (5 MeV), phân bố rộng (1–27); thời lượng 1.2–2.6 $T_0$ |
| Năng lượng xung đơn | $(0.42$–$0.48)\hbar\omega_0\approx4$ MeV | ✓ (TB 3.96–4.58 MeV trong đợt) |
| Năng lượng toàn đợt | 20–25 MeV, thang GDR | TB 18 MeV (1 MeV) / 55 MeV (5 MeV), trải 5–139 MeV. **Không tương thích GDR** (không phải một đỉnh, phụ thuộc nhiệt độ, là tổng của nhiều xung nhỏ) |

## Câu hỏi còn mở

- **Nguồn thật của thời gian lưu ~4.5–6 $T_0$** vẫn chưa xác định. Đã loại trừ ba nguồn: nhiệt độ ($k_BT$), hệ số ghép $g$, và nạp tất định. Còn lại các tham số $\tau$, $\gamma_x$ và $\Omega_0/\dot\theta_{\text{gãy}}$ (tỉ số đang đặt là 1.05). Bước kiểm tiếp theo: đo $\tau_{\text{lưu}}$ theo $\Omega_0/\dot\theta_{\text{gãy}}$.
- **Muốn so với phổ thực nghiệm**, mô hình cần một cơ chế để một biến cố phát ra **một** lượng tử (tổng năng lượng), thay vì nhiều xung riêng rẽ. Mô hình hiện tại không có cơ chế này.
