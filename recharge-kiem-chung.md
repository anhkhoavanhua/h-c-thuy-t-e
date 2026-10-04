# Chu kỳ nạp Van der Pol, tiên đề tự đồng dạng và "cận trên Carnot": giải và kiểm chứng số liệu

> Tiếp nối [kramers-A-kiem-chung.md](kramers-A-kiem-chung.md).
> Số liệu tính bằng [giai_so_recharge.py](giai_so_recharge.py) (chạy `python3 giai_so_recharge.py`, ≈ 4 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}$ ($\dot\theta_{\text{gãy}}=3.2048$), nhiễu 1 MeV. $T_0=2\pi$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| Phương trình trung bình $I_{\text{eff}}\dot\omega=\gamma_\theta(1-\omega^2/\Omega_0^2)\omega-\tfrac12\gamma_x\kappa^2\omega$, $I_{\text{eff}}=1.056$ | **Đúng dạng.** Đây là cân bằng công suất ở tài liệu trước |
| $T_{\text{recharge}}=\int_{\dot\theta_d}^{\dot\theta_{\text{gãy}}}d\dot\theta/\ddot\theta$ | **Tích phân không tồn tại.** $\dot\theta_d$ là điểm dừng **hút**. Giữa $\dot\theta_d$ và $\dot\theta_{\text{gãy}}$, vế phải **âm** ($F=-0.041$ tại $\dot\theta=2.918$), nên quỹ đạo tất định **không bao giờ** tự nạp lên tới ngưỡng. Điều này khớp với kết quả 0 lần gãy khi không nhiễu |
| $T_{\text{recharge}}\approx\tau_{\text{eff}}\ln\frac{\Omega_0-\dot\theta_d}{\Omega_0-\dot\theta_{\text{gãy}}}\approx6$–$7$ đv | $\ln(0.735/0.160)=1.523$, nhưng **$\tau_{\text{eff}}$ không được cho**. Muốn ra 6–7 đv thì cần $\tau_{\text{eff}}=3.9$–$4.6$, giá trị bị ngầm chọn |
| $T_{\text{cycle}}=6$–$7$ đv $+\,4$–$5$ đv $\approx$ "6.5–7.0 đv ≈ **6.7 $T_0$**" | **Lỗi đơn vị.** $6$–$7+4$–$5=$ **10–12 đv = 1.6–1.9 $T_0$**, suy ra $r=$ **0.52–0.63 mỗi $T_0$**, gấp 3.5–4 lần số đo. Còn $6.7\,T_0=42.1$ đv. "Khớp 0.149" có được là nhờ đổi đơn vị sai |
| $\dot\theta$ sau khóa $\approx\dot\theta_d\approx2.63$ | **Không đúng.** Đo được **2.15** (với $\tau=10$, $\gamma_x=0.7$) |
| $r\approx0.15$/$T_0$ là chu kỳ nạp của dao động tích thoát | **Không phải hằng số của mô hình.** $r$ thay đổi **8 lần** theo $(\tau,\gamma_x)$: từ 0.029 đến 0.231 mỗi $T_0$. Tỉ số $(1/r)/\tau_{\text{rel}}$ đổi từ 3.7 tới 21, nên không phải một thời gian hồi phục đơn giản. Hệ **chờ trong pha đồng bộ** trung bình 1.6–10.8 $T_0$ trước mỗi đợt gãy, chứ không phải "nạp 1 $T_0$ + ổn định 0.7 $T_0$" |
| Tiên đề tự đồng dạng $E_{\text{vt}}/\hbar\omega_0=3$ | **Nhất quán nội tại.** $E_{\text{vt}}=48.8$ / 36.0 / 20.8 MeV cho O / Ca / Pb là đúng số; $L\cdot\hbar\omega_0=\hbar c/3=65.8$ MeV·fm là đúng. $L/R_{\text{hạt nhân}}=1.34$ không đổi. Đây là **tiên đề mới**, thay cho "$E_{\text{vt}}=30$ MeV (độ sâu Woods–Saxon)" của vòng trước |
| $\Delta E_{\text{Carnot}}\approx(0.42$–$0.48)\hbar\omega_0\approx(17$–$20)A^{-1/3}$ MeV | **Số học đúng** (17.2–19.7) |
| Là "**cận trên** bao trùm các mức 2⁺, 3⁻ thấp nhất" | **Sai với ²⁰⁸Pb.** Mức 2⁺ thấp nhất là **4.085 MeV**, lớn hơn vùng Carnot 2.91–3.32. Với O và Ca thì bao trùm được. Ngoài ra, "định lý Carnot" ở đây là định lý mất động năng khi va chạm, **không phải** cận nhiệt động Carnot |
| Bảng cuối: $E_{\text{quay}}=51.35$ MeV | Chỉ đúng khi $\hbar\omega_0=10$ MeV. Với tiên đề tự đồng dạng: $E_{\text{quay}}=5.135\,\hbar\omega_0(A)$, tức 83.6 (O), 61.6 (Ca), 35.5 (Pb) MeV |
| Bảng cuối: "$r\approx0.15$/$T_0$, $T\approx6.7T_0$" | Chỉ đúng tại **đúng một** điểm $(\tau=10,\gamma_x=0.7)$. Không phải hệ quả tổng quát |

## 1. Phương trình trung bình

| $\tau$ | $\gamma_x$ | $\Omega_0$ | $\dot\theta_d$ | $F(\tfrac12(\dot\theta_d+\dot\theta_{\text{gãy}}))$ | $F'(\dot\theta_d)$ | $\tau_{\text{rel}}=I_{\text{eff}}/\lvert F'\rvert$ |
|---|---|---|---|---|---|---|
| 10 | 0.70 | 3.3650 | 2.6306 | **−0.04105** | −0.12222 | 8.64 = 1.375 $T_0$ |
| 10 | 0.789 | 3.3650 | 2.5219 | **−0.04650** | −0.11233 | 9.40 = 1.496 $T_0$ |

$F<0$ trên toàn đoạn từ $\dot\theta_d$ lên $\dot\theta_{\text{gãy}}$. Không có "nạp tới ngưỡng" theo nghĩa tất định. Mọi lần gãy đều cần nhiễu đẩy $\dot\theta$ vượt ngưỡng.

## 2. Chu kỳ gãy đo trực tiếp theo $(\tau,\gamma_x)$

**Thiết lập:** 8 quỹ đạo × $60T_0$ sau giai đoạn khởi động. "Đợt thoát" là lần gãy sau khi đã khóa ổn định ít nhất $T_0$.

| $\tau$ | $\gamma_x$ | $\dot\theta_d$ | $\tau_{\text{rel}}$ ($T_0$) | $r$ (mỗi $T_0$) | $1/r$ ($T_0$) | $(1/r)/\tau_{\text{rel}}$ | Chờ khóa→gãy TB ($T_0$) | $\dot\theta$ sau khóa TB |
|---|---|---|---|---|---|---|---|---|
| 5 | 0.50 | 3.1226 | 0.488 | 0.0979 | 10.21 | 20.9 | 1.58 | 2.821 |
| 5 | 0.70 | 3.0202 | 0.521 | 0.0896 | 11.16 | 21.4 | 1.84 | 2.824 |
| 10 | 0.50 | 2.8597 | 1.163 | 0.2313 | 4.32 | 3.7 | 3.08 | 2.726 |
| **10** | **0.70** | 2.6306 | 1.375 | **0.1563** | **6.40** | 4.7 | 5.09 | **2.153** |
| 15 | 0.50 | 2.5701 | 2.160 | 0.0771 | 12.97 | 6.0 | 8.87 | 2.569 |
| 15 | 0.70 | 2.1721 | 3.024 | 0.0292 | 34.29 | 11.3 | 10.81 | −0.557 |

**Đọc kết quả:**

1. **$r=0.15$ chỉ là giá trị tại $(\tau,\gamma_x)=(10,0.70)$.** Khi đổi tham số, $r$ đổi từ 0.029 tới 0.231.
2. **$r$ thay đổi không đơn điệu theo $\tau$** (0.09 → 0.23 → 0.08 khi $\tau$ đi 5 → 10 → 15, với $\gamma_x=0.5$). Mô hình "một thời gian nạp" không giải thích được.
   - Khi $\tau=5$, $\dot\theta_d$ nằm sát ngưỡng. Hệ gãy lại ngay sau khi khóa (chưa đủ $T_0$), nên các lần này không được tính là đợt thoát mới.
   - Khi $\tau=15$ và $\gamma_x=0.7$, động cơ thường **đảo chiều hoặc gần tắt**: $\dot\theta$ sau khóa trung bình là −0.56.
3. **Phần lớn chu kỳ là thời gian chờ trong pha đồng bộ** (5.1 $T_0$ tại điểm chuẩn), cho tới khi nhiễu đẩy vượt ngưỡng. Phần nạp lại chỉ chiếm một phần nhỏ.

## 3. Tiên đề tự đồng dạng và "cận trên Carnot"

| Hạt nhân | $\hbar\omega_0$ | $E_{\text{vt}}=3\hbar\omega_0$ | $L=\hbar c/E_{\text{vt}}$ | $1.2A^{1/3}$ | Carnot $(0.42$–$0.48)\hbar\omega_0$ | 2⁺ thấp nhất | 3⁻ thấp nhất | Bao trùm? |
|---|---|---|---|---|---|---|---|---|
| ¹⁶O | 16.271 | 48.81 | 4.04 fm | 3.02 fm | 6.83–7.81 | 6.917 | 6.130 | có |
| ⁴⁰Ca | 11.988 | 35.97 | 5.49 fm | 4.10 fm | 5.04–5.75 | 3.904 | 3.737 | có |
| ²⁰⁸Pb | 6.920 | 20.76 | 9.51 fm | 7.11 fm | 2.91–3.32 | **4.085** | 2.614 | **không** |

## 4. Bảng định thức đã sửa

| Thành phần | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Không gian cấu hình, liên kết, $R\kappa=\Delta E$ | ✓ | ✓ |
| Tiên đề thang | $E_{\text{vt}}/\hbar\omega_0=3$ | ✓ Nhất quán nội tại (là một lựa chọn tiên đề) |
| Điều kiện gãy | $\dot\theta_{\text{gãy}}=3.205\,\omega_0$, $E_{\text{quay}}=51.35$ MeV | ✓ $\dot\theta_{\text{gãy}}$; $E_{\text{quay}}=5.135\,\hbar\omega_0(A)$, tức 51.35 MeV chỉ khi $\hbar\omega_0=10$ MeV |
| Cơ chế nhịp gãy | Tích thoát Van der Pol, $r\approx0.15$, $T\approx6.7T_0$ | ✗ Không nạp tất định tới ngưỡng được. $r$ do nhiễu kích hoạt nhưng gần như không phụ thuộc $k_BT$, và thay đổi 0.03–0.23 theo $(\tau,\gamma_x)$. Suy luận giải tích "6.7 $T_0$" đến từ lỗi đơn vị |
| Phát xạ | $\approx0.42\,\hbar\omega_0$, "cận trên" | Giá trị đúng tại điểm chuẩn. "Cận trên" sai với ²⁰⁸Pb (2⁺ = 4.085 MeV) |
| Langevin, tái khóa 100%, $t\approx0.4T_0$ | ✓ | ✓ (khi $\gamma_x=0$; khi có $\gamma_x$ thì TB 1.0–1.7 đv, tối đa ~6 đv) |

## Câu hỏi còn mở

- **Cơ chế thật của nhịp gãy vẫn chưa được giải thích.** Không phải Kramers (không phụ thuộc $k_BT$), cũng không phải nạp tất định ($F<0$). Một hướng cần kiểm: phân bố thời gian chờ khóa→gãy có dạng mũ (quá trình Poisson) hay có đỉnh (quá trình gần tuần hoàn)? Nếu có dạng mũ, nhịp gãy là ngẫu nhiên thuần túy với tốc độ do hình học của vùng gần ngưỡng quyết định.
- **Bảng tham số cuối vẫn có $\tau=10$, $\gamma_x=0.7$ chọn tay.** Mọi con số định lượng ($r$, $\eta$, $\Delta E$) phụ thuộc mạnh vào hai giá trị này.
