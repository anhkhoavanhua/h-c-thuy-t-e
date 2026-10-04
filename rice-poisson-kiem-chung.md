# Nhịp gãy "Poisson băng ngưỡng (Rice)", cửa sổ tham số và cận Carnot đơn hạt: giải và kiểm chứng số liệu

> Tiếp nối [recharge-kiem-chung.md](recharge-kiem-chung.md).
> Số liệu tính bằng [giai_so_rice_poisson.py](giai_so_rice_poisson.py) (chạy `python3 giai_so_rice_poisson.py`, ≈ 6 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \tau=10,\ \gamma_x=0.70,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}$ ($\dot\theta_{\text{gãy}}=3.2048$). $T_0=2\pi$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| OU quanh $\dot\theta_d$: $\sigma_u^2=\gamma_\theta k_BT/(I_{\text{eff}}\lvert F'\rvert)$ | **Đúng dạng** khi chỉ có $\xi_\theta$. Còn thiếu phần nhiễu hạt chiếu qua $\kappa\cos\theta$: $+\tfrac12\gamma_x\kappa^2k_BT/(I_{\text{eff}}\lvert F'\rvert)$, làm $\sigma_u^2$ tăng 39% |
| Tốc độ Rice $r\propto\exp(-\Delta u^2/2\sigma_u^2)$ | **Bị bác bỏ, cùng lý do với Kramers.** Vì $\sigma_u^2\propto k_BT$, Rice dự đoán $r(1\text{ MeV})/r(5\text{ MeV})=$ **0.18** (đề xuất) hoặc **0.29** (có nhiễu hạt). **Đo được 1.12** |
| "Tăng $k_BT$ đẩy điểm dừng lên gần ngưỡng, nên $\Delta u/\sigma_u$ đổi chậm" | **Sai.** Đo ở pha đồng bộ không cho gãy: $\dot\theta$ tại mặt cắt $\theta\equiv0$ trung bình **2.857** (1 MeV), rồi **0.53 / −0.37 / 0.15** (3 / 5 / 8 MeV), với $\sigma\approx2.8$–$3.0$. Nhiễu mạnh làm rô-to **mất trật tự quay** (dừng, đảo chiều), chứ không dời lên gần ngưỡng |
| Thời gian chờ phân bố mũ (Poisson) | **Chỉ đúng một phần.** Phân bố **hai thành phần**: **68.7%** (1 MeV) và **82.0%** (5 MeV) là **gãy lại tức thì** (dưới $0.2T_0$) sau khi khóa. Hệ số biến thiên tổng CV = **2.5–3.9** (Poisson thì CV = 1). Riêng các lần chờ ≥ $T_0$: phần dư $w-T_0$ có **dạng mũ tốt** ở 1 MeV (CV 1.006; hàm sống 0.579/0.373/0.135/0.048 so với mũ 0.607/0.368/0.135/0.050) và gần mũ ở 5 MeV (CV 1.256) |
| $r\approx0.15$/$T_0$ do hình học không gian pha | **Không đúng như mô tả.** Tần suất vượt ngưỡng ứng suất ($\max\lvert\lambda\rvert/R\ge1$ mỗi vòng) ở pha đồng bộ **tăng gấp đôi** theo $k_BT$ (0.65 → 1.35 mỗi $T_0$) và **lớn hơn 4–9 lần** $r$ đo được. Tốc độ "đợt thoát" giữ ~0.15 vì phần lớn lần gãy rơi vào chuỗi gãy lại tức thì, không được tính là đợt mới |
| $(\tau\approx10,\gamma_x\approx0.7)$ là "vùng **duy nhất**" có chu trình | **Sai.** Bảng tài liệu trước: $(\tau,\gamma_x)=(10,0.5)$ có $r=0.231$, $(15,0.5)$ có $r=0.077$, $(5,0.5)$ có $r=0.098$, đều có chu trình. Một **khoảng** $5<\tau<15$ không cố định $\tau=10$, và $r$ đổi 0.08–0.23 trong khoảng đó |
| $\gamma_x\approx0.7$ "cố định trực tiếp từ $\tau_{\text{relax}}$" | **Không đúng.** $\tau_{\text{relax}}=8.34\times10^{-23}$ s cho $\gamma_x=0.789$ tại $\hbar\omega_0=10$ MeV, và 0.485–1.141 theo $A$. Giá trị 0.70 đến từ giả thiết $\lambda_{\text{mfp}}\propto A^{1/3}$ (không có căn cứ, xem tài liệu Kramers) |
| ²⁰⁸Pb là hạt nhân magic kép nên 2⁺ bị đẩy cao | **Sự kiện đúng** ($Z=82$, $N=126$). **Nhưng lập luận không phân biệt được:** ¹⁶O ($Z=N=8$) và ⁴⁰Ca ($Z=N=20$) **cũng magic kép**, mà "cận trên" vẫn đúng với chúng. Lập luận không giải thích được vì sao chỉ Pb vi phạm |
| Đổi tên thành "cận Carnot đơn hạt" | Chỉ là **đổi tên**, không phải kiểm chứng. Khi một cận trên bị vi phạm, việc thu hẹp phạm vi sau đó cần một tiêu chí độc lập |
| $t_{\text{khóa}}\approx1$–$2$ đv | ✓ Khớp đo (TB 1.03–1.69 đv khi có $\gamma_x$) |

## 1. Dự đoán Rice/OU theo công thức đề xuất

**Thông số:** $\dot\theta_d=2.6306$, $\Delta u=\dot\theta_{\text{gãy}}-\dot\theta_d=0.5742$, $\lvert F'\rvert=0.1222$, $I_{\text{eff}}=1.056$.

| $k_BT$ | Nhiễu | $\sigma_u^2$ | $\sigma_u$ | $\Delta u^2/2\sigma_u^2$ | $e^{-(\cdot)}$ |
|---|---|---|---|---|---|
| 1 MeV | chỉ $\xi_\theta$ (đề xuất) | 0.0775 | 0.278 | 2.127 | 0.119 |
| 1 MeV | $\xi_\theta+\xi_x$ | 0.1077 | 0.328 | 1.531 | 0.216 |
| 5 MeV | chỉ $\xi_\theta$ | 0.3876 | 0.623 | 0.425 | 0.654 |
| 5 MeV | $\xi_\theta+\xi_x$ | 0.5383 | 0.734 | 0.306 | 0.736 |

**Tỉ số $r(1)/r(5)$:** Rice dự đoán **0.182** hoặc **0.294**. Đo được **1.123**.

## 2. Pha đồng bộ khi không cho gãy

**Thiết lập:** 8 quỹ đạo × $60T_0$ sau khởi động. $M$ là $\max\lvert\lambda\rvert/R$ trong mỗi vòng quay.

| $k_BT$ | $\dot\theta(\theta\equiv0)$ TB | $\sigma$ | $(\dot\theta_{\text{gãy}}-\text{TB})/\sigma$ | $M$ TB | $\sigma_M$ | $P(M\ge1)$ mỗi vòng | Vòng mỗi $T_0$ | Vượt ngưỡng mỗi $T_0$ |
|---|---|---|---|---|---|---|---|---|
| 1 MeV | 2.857 | 0.294 | 1.18 | 0.877 | 0.189 | 0.263 | 2.47 | **0.650** |
| 3 MeV | 0.531 | 2.800 | 0.95 | 0.897 | 0.349 | 0.383 | 2.28 | **0.875** |
| 5 MeV | −0.366 | 2.912 | 1.23 | 0.978 | 0.402 | 0.451 | 2.31 | **1.042** |
| 8 MeV | 0.152 | 3.027 | 1.01 | 1.086 | 0.511 | 0.563 | 2.40 | **1.352** |

**Đọc bảng:**
- Từ 3 MeV trở lên, $\dot\theta$ tại mặt cắt trung bình về gần 0 với $\sigma\approx3$: rô-to quay hai chiều, không còn một điểm dừng $\dot\theta_d$ rõ ràng. Mô hình "OU quanh $\dot\theta_d$" chỉ còn áp được ở khoảng 1 MeV.
- Tần suất vượt ngưỡng ứng suất tăng gấp đôi theo $k_BT$. Nó phụ thuộc $k_BT$, trong khi tốc độ đợt thoát $r$ (0.140–0.158) thì không.

## 3. Phân bố thời gian chờ khóa → gãy (cho phép gãy)

**Thiết lập:** 16 quỹ đạo × $60T_0$.

| $k_BT$ | Số lần chờ | TB ($T_0$) | CV | Gãy lại tức thì (< $0.2T_0$) | Chờ ≥ $T_0$: số lần / TB / CV | Phần dư $w-T_0$: TB / CV |
|---|---|---|---|---|---|---|
| 1 MeV | 651 | 1.171 | **2.506** | **68.7%** | 126 / 5.535 / 0.824 | 4.535 / **1.006** |
| 5 MeV | 1 498 | 0.495 | **3.855** | **82.0%** | 146 / 4.510 / 0.977 | 3.510 / 1.256 |

**Hàm sống của phần dư** $P(w-T_0>k\cdot\text{TB})$, đo so với $e^{-k}$:

| $k$ | 1 MeV | 5 MeV | $e^{-k}$ |
|---|---|---|---|
| 0.5 | 0.579 | 0.548 | 0.607 |
| 1 | 0.373 | 0.342 | 0.368 |
| 2 | 0.135 | 0.137 | 0.135 |
| 3 | 0.048 | 0.068 | 0.050 |

**Kết luận về cấu trúc nhịp gãy:**

1. **Chuỗi gãy lại tức thì** (69–82%): vừa khóa lại thì $\lvert\lambda\rvert$ đã ≥ $R$, nên gãy ngay. Tỉ lệ này tăng theo $k_BT$.
2. **Pha đồng bộ dài**: thời gian lưu có dạng **mũ**, trung bình khoảng 3.5–4.5 $T_0$ sau $T_0$ đầu. Thành phần này đúng là **Poisson**, nhưng tốc độ của nó **không** theo công thức Rice.
3. Tốc độ "đợt thoát" ~0.15 mỗi $T_0$ thực chất là **nghịch đảo thời gian lưu trung bình của pha đồng bộ dài** (≈ 1 + 4.5 $T_0$). Nó gần như không phụ thuộc $k_BT$. Cơ chế đặt ra thời gian lưu này **vẫn chưa được xác định**.

## 4. Bảng chuẩn hóa đã sửa

| Khía cạnh | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Không gian trạng thái, pha đồng bộ, tiên đề tỉ lệ, E gãy, tái liên kết | ✓ | ✓ (như các tài liệu trước) |
| Động học chu kỳ gãy | Poisson băng ngưỡng Rice quanh $\dot\theta_d$, $r\approx0.15$ | **Hai thành phần:** chuỗi gãy lại tức thì (69–82%) cộng pha đồng bộ dài có thời gian lưu dạng mũ (CV ≈ 1). **Không theo Rice** (dự đoán 0.18–0.29 so với đo 1.12). $r=0.15$ chỉ tại $(\tau,\gamma_x)=(10,0.7)$ |
| Cửa sổ tham số | "$(\tau\approx10,\gamma_x\approx0.7)$ duy nhất" | Chu trình tồn tại trong một vùng rộng hơn ($\tau=5$–15, $\gamma_x=0.5$–0.7) với $r=0.03$–0.23 |
| Phát xạ | "Cận Carnot đơn hạt" $0.42\,\hbar\omega_0$ | Giá trị đúng tại điểm chuẩn. Vi phạm ở ²⁰⁸Pb vẫn chưa được giải thích, vì O và Ca cũng là hạt nhân magic kép |

## Câu hỏi còn mở

- **Thời gian lưu dạng mũ ~4.5 $T_0$ do đâu?** Nó không phụ thuộc $k_BT$ và lớn hơn nhiều $\tau_{\text{rel}}=1.4\,T_0$. Hướng kiểm tiếp: đo thời gian lưu theo $g$ (hệ số ghép pha–hạt), vì $g\sin\theta$ là nguồn cưỡng bức tất định duy nhất còn lại trong pha đồng bộ.
- **Có muốn tính các lần gãy lại tức thì như một phần của "một biến cố gãy" không?** Nếu có, mỗi "E gãy" thực chất là một **chuỗi** gãy và khóa liên tiếp, trung bình 651/126 ≈ **5.2** lần (1 MeV) và 1498/146 ≈ **10.3** lần (5 MeV). Khi đó năng lượng phát xạ mỗi biến cố bằng tổng Carnot của cả chuỗi, không phải của từng xung.
