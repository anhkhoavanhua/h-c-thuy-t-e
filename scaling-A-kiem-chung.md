# Thứ nguyên $\gamma_\theta$, chu trình với $\tau\sim10^3$ và luật tỉ lệ $\Delta E\propto A^{-1/3}$: giải và kiểm chứng số liệu

> Tiếp nối [langevin-day-du-kiem-chung.md](langevin-day-du-kiem-chung.md).
> Số liệu tính bằng [giai_so_scaling_A.py](giai_so_scaling_A.py) (chạy `python3 giai_so_scaling_A.py`, ≈ 3 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\gamma_\theta=\tfrac23\alpha\hbar^2/(m_pc^2)$ "đúng thứ nguyên J·s" | **Vẫn sai thứ nguyên:** $[\hbar^2/(mc^2)]=\text{J}^2\text{s}^2/\text{J}=$ **J·s²**. Khi rút gọn $\ell_0^2\omega_0^2=\hbar\omega_0/m$ bị **mất một thừa số $\omega_0$**. Dạng đúng: $\gamma_\theta=\tfrac23\alpha\hbar^2\omega_0/(m_pc^2)=5.47\times10^{-39}$ J·s |
| $\Gamma/\omega_0=\tfrac23\alpha\hbar\omega_0/(m_pc^2)$, $\tau\sim10^3$–$10^4$ | **Đúng.** Với dạng đúng: $\tau=19\,287$ ($\alpha=1/137$) |
| "$\tau\sim10^3$: nhịp gãy thưa, hạt nhân chủ yếu ở pha đồng bộ" | **Sai trong mô hình có $\gamma_x$.** Mô-men lái cực đại chỉ cỡ $\gamma_\theta=I/\tau$. Ở pha liên kết, hạt tiêu tán trung bình $\tfrac12\gamma_x\kappa^2\dot\theta^2$. Suy ra $\tau_c\approx2I/(\gamma_x\kappa^2)=$ **22.8**. Mô phỏng: với $\tau\ge20$, **động cơ tắt** sau lần gãy đầu. Nhịp gãy không "thưa" mà **dừng hẳn** |
| "$\tau=1000$: 96% tái khóa, 5.2 lần gãy" (từ Monte Carlo trước) | **Đó chỉ là giai đoạn chuyển tiếp** trước khi động cơ tắt. Mô phỏng tất định $\tau=1000$, $\gamma_x=0.789$: 1 lần gãy tại $t=0.39$, sau đó tắt ($\dot\theta\to0$) |
| Chu trình tuần hoàn khép kín (Van der Pol bù Carnot) khi có $\gamma_x$ | **Không đứng được, kể cả với $\tau=10$.** Không nhiễu: chỉ 2 lần gãy trong $50T_0$, **0 lần** trong nửa sau. Động cơ quay ổn định **dưới ngưỡng** ($\dot\theta=2.79<3.20$). Các lần gãy ở Monte Carlo trước (10.6–22.9 lần/quỹ đạo) là do **nhiễu** đẩy $\dot\theta$ vượt ngưỡng |
| $\Delta E=\eta\hbar\omega_0$, $\eta\approx0.37$–$0.41$ ⟹ $(15$–$17)A^{-1/3}$ MeV | **Số học đúng** ($0.37\times41=15.2$, $0.41\times41=16.8$). **Nhưng $\eta$ không phải hằng số:** $\gamma_x=m/\tau_{\text{relax}}$ tính theo đơn vị $1/\omega_0$ thay đổi theo $A$ (0.485 → 1.141), làm $\eta$ đổi từ **0.492 (¹⁶O) xuống 0.347 (²⁰⁸Pb)** |
| "Dự đoán độc lập không tham số tự do" | **Sai.** Có ít nhất 8 tham số đã chọn tay: $E_{\text{vt}}/\hbar\omega_0=3$, $\chi=1$, $g=0.1$, $\Omega_0=1.05\dot\theta_{\text{gãy}}$, $\tau=10$, $e=0$, $k_BT$, $\lambda_{\text{mfp}}/v_F$ |
| $\Delta E\propto A^{-1/3}$, bám mức 3⁻ thấp nhất | **Mô hình cho $\Delta E\propto A^{-0.47}$** (dốc hơn). So với 3⁻: ¹⁶O **8.01** so với 6.13 (+31%); ⁴⁰Ca **5.20** so với 3.74 (+39%); ²⁰⁸Pb **2.40** so với 2.61 (−8%). Số liệu 3⁻ (O, Pb) thì theo $A^{-0.33}$ |
| So sánh ⁵⁶Fe | **Vòng tròn:** đề xuất so "dự đoán" với chính kết quả mô phỏng dùng để lấy $\eta$, không phải với số liệu thực nghiệm |
| Mức 5⁻ của ⁴⁰Ca (4.49 MeV) gộp vào "3.74–4.49" | Khoảng này ghép **hai mức khác nhau** (3⁻ ở 3.737 và 5⁻ ở 4.49) để kéo vùng so sánh lại gần mô hình |

## 1. Thứ nguyên $\gamma_\theta$

| Biểu thức | Giá trị | Thứ nguyên |
|---|---|---|
| Đề xuất $\tfrac23\alpha\hbar^2/(m_pc^2)$ | $3.60\times10^{-61}$ | J·s² ✗ |
| Đúng $\tfrac23\alpha\hbar^2\omega_0/(m_pc^2)$ | $5.47\times10^{-39}$ | J·s ✓ |
| $I=\hbar/\omega_0$ | $6.94\times10^{-57}$ kg·m² | — |
| $\tau=I\omega_0/\gamma_\theta$ | **19 287** | không thứ nguyên |

## 2. Động cơ có duy trì được chu trình khi có cản hạt?

**Cân bằng công suất trong pha liên kết** (trung bình theo $\theta$):

$$\underbrace{\gamma_\theta\Big(1-\frac{\dot\theta^2}{\Omega_0^2}\Big)\dot\theta^2}_{\text{lái}}=\underbrace{\tfrac12\gamma_x\kappa^2\dot\theta^2}_{\text{cản hạt}}\ \Longrightarrow\ \dot\theta_{\text{dừng}}^2=\Omega_0^2\Big(1-\frac{\gamma_x\kappa^2}{2\gamma_\theta}\Big)$$

Từ đây có hai hệ quả:
- **Động cơ chỉ còn quay được** khi $\gamma_\theta>\tfrac12\gamma_x\kappa^2$, tức $\tau<\tau_c=2I/(\gamma_x\kappa^2)=22.8$.
- **Gãy tất định chỉ xảy ra** khi $\dot\theta_{\text{dừng}}>\dot\theta_{\text{gãy}}$. Với $\tau=10$: $\dot\theta_{\text{dừng}}\approx\Omega_0\sqrt{0.562}=2.52<3.20$, nên **không gãy**.

**Mô phỏng tất định** ($50T_0$, không nhiễu):

| $\tau$ | $\gamma_x$ | Gãy tổng | Gãy nửa sau ($t>25T_0$) | $t$ gãy cuối | $\dot\theta(50T_0)$ | Động cơ | Carnot/khóa |
|---|---|---|---|---|---|---|---|
| 10 | 0.789 | 2 | **0** | 3.70 | 2.791 | còn quay, **dưới ngưỡng** | 6.00 MeV |
| 15 | 0.789 | 2 | **0** | 3.78 | 1.932 | còn quay, dưới ngưỡng | 5.75 MeV |
| 20 | 0.789 | 1 | 0 | 0.39 | 0.050 | **tắt** | 6.52 MeV |
| 25 | 0.789 | 1 | 0 | 0.39 | −0.064 | **tắt** | 6.53 MeV |
| 30 | 0.789 | 1 | 0 | 0.39 | 0.001 | **tắt** | 6.53 MeV |
| 50 | 0.789 | 1 | 0 | 0.39 | 0.000 | **tắt** | 6.53 MeV |
| 100 | 0.789 | 1 | 0 | 0.39 | 0.000 | **tắt** | 6.54 MeV |
| 1000 | 0.789 | 1 | 0 | 0.39 | −0.000 | **tắt** | 6.54 MeV |
| 1000 | 0 | 3 | 1 | 300.17 | 2.997 | còn quay | 6.41 MeV |

**Đọc kết quả:**
- Với $\gamma_x=m/\tau_{\text{relax}}$, **không có** chu trình gãy tuần hoàn tất định ở bất kỳ $\tau$ nào đã thử.
- Ngưỡng tắt thực tế nằm giữa $\tau=15$ và $20$, gần với ước lượng 22.8.
- Chỉ khi bỏ cản hạt ($\gamma_x=0$), động cơ $\tau=1000$ mới giữ được nhịp gãy thưa: 3 lần trong $50T_0$, tức cỡ một lần mỗi ~100 $T_0$ ở giai đoạn sau.

## 3. $\eta(A)$ với $\gamma_x(A)$ đúng cho từng hạt nhân

**Cách tính:** $\gamma_x(A)=\hbar/(\tau_{\text{relax}}\,\hbar\omega_0(A))$, với $\tau_{\text{relax}}=8.34\times10^{-23}$ s và $\hbar\omega_0=41A^{-1/3}$ MeV. Monte Carlo 48 quỹ đạo, $\tau=10$, nhiễu 1 MeV, cửa sổ $10T_0$.

| Hạt nhân | $\hbar\omega_0$ (MeV) | $\gamma_x(A)$ | $\eta$ đo | $\Delta E$ mô hình | Công thức đề xuất $(15.2$–$16.8)A^{-1/3}$ | Mức 3⁻ thấp nhất (thực nghiệm) |
|---|---|---|---|---|---|---|
| ¹⁶O | 16.271 | 0.485 | 0.4920 | **8.01 MeV** | 6.03–6.67 | 6.130 |
| ⁴⁰Ca | 11.988 | 0.658 | 0.4335 | **5.20 MeV** | 4.44–4.91 | 3.737 |
| ⁵⁶Fe | 10.717 | 0.737 | 0.4286 | **4.59 MeV** | 3.97–4.39 | — |
| ²⁰⁸Pb | 6.920 | 1.141 | 0.3465 | **2.40 MeV** | 2.57–2.84 | 2.614 |

**Luật tỉ lệ thực của mô hình:** $\Delta E(^{16}\text{O})/\Delta E(^{208}\text{Pb})=3.34$ trên $A$ tăng 13 lần, nên $\Delta E\propto A^{-0.47}$.

Hai điểm 3⁻ thực nghiệm (O, Pb) cho $A^{-0.33}$. Mô hình **dốc hơn** số liệu, và lệch +31% / +39% / −8% ở O / Ca / Pb.

## 4. Tổng kết đã sửa

| Mục trong "tổng kết toàn diện" | Sau kiểm chứng |
|---|---|
| Hệ Euler–Lagrange + lực Rayleigh–Van der Pol | ✓ Về cấu trúc phương trình |
| E gãy: $\dot\theta_{\text{gãy}}\approx3.2\,\omega_0$, $R=2359$ N, vị trí và xung lượng liên tục | ✓ |
| Tái khóa Langevin $e=0$, xác suất 96–100% | ✓ Trong giai đoạn động cơ còn chạy. Kịp trước khi mất pha: 87–98% |
| $\Delta E\approx0.4\,\hbar\omega_0$ | ✗ $\eta$ thay đổi 0.35–0.49 theo $A$ (qua $\gamma_x$), và 0.76 khi $\gamma_x=0$ |
| "Tính tuần hoàn: chu kỳ ổn định" | ✗ Khi có $\gamma_x=m/\tau_{\text{relax}}$: không có chu kỳ tất định. Nhịp gãy chỉ do nhiễu (với $\tau\le15$), hoặc động cơ tắt (với $\tau\ge20$, bao gồm mọi $\tau$ bức xạ $10^3$–$10^4$) |
| Dự đoán $\Delta E\propto A^{-1/3}$ bám mức 3⁻ | ✗ Mô hình cho $A^{-0.47}$; lệch 8–39% |

## Câu hỏi còn mở

- **Mâu thuẫn giữa hai yêu cầu.** Muốn động cơ sống với cản hạt, cần $\tau<\tau_c\approx23$. Nhưng cản bức xạ cho $\tau\sim10^4$. Học thuyết cần một nguồn lái **mạnh hơn bức xạ khoảng 1000 lần**. Nguồn đó là gì?
- **Chu trình do nhiễu hay do động cơ?** Nếu chu kỳ gãy thật ra do nhiễu kích hoạt (như với $\tau=10$), thì tốc độ gãy phụ thuộc $k_BT$ theo kiểu Kramers, xấp xỉ $e^{-\Delta/k_BT}$. Muốn có dự đoán định lượng, cần đo tốc độ gãy theo $k_BT$.
- **Muốn có dự đoán thật không tham số tự do**, cần cố định cả 8 tham số từ số liệu độc lập **trước khi** so với phổ 3⁻. Nếu không, việc khớp chỉ là chỉnh tham số cho vừa.
