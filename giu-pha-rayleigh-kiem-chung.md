# Thang giữ pha tự thân $\tau_E$ và hàm tiêu tán Rayleigh: giải và kiểm chứng số liệu

> Tiếp nối [coherence-va-g-kiem-chung.md](coherence-va-g-kiem-chung.md).
> Số liệu tính bằng [giai_so_giu_pha_rayleigh.py](giai_so_giu_pha_rayleigh.py). Chạy `python3 giai_so_giu_pha_rayleigh.py A` (vài giây) và `… B` (≈ 20 phút).
> Hai hàm `constrained_acc` và `simulate` được thêm hệ số cản tùy chọn `gth`, `gx`. Khi không cản, kết quả cũ giữ nguyên từng chữ số.
> Bộ thông số: $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1$. Một đơn vị năng lượng bằng $\hbar\omega_0=10$ MeV. $T_0=2\pi$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\theta$ là "tọa độ xyclic xấp xỉ", $p_\theta=I\dot\theta$ gần bảo toàn | **Sai.** $L$ chứa $\theta$ tường minh qua $\Delta E\cos\theta$, $gx\sin\theta$ và qua liên kết, nên $\theta$ không xyclic. Đo được: $p_\theta$ dao động 2.38 → 3.04 (gợn 12%), 1.72 → 2.56 (20%), 0.31 → 1.92 (72%) |
| $I\gg m\kappa^2$ | **Không "≫":** $m\kappa^2=0.111\,I$, tức chỉ nhỏ hơn 9 lần |
| $\tau_E=\hbar/\delta E$ với $\delta E=2\%\cdot\Delta E=0.2$ MeV cho $\approx8T_0$ | **Số học đúng** ($3.29\times10^{-21}$ s $=7.96\,T_0$), **nhưng dùng sai số liệu đầu vào.** "2%" là **gợn vận tốc** $\delta\dot\theta/\dot\theta$ của một bộ thông số khác ở tài liệu đầu ($\Delta E=0.5$, $\Omega=5$: $\Delta E/I\Omega^2=0.020$), không phải độ bất định năng lượng |
| Áp đúng công thức $\hbar/\delta E$ cho bộ thông số hiện tại | Gợn động năng thực của động cơ là $\delta K=1.80$ đơn vị $=$ **18.0 MeV** (không phải 0.2 MeV), không đổi theo tốc độ quay. Suy ra $\hbar/\delta K=3.66\times10^{-23}$ s $=$ **0.088 $T_0$**, tức ngắn hơn **90 lần** so với con số đề xuất |
| Tính "độc lập với va chạm nucleon" | Công thức $\hbar/\delta E$ với $\delta E$ là gợn cổ điển nhầm lẫn hai khái niệm: hệ bảo toàn có $E$ **không đổi**, và gợn chỉ là trao đổi động năng–thế năng, không phải độ bất định lượng tử |
| Phương trình Lagrăng có hàm Rayleigh $\Phi=\tfrac12\gamma_x\dot x^2+\tfrac12\gamma_\theta\dot\theta^2$ | **Đúng dạng** (thêm $+\gamma_\theta\dot\theta$ và $+\gamma_x\dot x$ vào hai phương trình). Đã thêm vào mô phỏng, áp cho cả pha liên kết lẫn pha tự do |
| $\Phi$ mô hình hóa **mất pha (decoherence)** | **Không.** $\Phi$ là ma sát nhớt tất định: nó làm **mất năng lượng**, không làm **ngẫu nhiên hóa pha**. Mất pha cần thêm lực ngẫu nhiên (phương trình Langevin, kèm quan hệ thăng giáng–tiêu tán) |
| $\gamma_x$ "triệt tiêu nảy hỗn loạn, tái khóa nhanh hơn nhiều" | **Bị bác bỏ.** Chỉ bật $\gamma_x=m/10$: với $e=0.2$, **0/3** lần tái khóa trong $10T_0$ (không cản thì còn 1/4). Với $e=0$: đồng bộ ở 0.66 / 1.85 / 3.28 $T_0$, so với 0.63 / 2.13 / 3.64 $T_0$ khi không cản. Nhanh hơn khoảng 10%, không phải "nhiều" |
| $\gamma_\theta=I/\tau_{\text{coherence}}$: đồng bộ **bền vững**, động cơ **vẫn quay** | **Bị bác bỏ.** Ở cả 24/24 lần chạy có $\gamma_\theta$ (với $\tau$ từ $0.2T_0$ tới $8T_0$), **động cơ tắt** (năng lượng tụt dưới rào $2\Delta E$) trước $t=10T_0$. Hệ "đồng bộ" chỉ vì mọi chuyển động dừng lại |
| Khung cuối: "tái đồng bộ 2–9 va chạm, $0.6$–$3.6\,T_0$, gọn trong $\tau_E\approx8T_0$" | **Số va chạm và thời gian đúng** (khi $g=0.1$, $e=0$, không cản). Nhưng cận "$\tau_E\approx8T_0$" sai: áp đúng công thức thì $\tau_E=0.088\,T_0$ |
| $\dot\theta_{\text{gãy}}=3.205\,\omega_0$; $E_{\text{quay}}=51.35$ MeV $\approx0.055\,m_pc^2$; $R=2359.4$ N | **Đúng** (3.20479; 51.35 MeV $=0.0547\,m_pc^2$; 2359.4 N) |
| "Hệ số hồi phục **thực nghiệm** chuẩn xác $e=0$" | **Là lựa chọn mô hình**, không phải số đo thực nghiệm. Đây là giá trị bền vững nhất trong mô phỏng (4/4) |

## 1. Độ gợn thực của động cơ khi quay đồng bộ

Động cơ quay đồng bộ ngay dưới ngưỡng gãy, $\dot\theta_{\text{gãy}}=3.20479$. Mô phỏng trong $4T_0$.

| $\dot\theta_0$ | Khoảng $\dot\theta$ (tức $p_\theta$) | $\delta\dot\theta/\dot\theta$ | $\delta K$ (động năng động cơ) | $\hbar/\delta K$ |
|---|---|---|---|---|
| $0.95\,\dot\theta_{\text{gãy}}=3.0445$ | 2.3810 → 3.0445 | 12.23% | 1.800 = 18.0 MeV | $3.66\times10^{-23}$ s = 0.088 $T_0$ |
| $0.80\,\dot\theta_{\text{gãy}}=2.5638$ | 1.7243 → 2.5638 | 19.58% | 1.800 = 18.0 MeV | 0.088 $T_0$ |
| $0.60\,\dot\theta_{\text{gãy}}=1.9229$ | 0.3122 → 1.9229 | 72.07% | 1.800 = 18.0 MeV | 0.088 $T_0$ |
| *Đề xuất* | — | "2%" | 0.2 MeV | $3.29\times10^{-21}$ s = 7.96 $T_0$ |

$\delta K$ không phụ thuộc tốc độ quay. Nó bằng phần thế năng dao động trong một vòng: rào $2\Delta E=2$ trừ đi phần của hạt đi kèm.

## 2. Chu kỳ gãy/khóa có hàm tiêu tán Rayleigh

**Cửa sổ quan sát:** $10T_0=62.83$. **Tiêu chí động cơ còn quay:** $\tfrac12I\dot\theta^2-\Delta E\cos\theta>\Delta E$, tức vượt đỉnh rào.

| Chế độ cản | $\tau$ | $f$ | $e$ | Gãy | Va chạm | $t$ đồng bộ cuối ($T_0$) | $\dot\theta(10T_0)$ | Động cơ |
|---|---|---|---|---|---|---|---|---|
| $\gamma_\theta=I/\tau,\ \gamma_x=m/\tau$ | 1.267 ($\tau_{\text{relax}}$) | 1.05 | 0 / 0.2 | 0 / 0 | 0 / 0 | không gãy | 0.000 | **tắt** |
| | | 1.50 | 0 / 0.2 | 1 / 1 | 1 / 4 | 0.54 / 1.44 | 0.000 | **tắt** |
| | | 2.00 | 0 / 0.2 | 1 / 1 | 1 / 5 | 0.47 / 1.31 | 0.000 | **tắt** |
| | 10 ($\hbar/\Gamma_{1\text{MeV}}$) | 1.05 | 0 / 0.2 | 1 / 1 | 1 / 12 | 0.39 / 3.18 | 0.106 / −0.097 | **tắt** |
| | | 1.50 | 0 / 0.2 | 1 / 1 | 1 / 11 | 0.50 / 3.22 | −0.111 / −0.092 | **tắt** |
| | | 2.00 | 0 / 0.2 | 3 / 1 | 3 / 8 | 1.07 / 1.17 | −0.105 / 0.166 | **tắt** |
| | 17.34 ($\tau_{\text{Pauli}}$) | 1.05 | 0 / 0.2 | 1 / 1 | 1 / 16 | 0.39 / 5.10 | 0.555 / −0.175 | **tắt** |
| | | 1.50 | 0 / 0.2 | 2 / 1 | 2 / 9 | 0.86 / 2.08 | −0.181 / 0.050 | **tắt** |
| | | 2.00 | 0 / 0.2 | 2 / 1 | 2 / 14 | 1.04 / 4.94 | −0.180 / −0.278 | **tắt** |
| | 50.27 ($\tau_E$ đề xuất, $8T_0$) | 1.05 | 0 / 0.2 | 1 / 1 | 1 / 23 | 0.39 / 6.57 | −1.527 / −0.504 | **tắt** |
| | | 1.50 | 0 / 0.2 | 3 / 1 | 3 / 21 | 1.36 / 7.20 | −1.715 / 0.949 | **tắt** |
| | | 2.00 | 0 / 0.2 | 4 / 1 | 4 / 18 | 1.85 / 6.46 | 1.338 / 0.121 | **tắt** |
| chỉ $\gamma_x=m/\tau$ | 10 | 1.05 | 0 / 0.2 | 2 / 1 | 2 / 24 | 0.66 / **không khóa** | 2.418 / −0.742 | còn quay / tắt |
| | | 1.50 | 0 / 0.2 | 4 / 2 | 4 / 28 | 1.85 / **không khóa** | 2.203 / 1.039 | còn quay / tắt |
| | | 2.00 | 0 / 0.2 | 8 / 2 | 8 / 30 | 3.28 / **không khóa** | 2.441 / 0.123 | còn quay / tắt |

**Đọc kết quả:**

1. **Có $\gamma_\theta$ thì động cơ luôn tắt** (24/24), kể cả với $\tau=8T_0$. Với $\gamma_\theta=I/\tau$, vận tốc giảm theo $e^{-t/\tau}$, nên động năng quay giảm theo $e^{-2t/\tau}$. Sau $10T_0=62.8$ với $\tau=50.3$, động năng chỉ còn cỡ $e^{-2.5}\approx8\%$, không đủ vượt rào $2\Delta E$. Việc "tái đồng bộ" trong các lần chạy này là trạng thái **đứng yên**, không phải đồng bộ quay. Muốn động cơ còn quay khi có $\gamma_\theta$, mô hình cần một **nguồn cấp năng lượng** (mô-men lái) bù cho tiêu tán.
2. **Chỉ $\gamma_x$** thì động cơ còn quay với $e=0$, đồng bộ ở 0.66–3.28 $T_0$, chỉ nhanh hơn không cản khoảng 10%. Với $e=0.2$, cản $\gamma_x$ **làm tệ hơn**: 0/3 lần khóa, so với 1/4 khi không cản. Cản làm giảm vận tốc tiếp cận $u$, nhưng lực $g\sin\theta$ vẫn kích hạt liên tục, nên các lần nảy không hội tụ.
3. **Kết quả chung:** hàm Rayleigh không cứu được mâu thuẫn thời gian.
   - Cản đủ mạnh ($\tau\lesssim T_0$) thì tái khóa nhanh, nhưng động cơ chết.
   - Cản yếu ($\tau\sim8T_0$) thì động cơ vẫn chết trong $10T_0$, và tái khóa với $e=0.2$ mất tới 6–7 $T_0$.

## 3. Khung mô hình cuối cùng, phần đã sửa

| Mục | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Pha đồng bộ | $\kappa=0.679$ fm, $k\kappa/R=0.1111$, giữ 180° | ✓ |
| E gãy | $\dot\theta_{\text{gãy}}=3.205\,\omega_0$, 51.35 MeV, $R=2359.4$ N, $\Delta\ddot x=-\lambda/m$ | ✓ |
| Va chạm tái khóa | $e=0$, $\Delta E=\tfrac12m_{\text{eff}}u^2$ | ✓ về công thức; $e=0$ là **lựa chọn mô hình**, không phải số đo |
| Chu trình tự ổn định | 2–9 va chạm, $0.6$–$3.6\,T_0$ | ✓ (khi $g=0.1$, $e=0$, không cản) |
| Thời gian giữ pha | $\tau_E\approx8T_0$ | ✗. Áp đúng công thức: **0.088 $T_0$**. Các ước lượng hạt nhân: 0.2–2.8 $T_0$ |
| Cơ chế mất pha | Hàm Rayleigh $\Phi$ | ✗. $\Phi$ là tiêu tán năng lượng, không phải mất pha. Với $\gamma_\theta=I/\tau$, động cơ tắt trong $10T_0$ ở 24/24 lần chạy |
| "Duy trì trạng thái đồng bộ bền vững" | Động cơ tự hạ tốc rồi quay đồng bộ mãi | ✓ khi **không có** tiêu tán liên tục. ✗ khi có $\gamma_\theta$ (động cơ dừng hẳn) |

## Câu hỏi còn mở

- **Nguồn năng lượng giữ động cơ quay.** Nếu động cơ $\mathcal E_{\text{kép}}$ có tiêu tán, cái gì cấp năng lượng để nó tiếp tục quay? Có thể mô hình hóa bằng mô-men lái $\tau_{\text{drive}}$ trong phương trình (1), cân bằng với $\gamma_\theta\dot\theta$. Khi đó sẽ có một trạng thái quay dừng để kiểm tra lại toàn bộ chu trình.
- **Mất pha thực sự** cần lực ngẫu nhiên $\xi(t)$ với $\langle\xi\xi\rangle=2\gamma k_BT\,\delta$ (Langevin). Nhiệt độ hoặc cường độ nhiễu của "chân không hạt nhân" lấy bằng bao nhiêu?
