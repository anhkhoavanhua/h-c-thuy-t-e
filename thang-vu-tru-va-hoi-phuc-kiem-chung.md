# Hợp nhất thang $\ell_0$–$L_{\text{AdS}}$, hệ số hồi phục $e$ và $E_{\text{kép}}=\Delta E$: giải và kiểm chứng số liệu

> Tiếp nối [va-cham-tai-lien-ket-kiem-chung.md](va-cham-tai-lien-ket-kiem-chung.md).
> Số liệu tính bằng [giai_so_thang_vu_tru_va_hoi_phuc.py](giai_so_thang_vu_tru_va_hoi_phuc.py) (hằng số SI CODATA, H₀ Planck 2018; RK4 $dt=10^{-4}$; chạy `python3 giai_so_thang_vu_tru_va_hoi_phuc.py` ≈ 4 phút).
> [giai_so_va_cham_tai_lien_ket.py](giai_so_va_cham_tai_lien_ket.py) được mở rộng thêm hệ số hồi phục $e$. Với $e=0$, kết quả cũ giữ nguyên từng chữ số.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\kappa=\ell_0E_{\text{kép}}/E_{\text{vt}}=\dfrac{L_{\text{AdS}}E_{\text{kép}}}{c\sqrt{\hbar m\omega_0}}$ | **Đúng đại số và thứ nguyên** (hai dạng khớp tới $10^{-9}$) |
| Chọn $\ell=\ell_0$ vì "$\kappa$ không vượt vùng định vị $\ell_0$" | **Tự mâu thuẫn khi $L$ là thang vũ trụ.** Với $L=c/H_0$: $\kappa=4.5\times10^{24}$ m, tức $\kappa/\ell_0=7\times10^{38}$, vi phạm chính điều kiện đó 38 bậc độ lớn. Điều kiện $\kappa\le\ell_0$ đòi $L\le\hbar c/E_{\text{kép}}=1.97\times10^{-13}$ m (cỡ hạt nhân) |
| Liên kết đồng bộ tồn tại được với $E_{\text{vt}}$ vũ trụ học | **Không.** Lực đàn hồi tĩnh vượt ngưỡng: $k\kappa/R=4.8\times10^{77}$ (Hubble), $2.0\times10^{17}$ (năng lượng tối), nên liên kết chỉ giữ được khi pha lệch khỏi $\theta=0$ không quá cỡ $R/(k\kappa)\sim10^{-77}$–$10^{-17}$ rad, tức gần như không thể đồng bộ (đã chính xác hóa ở [chi-va-dai-hoi-phuc-kiem-chung.md](chi-va-dai-hoi-phuc-kiem-chung.md)) |
| $e=0$ là lựa chọn đúng; $e>0$ làm hạt "nảy hỗn loạn, phá đồng bộ" | **Đúng một phần.** $e=0$ khóa nhanh nhất ($t=5.54$). $e=0.3$ **vẫn khóa** (ở $t=46.6$, sau 39 lần nảy). Chỉ $e\ge0.6$ là không khóa trong $t\le200$ |
| $e=0$ là bắt buộc để có tính bất khả nghịch (mũi tên thời gian) | **SAI.** Mọi $e<1$ đều tiêu tán. $e=0.3$–$0.9$ tiêu tán **nhiều hơn** $e=0$ (62–87% so với 38.6%). Chỉ $e=1$ là thuận nghịch (mất 0.000000) |
| Năng lượng mất khi va chạm | **Tổng quát hóa Carnot:** $\Delta E=\tfrac12(1-e^2)m_{\text{eff}}u^2$. Khớp mô phỏng tới $\le1.4\times10^{-14}$ với mọi $e$ |
| $R\kappa=\Delta E$ là "công cực đại bằng đúng chiều cao rào thế" | **SAI hệ số 2.** Rào của $-\Delta E\cos\theta$ (từ $\theta=0$ tới $\pi$) cao $2\Delta E$. Vậy $R\kappa$ bằng **nửa** rào thế |
| Đặt $E_{\text{kép}}=\Delta E$ "không mất tổng quát" | **Không chính xác.** Nó **bớt đi một tham số tự do**: $R$ bị buộc theo $\kappa$. Đây là một lựa chọn mô hình hợp lệ, nhưng có mất tổng quát |
| Sơ đồ 4 trạng thái | **Đúng với $e=0$**, kèm 2 lưu ý: $\lambda\propto\dot\theta^2$ chỉ đúng tiệm cận, và chu trình **tự dừng** (không lặp mãi) |

## 1. Thang độ dài: số liệu SI thật

**Giả thiết:** proton, $\hbar\omega_0=1$ MeV, $E_{\text{kép}}=1$ MeV. Khi đó:
- $\omega_0=1.5193\times10^{21}$ rad/s
- $\ell_0=6.4420\times10^{-15}$ m
- $k=m\omega_0^2=3.8607\times10^{15}$ N/m

Các đại lượng dùng trong bảng:
- $E_{\text{vt}}=\hbar c/L$
- $R=E_{\text{vt}}/\ell_0$
- $k\kappa/R=\hbar\omega_0E_{\text{kép}}/E_{\text{vt}}^2$. Nếu tỉ số này lớn hơn 1 thì liên kết chỉ giữ được hạt đứng yên trong một góc pha nhỏ quanh $\theta=0$, cỡ $R/(k\kappa)$ (xem tài liệu tiếp theo).
- $\dot\theta_{\text{gãy}}/\omega_0=E_{\text{vt}}/\sqrt{\hbar\omega_0E_{\text{kép}}}$

| Thang $L$ | $L$ (m) | $E_{\text{vt}}$ (eV) | $\kappa$ (m) | $\kappa/\ell_0$ | $R$ (N) | $k\kappa/R$ | $\dot\theta_{\text{gãy}}/\omega_0$ |
|---|---|---|---|---|---|---|---|
| Bán kính Hubble $c/H_0$ | $1.372\times10^{26}$ | $1.438\times10^{-33}$ | $4.481\times10^{24}$ | $6.96\times10^{38}$ | $3.576\times10^{-38}$ | $4.84\times10^{77}$ | $1.44\times10^{-39}$ |
| Năng lượng tối $\rho_\Lambda^{1/4}\approx2.24$ meV | $8.809\times10^{-5}$ | $2.240\times10^{-3}$ | $2.876\times10^{-6}$ | $4.46\times10^{8}$ | $5.571\times10^{-8}$ | $1.99\times10^{17}$ | $2.24\times10^{-9}$ |
| $L=1$ m | $1$ | $1.973\times10^{-7}$ | $3.265\times10^{-2}$ | $5.07\times10^{12}$ | $4.908\times10^{-12}$ | $2.57\times10^{25}$ | $1.97\times10^{-13}$ |
| $L=\hbar c/E_{\text{kép}}$ (biên) | $1.973\times10^{-13}$ | $1.000\times10^{6}$ | $6.442\times10^{-15}$ | $1.00$ | $2.487\times10^{1}$ | $1.00$ | $1.00$ |
| $L=1$ fm | $10^{-15}$ | $1.973\times10^{8}$ | $3.265\times10^{-17}$ | $5.07\times10^{-3}$ | $4.908\times10^{3}$ | $2.57\times10^{-5}$ | $1.97\times10^{2}$ |

**Hệ quả:**

1. Điều kiện "$\kappa\le\ell_0$" tương đương $E_{\text{vt}}\ge E_{\text{kép}}$. Hai điều kiện "$k\kappa<R$" và "$\dot\theta_{\text{gãy}}>\omega_0$" cùng tương đương $E_{\text{vt}}>\sqrt{\hbar\omega_0E_{\text{kép}}}$. Với giả thiết $E_{\text{kép}}=\hbar\omega_0$ ở đây, cả ba gặp nhau tại **cùng một biên** $E_{\text{vt}}=1$ MeV. Nghĩa là mô hình đồng bộ chỉ có ý nghĩa khi $E_{\text{vt}}\gtrsim E_{\text{kép}}$, tức $L\lesssim\hbar c/E_{\text{kép}}\approx197$ fm·MeV$/E_{\text{kép}}$.
2. Mọi thang vũ trụ học (Hubble, năng lượng tối) cho $E_{\text{vt}}$ nhỏ hơn MeV từ 9 tới 39 bậc độ lớn. Khi đó liên kết **không thể tồn tại**, và hạt luôn ở trạng thái tự do ($n=2$).
3. Vũ trụ quan sát có $\Lambda>0$ (de Sitter), nên không có bán kính AdS vật lý. Bảng trên dùng $c/H_0$ làm thang thay thế. Muốn giữ $E_{\text{vt}}=\hbar c/L_{\text{AdS}}$ thì $L_{\text{AdS}}$ phải là một thang **cỡ hạt nhân** của "vũ trụ con" trong học thuyết, không phải bán kính vũ trụ.

## 2. Hệ số hồi phục $e$

### 2.1. Phương trình va chạm tổng quát

Với vận tốc tiếp cận $u=\kappa\cos\theta\,\dot\theta^--v^-$, quy tắc Newton cho va chạm là $\dot f^+=-e\,\dot f^-$. Từ đó:

$$\Lambda=(1+e)\,m_{\text{eff}}\,u,\qquad \Delta E_{\text{mất}}=\tfrac12(1-e^2)\,m_{\text{eff}}u^2$$

- $e=0$: khóa, về lại liên kết.
- $e>0$: hạt nảy, tiếp tục tự do. Ta coi là khóa khi $\lvert u\rvert<10^{-3}$, để cắt chuỗi va chạm Zeno (vô số lần nảy ngày càng nhỏ).

### 2.2. Quét $e$ (với $E_{\text{kép}}=\Delta E$)

**Thông số:** $I=1,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.5,\ \kappa=0.2,\ R=2.5$, sao cho $R\kappa=0.5=\Delta E$. Điều kiện đầu $\dot\theta(0)=5$, $E_0=12.5$. Thời gian $t\in[0,200]$. Tốc độ tới hạn chính xác là 4.04819.

| $e$ | Gãy | Nảy | Khóa | $E$ mất | % $E_0$ | Σ Carnot $(1-e^2)$ | Lệch max | $t$ khóa cuối | Cuối | $\dot\theta$ cuối |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0 | 4 | 0 | 4 | 4.824994 | 38.60% | 4.824994 | $2.2\times10^{-15}$ | 5.5449 | **Đồng bộ** | 3.86598 |
| 0.3 | 1 | 39 | 1 | 9.922382 | 79.38% | 9.922382 | $3.8\times10^{-15}$ | 46.6186 | **Đồng bộ** | 2.33272 |
| 0.6 | 1 | 122 | 0 | 10.922454 | 87.38% | 10.922454 | $1.5\times10^{-15}$ | — | Tự do | 1.85944 |
| 0.9 | 1 | 122 | 0 | 9.708514 | 77.67% | 9.708514 | $2.3\times10^{-15}$ | — | Tự do | 2.08168 |
| 1.0 | 1 | 130 | 0 | 0.000000 | 0.00% | 0.000000 | $1.4\times10^{-14}$ | — | Tự do | −4.29655 |

**Đọc kết quả:**

- **$e=0.3$ vẫn khóa.** Vận tốc tiếp cận $\lvert u\rvert$ ở 5 lần nảy cuối giảm theo cấp số nhân: 0.092 → 0.0232 → 0.0055 → 0.0012 → 0.0002. Đây là chuỗi Zeno hội tụ, nên hạt bị bắt lại. Nhận định "$e>0$ luôn phá đồng bộ" vì vậy không đúng với $e$ nhỏ.
- **$e\ge0.6$ không khóa.** $\lvert u\rvert$ dao động quanh 0.12–0.41 mà không hội tụ. Trong khoảng này, nhận định "nảy hỗn loạn" được xác nhận trong cửa sổ $t\le200$.
- **Tiêu tán không cần $e=0$.** Với $e=0$, hệ khóa sớm (sau 4 chu kỳ) nên chỉ mất 38.6%. Với $e$ trung gian, hạt va chạm hàng trăm lần và mất nhiều hơn (tới 87.4%). Chỉ $e=1$ bảo toàn năng lượng tuyệt đối. Ở đó pha còn **đảo chiều** ($\dot\theta=-4.30$) do trao đổi năng lượng với hạt.
- **Carnot tổng quát khớp mô phỏng** tới $10^{-14}$ ở cả 5 giá trị $e$.

**Kết luận về $e$:** chọn $e=0$ là hợp lệ. Nhưng lý do đúng là **$e=0$ khóa nhanh và chắc chắn nhất**, không phải "chỉ $e=0$ mới bất khả nghịch". Còn việc diễn giải va chạm dẻo là "phép đo, sụp đổ trạng thái lượng tử" thì nằm ngoài cơ học cổ điển, mô phỏng không kiểm chứng được.

## 3. $E_{\text{kép}}=\Delta E$

- **Rào thế:** $\Pi_{\mathcal E}=-\Delta E\cos\theta$ có đáy $-\Delta E$ (tại $\theta=0$) và đỉnh $+\Delta E$ (tại $\theta=\pi$). Chiều cao rào vì thế là $2\Delta E$. Hệ thức $R\kappa=E_{\text{kép}}=\Delta E$ nghĩa là công cực đại của lực liên kết bằng **một nửa** rào thế.
- **Tính tổng quát:** đồng nhất $E_{\text{kép}}=\Delta E$ ràng buộc $R=\Delta E/\kappa$, tức bớt một tham số tự do. Đây là lựa chọn tối giản hợp lệ, nhưng có làm mất tính tổng quát.
- **Hệ quả số:** mô phỏng ở tài liệu trước dùng $R=3,\ \kappa=0.2$, tức $R\kappa=0.6\neq\Delta E$, nên chưa tuân theo quy ước này. Chạy lại với $R=2.5$ (để $R\kappa=\Delta E$), dòng $e=0$ ở bảng 2.2 cho:
  - Tốc độ tới hạn 4.04819, thay cho 4.33634.
  - 4 chu kỳ gãy/khóa, thay cho 3.
  - Mất 38.60% năng lượng, thay cho 28.51%.

## 4. Sơ đồ chu trình đã hiệu chỉnh

$$\begin{matrix}
\textbf{1. Đồng bộ} & \xrightarrow{\ \lvert\lambda\rvert=R=E_{\text{vt}}/\ell\ } & \textbf{2. Gãy} \\
\left[n=1,\ E=\text{const},\ \lambda\to m\kappa\dot\theta^2\ \text{(tiệm cận)}\right] & & \left[\Delta x=0,\ \Delta p=0,\ \Delta\ddot x=-\lambda/m\right] \\
\uparrow\ (e=0\ \text{hoặc}\ e\ \text{nhỏ, Zeno}) & & \downarrow \\
\textbf{4. Va chạm tại } f=0 & \xleftarrow{\ f(\theta,x)=0\ } & \textbf{3. Tự do} \\
\left[\Lambda=(1+e)m_{\text{eff}}u,\ \Delta E=\tfrac12(1-e^2)m_{\text{eff}}u^2\right] & & \left[n=2,\ \omega_0\right]
\end{matrix}$$

**Điều kiện tồn tại chu trình:** $E_{\text{vt}}\gtrsim E_{\text{kép}}$ (mục 1). Nếu không, hệ gần như luôn ở trạng thái 3.

**Điều kiện dừng:** chu trình kết thúc khi $\dot\theta$ tụt dưới $\dot\theta_{\text{gãy}}$. Khi đó hệ nằm lại ở trạng thái 1. Với $e\ge0.6$, hệ có thể dao động giữa trạng thái 3 và 4 rất lâu mà không về được trạng thái 1.

## Câu hỏi còn mở

- $L$ trong $E_{\text{vt}}=\hbar c/L$ là thang nào của học thuyết? Để mô hình đồng bộ tồn tại, $L$ phải $\lesssim1.97\times10^{-13}$ m (với $E_{\text{kép}}=1$ MeV). Thang vũ trụ học bị loại.
- $E_{\text{kép}}$ và $\hbar\omega_0$ có quan hệ gì? Bảng mục 1 đang giả thiết cả hai bằng 1 MeV. Đổi tỉ số này sẽ dịch biên $E_{\text{vt}}$ theo $\sqrt{\hbar\omega_0E_{\text{kép}}}$.
