# Toán tử nhảy chung $L=\sqrt\Gamma\,e^{-i\theta}a$, bộ tiêu tán Caldeira–Leggett và ranh giới $k_BT\lesssim1$ MeV: giải và kiểm chứng số liệu

> Tiếp nối [winding-lindblad-kiem-chung.md](winding-lindblad-kiem-chung.md).
> Số liệu tính bằng [giai_so_lindblad_joint.py](giai_so_lindblad_joint.py) (cần numpy). Chạy `python3 giai_so_lindblad_joint.py`, ≈ 1.5 phút. Script **không sửa** file mô phỏng cũ.
> **Đơn vị:** $\hbar=m=\omega_0=1$, $I=1$, $\kappa=1/3$, $\gamma=\gamma_x=0.7$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $K_{\text{rotor}}=\tfrac12I\dot\theta_d^2\approx3.46$ đơn vị (34.6 MeV), $K/k_BT\approx35$ ở 1 MeV | ✓ Đúng số |
| Ở 1 MeV, "θ̇ > 0 với xác suất trên **91.6%**" | **Con số gần đúng, nhưng suy ra sai.** 91.6% là 1 − 8.4%, trong đó 8.4% là tỉ lệ **bước đếm N lùi**, không phải xác suất θ̇ < 0. Đo trực tiếp: tỉ lệ **thời gian** dΦ/dt > 0 là **91.3%** |
| "Tốc độ trôi ròng **luôn** dương" ở 1 MeV | **Sai.** Tốc độ ròng 8 quỹ đạo: +2.58, +2.40, +2.29, **−0.60**, +2.52, +2.32, +2.46, +2.47 vòng/$T_0$. Vậy **1/8 quỹ đạo quay ngược** trong suốt $60T_0$. Nên viết "trôi ròng dương ở 7/8 quỹ đạo" |
| Chọn mũi tên thời gian = $\Delta E_{\text{diss}}>0$ và $D(\rho\Vert\rho_{\text{ss}})$ giảm | ✓ Đúng hướng, khớp kiểm chứng trước |
| $L_{\text{impact}}=\sqrt{\Gamma}e^{-i\theta}a$, với $L^\dagger L=\Gamma a^\dagger a$ | ✓ Đại số đúng ($e^{-i\theta}$ unita). **Nhưng năng lượng mỗi bước nhảy sai bậc** (xem dòng dưới) |
| $e^{-i\theta}$ mô tả phản xung $\Delta\dot\theta=-\Lambda\kappa\cos\theta/I$ | **Sai độ lớn.** Va chạm cổ điển chuyển $\Delta\ell=\Lambda\kappa\lvert\cos\theta\rvert\le\kappa\sqrt{2m_{\text{eff}}\Delta E}\le$ **0.30 ħ**. $e^{-i\theta}$ chuyển **đúng 1 ħ**, gấp hơn 3 lần. Vì mô-men rô-to lượng tử hóa theo ħ, không có toán tử $e^{-i\theta}$ nào chuyển được 0.3ħ |
| Năng lượng mỗi bước nhảy | **Gấp 3.8–13.8 lần Carnot.** $\lvert\ell,1\rangle\to\lvert\ell-1,0\rangle$ lấy đi $1+(2\ell-1)/2I$. Tại $\ell=3$ (tốc độ của mô hình) là **3.5 $\hbar\omega_0$ = 35 MeV**, so với Carnot trung bình 0.40 $\hbar\omega_0$ (**gấp 8.8 lần**). Với $\ell\le0$ (rô-to quay ngược), rô-to lại **nhận** năng lượng, nên $L$ không còn là toán tử tiêu tán |
| $\Gamma_{\text{eff}}\approx\langle\Delta E\rangle/\hbar\cdot r_{\text{burst}}\cdot\langle N\rangle$ | **Sai thứ nguyên:** [1/s]·[1/s] = **1/s²**. Tốc độ va chạm cần tái lập là $r_{\text{burst}}\langle N\rangle=0.778$ mỗi $T_0$ = 0.124 mỗi đơn vị thời gian, suy ra $\Gamma=0.124/\langle a^\dagger a\rangle\approx1.74$. Nhưng với 3.5 $\hbar\omega_0$ mỗi nhảy, công suất tiêu tán sẽ **gấp 8.8 lần** công suất Carnot cổ điển |
| $H_{\text{renorm}}$ với $-g\sqrt{\hbar/2m\omega_0}(a+a^\dagger)\sin\theta$ | ✓ Đúng (vì $x=\sqrt{\hbar/2m\omega_0}(a+a^\dagger)$). Nhưng **không có số hạng liên kết**, và toán tử $\hat\lambda$ ("toán tử ứng suất") ở bảng cuối **chưa được định nghĩa**. Kiểm chứng S2 trước đã cho thấy liên kết cứng không lượng tử hóa sắc nét được |
| Caldeira–Leggett $-\tfrac{i\gamma_x}{2\hbar}[x,\{p,\rho\}]-\tfrac{\gamma_xmk_BT}{\hbar^2}[x,[x,\rho]]$ | **Hệ số nhất quán** (khi $\gamma_x$ hiểu là tốc độ tắt dần; lưu ý bài báo lại định nghĩa $\gamma_x$ là hệ số lực). **Nhưng phương trình Caldeira–Leggett chỉ hợp lệ khi $k_BT\gg\hbar\omega_0$**, trong khi ở đây $k_BT/\hbar\omega_0=0.1$. Mô phỏng ở 1 MeV: $\langle x^2\rangle\langle p^2\rangle=$ **0.010** (**vi phạm Heisenberg**, cần ≥ 0.25); **ρ có trị riêng âm −1.11** (không còn là ma trận mật độ) |
| "Phương trình Lindblad trên … thỏa định lý Spohn" | **Không áp dụng được.** Caldeira–Leggett **không có dạng Lindblad** (không bảo đảm dương hoàn toàn). Định lý Spohn đòi hỏi bộ sinh Lindblad. Mô phỏng cho thấy ρ đã âm, nên $D(\rho\Vert\rho_{\text{ss}})$ không còn xác định |
| Tái lập bội số hình học với $p_{\text{lock}}\approx0.195$ "khớp tham số cổ điển" | **Chưa chứng minh.** $p_{\text{lock}}$ được **gán** bằng giá trị cổ điển, không suy ra từ phương trình master. Muốn kiểm, cần định nghĩa "không gian con bị khóa" rồi mô phỏng quỹ đạo nhảy lượng tử |

## 1. Chiều quay ở 1 MeV

8 quỹ đạo × $60T_0$ (cùng hạt giống với [winding-lindblad-kiem-chung.md](winding-lindblad-kiem-chung.md)):

| Chỉ số | Giá trị |
|---|---|
| Tỉ lệ thời gian dΦ/dt > 0 | **91.3%** |
| Tốc độ ròng các quỹ đạo | +2.58, +2.40, +2.29, **−0.60**, +2.52, +2.32, +2.46, +2.47 vòng/$T_0$ |
| Số quỹ đạo lùi ròng | **1/8** |

## 2. Bước nhảy $e^{-i\theta}a$ so với va chạm Carnot

| $\ell$ | $\Delta E_{\text{hạt}}$ | $\Delta E_{\text{rô-to}}=((\ell-1)^2-\ell^2)/2I$ | Tổng | So với Carnot 0.40 |
|---|---|---|---|---|
| −2 | −1 | **+2.50** | +1.50 | rô-to nhận năng lượng |
| −1 | −1 | **+1.50** | +0.50 | rô-to nhận năng lượng |
| 0 | −1 | +0.50 | −0.50 | 1.2× |
| 1 | −1 | −0.50 | −1.50 | 3.8× |
| 2 | −1 | −1.50 | −2.50 | 6.2× |
| **3** | −1 | −2.50 | **−3.50** | **8.8×** |
| 4 | −1 | −3.50 | −4.50 | 11.2× |
| 5 | −1 | −4.50 | −5.50 | 13.8× |

Cổ điển, mỗi va chạm: $\Delta\ell\le\kappa\sqrt{2\cdot0.40}=$ **0.30 ħ**.

**Hệ quả:** một toán tử nhảy lượng tử hạ **nguyên một** lượng tử mô-men của rô-to không thể mô phỏng va chạm Carnot. Va chạm cổ điển chuyển ít hơn 1ħ mô-men. Có hai hướng khả dĩ:
- **(i)** Toán tử nhảy là hàm của **vận tốc tương đối** $\hat u=\kappa\cos\theta\,\hat p_\theta/I-\hat p_x$, ví dụ $L\propto\hat u$ theo kiểu tiêu tán Caldeira–Leggett trên tọa độ tương đối. Toán tử này tiêu tán liên tục và chỉ lấy năng lượng của chuyển động tương đối.
- **(ii)** Giữ $e^{-i\theta}a$ nhưng khi đó phải chấp nhận mỗi bước nhảy lấy đi $\sim(1+\ell)\hbar\omega_0$, và **tính lại** toàn bộ thống kê năng lượng. Kết quả sẽ khác mô hình cổ điển gần một bậc độ lớn.

## 3. Caldeira–Leggett so với dạng Lindblad

Xuất phát từ $\lvert0\rangle$, tiến hóa tới $15$ đơn vị thời gian:

| $k_BT$ | Bộ tiêu tán | $\langle x^2\rangle$ (đúng: cân bằng lượng tử) | $\langle x^2\rangle\langle p^2\rangle$ (≥ 0.25) | Trị riêng nhỏ nhất của ρ |
|---|---|---|---|---|
| 1 MeV | Caldeira–Leggett (đề xuất) | 0.100 (0.500) | **0.010** ✗ | **−1.11** ✗ |
| 1 MeV | Lindblad $L=Ax+iBp$ (dạng Diósi) | 0.782 (0.500), **+56%** | 0.522 ✓ | ≈ 0 ✓ |
| 5 MeV | Caldeira–Leggett | 0.500 (0.657), −24% | 0.250 (vừa chạm biên) | 0.00 |
| 5 MeV | Lindblad (Diósi) | 0.584 (0.657), −11% | 0.389 ✓ | ≈ 0 ✓ |

Kết hợp với kiểm chứng trước (B2 ở [winding-lindblad-kiem-chung.md](winding-lindblad-kiem-chung.md)), có một **đánh đổi** ở $k_BT<\hbar\omega_0$:

| Bộ tiêu tán | Cân bằng nhiệt | Công suất tiêu tán tại $\Omega=2.63$ | Dạng Lindblad (dùng được Spohn) |
|---|---|---|---|
| Nhảy nhiệt $a$, $a^\dagger$ (xấp xỉ sóng quay) | **Đúng** | Lệch **−39%** | ✓ |
| Caldeira–Leggett | **Sai** (vi phạm Heisenberg) | Đúng (ma sát vận tốc) | ✗ |
| Caldeira–Leggett + Diósi | Lệch +56% (1 MeV) | Gần đúng | ✓ |

Không bộ tiêu tán Markov nào đơn giản đạt được cả ba cột ở $k_BT=0.1\,\hbar\omega_0$. Bài Part II cần một **mô hình bể không Markov** (chuyển động Brown lượng tử với mật độ phổ đầy đủ), hoặc phải nêu rõ đang chọn đánh đổi nào.

## 4. Bảng định vị Part I / Part II đã sửa

| Cấp độ | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Cơ chế gãy (Part II) | "Toán tử ứng suất $\hat\lambda$ vượt ngưỡng" | **Chưa định nghĩa.** Liên kết cứng không lượng tử hóa sắc nét được (S2), và "gãy" lượng tử là xác suất tăng liên tục (S3) |
| Tái khóa (Part II) | $L=\sqrt\Gamma e^{-i\theta}a$ | Sai độ lớn: 1ħ so với ≤ 0.3ħ; năng lượng 3.5 so với 0.40 $\hbar\omega_0$. Đề xuất $L\propto\hat u$ (vận tốc tương đối) |
| Bể nhiệt (Part II) | Caldeira–Leggett | Không hợp lệ ở $k_BT=1$ MeV (ρ âm). Cần bể không Markov, hoặc Lindblad dạng Diósi kèm ghi rõ sai số |
| Mũi tên thời gian (Part II) | Spohn | ✓ **chỉ khi** bộ sinh có dạng Lindblad |
| Phạm vi (Part I) | $k_BT\lesssim1$ MeV, "trôi ròng luôn dương" | Giới hạn nhiệt độ hợp lý. Sửa thành "91% thời gian quay thuận; 7/8 quỹ đạo trôi ròng dương" |

## Câu hỏi còn mở

- **Toán tử vận tốc tương đối.** $\hat u$ phụ thuộc $\cos\theta\,\hat p_\theta$. Cần đối xứng hóa, ví dụ $\tfrac12\{\cos\theta,\hat p_\theta\}$, để $\hat u$ Hermite. Sau đó kiểm xem $L\propto\hat u$ có tái lập được thống kê bội số hình học và thời gian lưu mũ hay không. Đây là phép kiểm đầu tiên nên làm cho Part II.
- **Bể không Markov.** Cần một mật độ phổ $J(\omega)$ cho "chân không hạt nhân". Nếu lấy độ rộng GDR (~4–5 MeV) làm thang tần số cắt, có thể dùng mô hình Drude–Lorentz kèm phương pháp HEOM. Việc này tốn tính toán nhưng giải quyết được đánh đổi ở mục 3.
