# Thời gian giữ pha $\tau_{\text{coherence}}$ và hằng số ghép cặp $g$: giải và kiểm chứng số liệu

> Tiếp nối [bo-thong-so-hat-nhan-kiem-chung.md](bo-thong-so-hat-nhan-kiem-chung.md).
> Số liệu tính bằng [giai_so_coherence_va_g.py](giai_so_coherence_va_g.py) (chạy `python3 giai_so_coherence_va_g.py`, ≈ 4 phút).
> [giai_so_bo_thong_so_hat_nhan.py](giai_so_bo_thong_so_hat_nhan.py) được thêm tham số $g$ và cửa sổ thời gian; mặc định ($g=0$) giữ nguyên kết quả cũ.

## Tóm tắt kết luận

**Hằng số ghép cặp $g$:** phần số học đúng.

**Bảng "bộ thông số chính thức" có hai chỗ không đứng được:**
1. **$\tau_{\text{coherence}}$ phía trên bị sai 10 lần.** Tính đúng, $\hbar/\Gamma$ với $\Gamma=1$–$3$ MeV chỉ cho $0.53$–$1.59\,T_0$, không phải $5$–$15\,T_0$. Còn $\tau_{\text{Pauli}}$ thì đếm hiệu ứng chặn Pauli hai lần.
2. **Sau khi bật $g=0.1$, kết quả "$e\le0.2$ khóa 12/12" không còn đúng.** Mô phỏng lại cho **8/12** trong $10T_0$, và chỉ **4/12** trong $\tau_{\text{Pauli}}$.

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\tau_{\text{Pauli}}=(\varepsilon_F/(E-\varepsilon_F))^2\tau_{\text{relax}}=13.7\,\tau_{\text{relax}}$ | **Số học đúng:** $1.142\times10^{-21}$ s $=2.76\,T_0$ |
| Chặn Pauli kéo dài $\tau$ so với $\lambda_{\text{mfp}}/v_F$ | **Đếm hai lần.** $\lambda_{\text{mfp}}\approx5$ fm là quãng đường tự do **trong hạt nhân**, suy ra từ phần ảo của thế quang học, nên **đã bao gồm** chặn Pauli. Nhân thêm $(\varepsilon_F/(E-\varepsilon_F))^2$ là tính hiệu ứng lần thứ hai. Ngoài ra $(E-\varepsilon_F)/\varepsilon_F=0.27$, không nhỏ, nên khai triển Landau chỉ gần đúng |
| $\tau_{\text{coherence}}\sim\hbar/\Gamma_{\text{giant}}$, $\Gamma=1$–$3$ MeV, cho $(2$–$6)\times10^{-21}$ s $\approx(5$–$15)T_0$ | **Sai 10 lần.** $\hbar/\Gamma=6.58\times10^{-22}/(1$–$3)=$ **$(2.2$–$6.6)\times10^{-22}$ s $=0.53$–$1.59\,T_0$**. Với độ rộng cộng hưởng lưỡng cực khổng lồ thực tế ở hạt nhân nặng ($\approx4$–$5$ MeV) thì chỉ còn $\approx0.35\,T_0$ |
| Mode tập thể tại cộng hưởng "bảo vệ" pha | **Không có trong mô hình.** Hệ (1)–(2) bảo toàn năng lượng, không có cơ chế mất pha. Cộng hưởng chỉ tách hai mode: $\hbar\Delta\omega=1.00$ MeV, chu kỳ phách $9.99\,T_0$. Đó là trao đổi năng lượng qua lại giữa pha và hạt, không phải bảo vệ kết hợp. Việc gán nó với cộng hưởng khổng lồ của hạt nhân là một **giả thuyết đồng nhất**, không phải kết quả tính |
| $k=mω_0^2=2.41$ MeV/fm²; $g_{\max}=\sqrt{k\Delta E}=4.91$ MeV/fm $=786$ N | **Đúng:** 2.4097; 4.9088 MeV/fm; 786.5 N. Đáng chú ý: $g_{\max}$ bằng **đúng 1 đơn vị lực** $\hbar\omega_0/\ell_0$ |
| $g=0.1g_{\max}=0.49$ MeV/fm $\approx78.6$ N | **Đúng:** 0.4909 MeV/fm = 78.65 N = 0.100 đơn vị |
| $d\lambda/d\theta\vert_0\approx-0.1/1.11\approx-0.09$ | **Gần đúng, bỏ sót số hạng $+gm\kappa^2=+0.0111$.** Đầy đủ: $(0-0.1+0.0111)/1.1111=$ **−0.0800** |
| Với $g=0.1$: giữ tĩnh, gãy ở ~53 MeV | **Đúng:** giữ tĩnh 180°; $\dot\theta_{\text{gãy}}=3.20479$ (so với 3.25544 khi $g=0$); năng lượng quay tới hạn **51.35 MeV** (so với 52.99) |
| Bảng chính thức: "$0\le e\le0.2$, khóa 12/12" | **Chỉ đúng khi $g=0$.** Khi $g=0.1$: **8/12** trong $10T_0$. $e=0.2$ hỏng 3/4, $e=0.1$ hỏng 1/4, $e=0$ vẫn đạt 4/4 |
| Bảng chính thức: $\tau_{\text{coherence}}\sim(1.2$–$4.1)\times10^{-21}$ s $=3$–$10\,T_0$ | **Không có cơ sở cho đầu trên.** $4.1\times10^{-21}$ s đến từ phép tính $\hbar/\Gamma$ bị sai 10 lần. Các ước lượng đã sửa nằm trong khoảng $0.2$–$2.8\,T_0$ |

## 1. Các thang thời gian ($\hbar\omega_0=10$ MeV, $T_0=4.136\times10^{-22}$ s)

| Ước lượng | Giá trị (s) | Theo $T_0$ | Đề xuất |
|---|---|---|---|
| $\tau_{\text{relax}}=5\text{ fm}/0.2c$ | $8.34\times10^{-23}$ | 0.20 | — |
| $\hbar/\Gamma$, $\Gamma=4.5$ MeV (cộng hưởng khổng lồ, hạt nhân nặng) | $1.46\times10^{-22}$ | 0.35 | — |
| $\hbar/\Gamma$, $\Gamma=3$ MeV | $2.19\times10^{-22}$ | 0.53 | "$6\times10^{-21}$, $15T_0$" |
| $\hbar/\Gamma$, $\Gamma=1$ MeV | $6.58\times10^{-22}$ | 1.59 | "$2\times10^{-21}$, $5T_0$" |
| $\tau_{\text{Pauli}}$ (đếm hai lần chặn Pauli, tức cận trên lạc quan) | $1.14\times10^{-21}$ | 2.76 | "$1.2\times10^{-21}$, $3T_0$" ✓ |
| Thời gian mô hình cần để tái đồng bộ ($e=0$) | — | 0.6–3.6 | — |

## 2. Hằng số ghép cặp $g$

| Đại lượng | Giá trị |
|---|---|
| $k=m_pc^2(\hbar\omega_0)^2/(\hbar c)^2$ | 2.4097 MeV/fm² |
| $g_{\max}=\sqrt{k\Delta E}$ | 4.9088 MeV/fm = 786.5 N = **1.0000** đơn vị $\hbar\omega_0/\ell_0$ |
| $g=0.1g_{\max}$ | 0.4909 MeV/fm = 78.65 N = 0.100 đơn vị |
| $d\lambda/d\theta\vert_{0}$ (đủ số hạng) | −0.0800 |
| Hai mode tự do tại cộng hưởng, $\omega_\pm=\sqrt{1\pm g}$ | 1.04881, 0.94868 |
| $\hbar(\omega_+-\omega_-)$ | 1.0013 MeV |
| Chu kỳ phách | 62.753, tức $9.99\,T_0$ |

Điều kiện $g<\sqrt{k\Delta E}$ đến từ tuyến tính hóa hệ **tự do** (1)–(2) quanh $(0,0)$. Nó bảo đảm pha tự do sau khi gãy vẫn ổn định. Với $g=0.1$, điều kiện này thỏa rộng rãi: $g^2=0.01<1$.

## 3. Mô phỏng lại chu kỳ gãy/khóa với $g=0.1$

**Thông số:** $I=1$, $\Delta E=1$, $\kappa=1/3$, $R=3$, $g=0.1$ (đơn vị $m=\omega_0=\hbar=\ell_0=1$). Cửa sổ quan sát $10T_0$.

**Kết quả tĩnh:**
- Giữ tĩnh 180°.
- $\dot\theta_{\text{gãy}}=3.20479$.
- Năng lượng quay tới hạn 51.35 MeV.

**Các cửa sổ** (đơn vị $1/\omega_0$): $\tau_{\text{relax}}=1.267$, $\hbar/\Gamma_{3\text{MeV}}=3.333$, $\hbar/\Gamma_{1\text{MeV}}=10.000$, $\tau_{\text{Pauli}}=17.344$, $10T_0=62.832$.

| $f=\dot\theta_0/\dot\theta_{\text{gãy}}$ | $e$ | Gãy | Va chạm | $t$ đồng bộ cuối | $/T_0$ | $\le\tau_{\text{relax}}$ | $\le\hbar/\Gamma_3$ | $\le\hbar/\Gamma_1$ | $\le\tau_{\text{Pauli}}$ | $\le10T_0$ |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.05 | 0.0 | 2 | 2 | 3.946 | 0.63 | – | – | ✓ | ✓ | ✓ |
| 1.05 | 0.1 | 1 | 9 | 12.996 | 2.07 | – | – | – | ✓ | ✓ |
| 1.05 | 0.2 | 1 | 24 | **không khóa** | — | – | – | – | – | **✗** |
| 1.20 | 0.0 | 3 | 3 | 7.653 | 1.22 | – | – | ✓ | ✓ | ✓ |
| 1.20 | 0.1 | 2 | 17 | 31.520 | 5.02 | – | – | – | – | ✓ |
| 1.20 | 0.2 | 1 | 26 | **không khóa** | — | – | – | – | – | **✗** |
| 1.50 | 0.0 | 5 | 5 | 13.384 | 2.13 | – | – | – | ✓ | ✓ |
| 1.50 | 0.1 | 2 | 28 | **không khóa** | — | – | – | – | – | **✗** |
| 1.50 | 0.2 | 3 | 24 | 25.194 | 4.01 | – | – | – | – | ✓ |
| 2.00 | 0.0 | 9 | 9 | 22.863 | 3.64 | – | – | – | – | ✓ |
| 2.00 | 0.1 | 2 | 22 | 20.751 | 3.30 | – | – | – | – | ✓ |
| 2.00 | 0.2 | 2 | 29 | **không khóa** | — | – | – | – | – | **✗** |

### Tổng hợp theo cửa sổ

Số dòng $g=0$ lấy từ bảng ở tài liệu trước. $t/T_0$ tương ứng: 0.81, 0.59, 5.78, 1.19, 6.16, 1.17, 2.19, 4.08, 3.28, 3.03, 3.23, 8.47.

| Cửa sổ | $g=0$ | $g=0.1$ | Riêng $e=0$, $g=0.1$ |
|---|---|---|---|
| $\tau_{\text{relax}}$ (0.20 $T_0$) | 0/12 | 0/12 | 0/4 |
| $\hbar/\Gamma$, 3 MeV (0.53 $T_0$) | 0/12 | 0/12 | 0/4 |
| $\hbar/\Gamma$, 1 MeV (1.59 $T_0$) | 4/12 | 2/12 | 2/4 |
| $\tau_{\text{Pauli}}$ (2.76 $T_0$) | 5/12 | 4/12 | 3/4 |
| $10T_0$ | 12/12 | **8/12** | 4/4 |

**Đọc kết quả:**

1. **$g$ làm khó việc tái khóa.** Lực cưỡng bức $g\sin\theta$ tiếp tục kích hạt trong pha tự do, nên các lần nảy khó hội tụ về $\lvert u\rvert\to0$. Với $e>0$, 4 trong 8 lần chạy không khóa được trong $10T_0$.
2. **Chỉ $e=0$ bền vững** với $g=0.1$: 4/4 trong $10T_0$ và 3/4 trong $\tau_{\text{Pauli}}$. Dải "$0\le e\le0.2$" nên thu hẹp thành **$e\approx0$** (hoặc $e\le0.1$ nếu chấp nhận hỏng 1/4).
3. **Ngay cả với ước lượng lạc quan nhất** ($\tau_{\text{Pauli}}$, vốn đếm hai lần chặn Pauli), mô hình vẫn hỏng: chỉ 4/12 lần chạy tái đồng bộ kịp. Với các ước lượng từ độ rộng thực nghiệm (≤ 1.6 $T_0$), còn 2/12. Với $\tau_{\text{relax}}$ là 0/12.

## 4. Bảng thông số đã sửa

| Đại lượng | Giá trị | Không thứ nguyên | Trạng thái |
|---|---|---|---|
| $\hbar\omega_0$ | 10 MeV | 1 | — |
| $T_0$ | $4.136\times10^{-22}$ s | $2\pi$ | ✓ |
| $\Delta E=E_{\text{kép}}$ ($\chi=1$) | 10 MeV | 1 | ✓ |
| $I$ | $6.941\times10^{-57}$ kg·m² | 1 | ✓ |
| $E_{\text{vt}}$ ($L=6.58$ fm) | 30 MeV | 3 | ✓ |
| $\ell_0$ | 2.037 fm | 1 | ✓ |
| $\kappa$ | 0.679 fm | 1/3 | ✓ |
| $R$ | 2359.4 N (14.73 MeV/fm) | 3 | ✓ |
| $g$ | 0.4909 MeV/fm = 78.65 N | 0.100 | ✓ |
| $d\lambda/d\theta\vert_0$ | — | **−0.0800** | sửa từ −0.09 |
| $E_{\text{quay}}^{\text{tới hạn}}$ | **51.35 MeV** (khi $g=0.1$) | 5.135 | sửa từ 53.0, vì 53.0 là giá trị khi $g=0$ |
| $e$ | **$e\approx0$** (4/4); $e\le0.1$ đạt 7/8 | — | sửa từ $0\le e\le0.2$, vì $e=0.2$ chỉ đạt 1/4 khi $g=0.1$ |
| $\tau_{\text{coherence}}$ | **$(0.8$–$11)\times10^{-22}$ s $=0.2$–$2.8\,T_0$** | 1.3–17 | sửa từ $3$–$10\,T_0$; đầu trên là cận lạc quan |

## Câu hỏi còn mở

- **Mâu thuẫn thời gian vẫn còn.** Với $e=0$, $g=0.1$, việc tái đồng bộ cần 0.6–3.6 $T_0$. Mọi ước lượng dựa trên số liệu thực nghiệm cho thời gian giữ pha chỉ 0.2–1.6 $T_0$. Muốn khép kín mô hình, cần một cơ chế giữ pha có số đo **độc lập**, không suy ra từ $\lambda_{\text{mfp}}$ hay $\Gamma$ của hạt nhân.
- **Cơ chế mất pha chưa có trong mô hình.** Mô hình hiện không chứa số hạng mất pha: phương trình (1)–(2) bảo toàn năng lượng. Nếu muốn kiểm chứng "bảo vệ pha" bằng số, cần thêm vào phương trình một số hạng ngẫu nhiên hoặc ma sát với thang $\tau$, rồi đo tỉ lệ tái đồng bộ.
