# Kiểm mô hình $\mathcal E_{\text{kép}}$ bằng Schrödinger, Dirac, Page–Wootters và Wheeler–DeWitt

> Số liệu tính bằng [giai_so_luong_tu.py](giai_so_luong_tu.py) (cần numpy).
> Chạy `python3 giai_so_luong_tu.py`: các phần S1, S2, D, PW, PW2 mất ≈ 20 s; phần S3 mất ≈ 45 s.
> **Đơn vị:** $\hbar=m=\omega_0=\ell_0=1$; một đơn vị năng lượng bằng $\hbar\omega_0=10$ MeV; $T_0=2\pi$. Phần Dirac dùng MeV–fm.
> **Tham số mô hình:** $I=1$, $k=1$, $\Delta E=1$, $g=0.1$, $\kappa=1/3$, $R=3$.

## Tóm tắt

| Phép kiểm | Kết quả chính | Hệ quả cho học thuyết |
|---|---|---|
| **Schrödinger: phổ** của $L$ (1)–(2) | Kích thích thấp nhất **0.866** và **1.011**, trong khi mode cổ điển là $\omega_\mp=0.949$ và 1.049. Con lắc lượng tử riêng: **0.878** thay vì 1 | Ở thang ħ, thế $-\Delta E\cos\theta$ **phi điều hòa mạnh** (lệch 9–12%). Phân tích dao động nhỏ cổ điển (tài liệu đầu tiên, trường hợp B) chỉ là gần đúng |
| **Schrödinger: liên kết** $x=\kappa\sin\theta$ với lực cực đại $R=3$ | Mọi giếng "bond" có lực cực đại $R$ chỉ giữ được hạt với độ nhòe **Δy ≥ 0.52 = 1.57κ**, ở bề rộng giếng $w\approx1$. Muốn nén xuống $\Delta y\le\kappa/2$ cần năng lượng điểm không **~180 MeV**. Dao động tử tự do có $\Delta x_0=0.707>\kappa$ | **Biên độ đồng bộ $\kappa=0.68$ fm nhỏ hơn độ nhòe lượng tử (≈1.1–1.4 fm).** Liên kết "cứng" $x=\kappa\sin\theta$ không thực hiện được trong cơ học lượng tử; chỉ có thể hiểu như một trung bình |
| **Schrödinger phụ thuộc thời gian 2D** $(\theta,x)$: "E gãy" | Xác suất hạt thoát khỏi lõi liên kết tăng **dần** theo $\ell$: 0.03%, 0.2%, **1.4%** ($\ell=3$), 3.0%, 3.2%, **14.3%** ($\ell=6$; tối đa 32%). Cổ điển với cùng thế: chỉ gãy khi $\ell=6$, lúc $t=8.6\,T_0$ | Trong cơ học lượng tử, "gãy" là **xác suất tăng liên tục**, không có ngưỡng sắc. Rò rỉ bắt đầu thấy được quanh $\ell\approx3$, gần với $\dot\theta_{\text{gãy}}=3.2$ của mô hình cứng |
| **Dirac** (proton, thế vô hướng dao động tử) | Mức năng lượng lệch **−0.23% tới −1.36%** so với Schrödinger ($\hbar\omega_0=6.9$–16.3 MeV), đúng cỡ $\hbar\omega_0/m_pc^2=0.7$–1.7% | Hiệu chỉnh tương đối tính **không đáng kể**. Dùng Schrödinger là đủ (và khớp với $v/c\approx0.1$) |
| **Page–Wootters**: rô-to làm đồng hồ | Trạng thái điều kiện **tuần hoàn**: $\lvert\langle\chi(\theta)\vert\chi(\theta+2\pi)\rangle\rvert^2=1.000000$ | Số đọc góc của rô-to chỉ xác định thời gian **modulo một vòng** ($T_{\text{rot}}\approx0.38\,T_0$ tại $\dot\theta_d$). Một chu trình ba pha ($\tau_{\text{lưu}}\approx4.5\,T_0$) kéo dài **~12 vòng**, nên **góc θ không phân biệt được quá khứ, hiện tại, tương lai**. Cần thêm một bộ đếm vòng |
| **Page–Wootters**: thời gian có "nổi lên" không? (đồng hồ chậm Ω = 0.3, g = 1) | Đồng hồ **30ħ**, cửa sổ 1 vòng: F = **0.986–0.999** (F tĩnh 0.39–0.87), tức **tái lập đúng Schrödinger**. Đồng hồ **3ħ** (như mô hình): F = **0.946–0.975** | Thời gian nổi lên tốt với đồng hồ bán cổ điển. Với đồng hồ ~3ħ, sai lệch 2.5–5.4% ngay trong một vòng |
| **Wheeler–DeWitt** ($H\Psi=E\Psi$, cửa sổ năng lượng rất hẹp) | Đồng hồ 30ħ: F ≈ **0.95**. Đồng hồ 3ħ: F giảm còn **0.69–0.91**. F tĩnh **đối xứng** $t\leftrightarrow\tfrac12T_{\text{rot}}-t$ (0.79 / 0.57 / 0.57 / 0.79) | Trạng thái vượt thời gian trộn các lần quay qua cùng một số đọc θ, nên dữ liệu **không có mũi tên thời gian**. Với đồng hồ ~3ħ, thời gian nổi lên kém rõ |

## 1. Schrödinger: phổ (S1)

**Cơ sở:** $e^{i\ell\theta}$ với $\ell\in[-15,15]$, nhân với các mức dao động tử $n\le30$. Hamiltonian:

$$H=\frac{p_\theta^2}{2I}-\Delta E\cos\theta+\tfrac12p_x^2+\tfrac12x^2-g\,x\sin\theta$$

| | Giá trị |
|---|---|
| 6 mức thấp nhất | −0.0359, 0.8303, 0.9752, 1.3508, 1.8206, 1.9847 |
| $E_1-E_0$, $E_2-E_0$ | 0.8662, 1.0111 |
| Cổ điển $\omega_-$, $\omega_+$ | 0.9487, 1.0488 |
| Con lắc lượng tử riêng ($g=0$): kích thích đầu | 0.8784 (điều hòa là 1) |

## 2. Liên kết lượng tử hóa (S2)

Liên kết $x=\kappa\sin\theta$ với ngưỡng gãy $R$ được thay bằng một giếng Gauss theo $y=x-\kappa\sin\theta$:

$$V_b=-D\,e^{-y^2/2w^2},\qquad D=R\,w\sqrt e$$

Cách chọn $D$ này bảo đảm lực cực đại của giếng đúng bằng $R=3$.

| $w$ | $D$ ($\hbar\omega_0$) | Số trạng thái liên kết | $E_0$ | $\Delta y$ rms | $\Delta y/\kappa$ |
|---|---|---|---|---|---|
| 0.05 | 0.247 | 1 | −0.0002 | 12.18 | 36.5 |
| 0.10 | 0.495 | 1 | −0.0075 | 5.73 | 17.2 |
| 0.20 | 0.989 | 1 | −0.101 | 1.58 | 4.75 |
| 0.333 | 1.649 | 1 | −0.495 | 0.763 | 2.29 |
| 0.50 | 2.473 | 1 | −1.272 | 0.565 | 1.69 |
| **1.00** | 4.946 | 4 | −3.929 | **0.522** | **1.57** |
| 2.00 | 9.892 | 10 | −9.130 | 0.582 | 1.74 |

**Đọc bảng:**
- Giếng hẹp ($w$ nhỏ) thì **nông** (vì $D\propto w$), nên trạng thái liên kết lan rất rộng.
- Giếng rộng thì bản thân nó đã rộng.
- Tốt nhất đạt được là $\Delta y=0.52\approx1.6\kappa$. Vậy **không có cấu hình nào** giữ hạt sắc nét trong phạm vi κ.

## 3. "E gãy" lượng tử: Schrödinger phụ thuộc thời gian 2D (S3)

**Thiết lập:**
- Giải bằng phương pháp tách toán tử trên lưới $64\times256$ ($\theta\in[0,2\pi)$, $x\in[-12,12]$), bước $dt=0.005$, tới $10T_0$.
- Giếng bond: $w=1$, $D=4.946$.
- Trạng thái đầu: $e^{i\ell\theta}\times$ gói Gauss quanh $x=\kappa\sin\theta$, mang vận tốc đồng bộ.

| $\ell$ | $P_{\text{ngoài lõi}}$ TB (nửa sau) | $P_{\text{ngoài}}$ max | Cổ điển (cùng thế) | $t$ gãy cổ điển |
|---|---|---|---|---|
| 0 | 0.0004 | 0.0008 | không | — |
| 1 | 0.0003 | 0.0006 | không | — |
| 2 | 0.0017 | 0.0034 | không | — |
| 3 | 0.0141 | 0.0327 | không | — |
| 4 | 0.0300 | 0.0650 | không | — |
| 5 | 0.0323 | 0.0675 | không | — |
| 6 | **0.1431** | **0.3232** | **có** | 8.60 $T_0$ |

($P_{\text{ngoài}}$ là xác suất $\lvert x-\kappa\sin\theta\rvert>2w$.)

Ngưỡng cổ điển với giếng mềm ($\ell\approx6$) cao hơn ngưỡng của liên kết cứng ($\dot\theta_{\text{gãy}}=3.2$). Lý do: giếng mềm rộng hơn nhiều so với κ.

## 4. Dirac (D)

**Phương trình:** $H=\hbar c\,\sigma_x(-i\partial_x)+\sigma_z(m_pc^2+\tfrac12kx^2)$, với đạo hàm phổ và 256 điểm.

| $\hbar\omega_0$ (MeV) | Dirac $E-mc^2$ (3 mức đầu) | Schrödinger | Độ lệch | $\hbar\omega_0/m_pc^2$ |
|---|---|---|---|---|
| 6.920 (²⁰⁸Pb) | 3.4521, 10.3405, 17.1979 | 3.460, 10.380, 17.300 | −0.23%, −0.38%, −0.59% | 0.74% |
| 10.000 | 4.9835, 14.9180, 24.7883 | 5.000, 15.000, 25.000 | −0.33%, −0.55%, −0.85% | 1.07% |
| 16.271 (¹⁶O) | 8.0920, 24.1914, 40.1246 | 8.136, 24.407, 40.678 | −0.53%, −0.88%, −1.36% | 1.73% |

## 5. Page–Wootters và Wheeler–DeWitt (PW, PW2)

**Cách dựng:**
- **Trạng thái toàn cục:** $\Psi\propto\sum_je^{-(E_j-E_0)^2\tau_w^2/2}\langle j\vert\Phi_0\rangle\lvert j\rangle$, với $\Phi_0$ = (gói rô-to quanh $\theta=0$, mô-men $\ell_0$) ⊗ (hạt ở trạng thái cơ bản).
  - $\tau_w\to\infty$ cho nghiệm Wheeler–DeWitt $H\Psi=E_0\Psi$ (ở đây lấy $\tau_w=50\,T_{\text{rot}}$).
  - $\tau_w$ cỡ một vòng quay cho trạng thái Page–Wootters.
- **Trạng thái điều kiện:** $\chi(\theta)=\langle\theta\vert\Psi\rangle$, là trạng thái của hạt khi đồng hồ chỉ góc θ.
- **So sánh:** $\chi(\theta)$ được so với nghiệm Schrödinger của hạt khi bị lái bởi $\theta=\Omega t$. $F$ là độ trùng khớp; "tĩnh" là độ trùng khớp với trạng thái ban đầu (để biết hạt đã thay đổi bao nhiêu).
- Đồng hồ sạch: $\Delta E=0$.

### 5.1. Đồng hồ nhanh như mô hình ($\Omega=3$, $g=0.3$)

Trong ≤ 0.2 vòng, hạt gần như không đổi (F tĩnh ≥ 0.998), nên phép thử **không phân biệt được**: F luôn ≈ F tĩnh. Chỉ rút ra được kết luận về tính tuần hoàn: $\lvert\langle\chi(\theta)\vert\chi(\theta+2\pi)\rangle\rvert^2=1$.

### 5.2. Đồng hồ chậm ($\Omega=0.3$, $g=1$): phép thử phân biệt được

| Đồng hồ | Cửa sổ | $t=0.1\,T_{\text{rot}}$ | $0.2$ | $0.3$ | $0.4$ |
|---|---|---|---|---|---|
| 30ħ ($I=100$) | 1 vòng (0.35 $T_{\text{rot}}$) | **0.9986** (tĩnh 0.84) | **0.9987** (0.39) | **0.9969** (0.55) | **0.9857** (0.87) |
| 30ħ | 4 vòng | 0.955 | 0.954 | 0.954 | 0.952 |
| 30ħ | **WdW** | 0.955 | 0.954 | 0.954 | 0.952 |
| 3ħ ($I=10$) | 1 vòng | 0.966 (0.87) | 0.946 (0.59) | 0.962 (0.58) | 0.975 (0.84) |
| 3ħ | 4 vòng | 0.882 | 0.907 | 0.918 | 0.816 |
| 3ħ | **WdW** | **0.694** (0.90) | 0.901 (0.63) | 0.908 (0.63) | **0.723** (0.90) |

**Kết luận:**

1. **Page–Wootters hoạt động** khi đồng hồ bán cổ điển và cửa sổ không quá một vòng: trạng thái điều kiện tái lập Schrödinger với sai lệch ≤ 1.4%.
2. **Đồng hồ ~3ħ** (đúng cỡ rô-to của mô hình, $I\dot\theta_{\text{gãy}}=3.2\hbar$): sai lệch 2.5–5.4% ngay trong một vòng, và tới **31%** ở giới hạn Wheeler–DeWitt.
3. **Ở giới hạn WdW,** trạng thái điều kiện tại số đọc $t$ và $\tfrac12T_{\text{rot}}-t$ gần như trùng nhau (F tĩnh đối xứng). Phương trình vượt thời gian không phân biệt "trước" và "sau".
4. **Mũi tên thời gian của mô hình** đến từ **tiêu tán**: va chạm dẻo Carnot, cản $\gamma_x$, nhiễu. Page–Wootters và WdW là các khuôn khổ **đơn nhất, đóng**, nên không chứa tiêu tán. Muốn mô tả ba pha có mũi tên, cần một khuôn khổ **hệ mở** (ví dụ phương trình Lindblad, hoặc Page–Wootters mở rộng với môi trường).

## 6. Hệ quả cho bài báo

- **Thêm vào Mục "Semiclassical Domain of Validity":**
  - Biên độ đồng bộ κ nhỏ hơn độ nhòe điểm không (S2).
  - "Gãy" lượng tử là xác suất liên tục, không phải ngưỡng sắc (S3).
  - Hiệu chỉnh Dirac ≤ 1.4% (D).
- **Nếu bài muốn gắn ba pha quá khứ–hiện tại–tương lai với "thời gian nổi lên" (Page–Wootters/WdW),** cần nói rõ ba điều:
  - Rô-to chỉ đo được thời gian **trong một vòng**, nên cần thêm bộ đếm vòng.
  - Với ~3ħ, thời gian nổi lên kém (độ trùng khớp 0.69–0.97).
  - Mũi tên thời gian không đến từ WdW mà đến từ tiêu tán.

## Câu hỏi còn mở

- **Bộ đếm vòng.** Học thuyết có một biến nào đóng vai trò đếm số vòng của $\mathcal E_{\text{kép}}$ không? Nếu có, có thể đưa biến đó vào Page–Wootters như một đồng hồ không tuần hoàn. Nếu không, "quá khứ/tương lai" chỉ phân biệt được trong một vòng quay.
- **Hệ mở.** Bước kiểm tiếp theo là giải phương trình Lindblad cho mô hình, với các toán tử nhảy $\hbar\omega_0$ thay cho Carnot. Làm vậy vừa có mũi tên thời gian, vừa có lượng tử hóa xung (vấn đề "0.4ħω₀").
