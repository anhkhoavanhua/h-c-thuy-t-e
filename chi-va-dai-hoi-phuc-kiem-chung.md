# Thang $L$ hạt nhân, tỉ số $\chi=E_{\text{kép}}/\hbar\omega_0$ và dải hồi phục $e$: giải và kiểm chứng số liệu

> Tiếp nối [thang-vu-tru-va-hoi-phuc-kiem-chung.md](thang-vu-tru-va-hoi-phuc-kiem-chung.md).
> Số liệu tính bằng [giai_so_chi_va_dai_hoi_phuc.py](giai_so_chi_va_dai_hoi_phuc.py). Chạy `python3 giai_so_chi_va_dai_hoi_phuc.py A B` (≈ 4 s) và `… C` (≈ 18 phút).
> Độ cứng tĩnh của liên kết được kiểm bằng sympy.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $L$ cỡ hạt nhân, $E_{\text{vt}}=\hbar c/L$ | **Hướng đúng**, khớp với biên đã tìm ở tài liệu trước |
| "$L\approx10^{-15}$–$10^{-13}$ m ($\sim0.2$–$200$ fm)" và "$E_{\text{vt}}\approx1$–$200$ MeV" | **Sai đổi đơn vị.** $10^{-15}$ m = **1 fm**, $10^{-13}$ m = **100 fm**. Khoảng $E_{\text{vt}}=1$–$200$ MeV ứng với $L\approx1$–$197$ fm. Với $\chi=1$, đầu trên $L=197$–$200$ fm nằm **đúng trên hoặc ngoài** biên ($k\kappa/R=1.00$–$1.03$) |
| Hướng B: chiều phụ compact $L\lesssim10^{-13}$ m "không vi phạm thực nghiệm" | **Chỉ đúng nếu chỉ hấp dẫn đi vào chiều phụ** (kiểu ADD, giới hạn thực nghiệm cỡ chục µm). Nếu trường của $\mathcal E_{\text{kép}}$ hay vật chất đi vào chiều phụ, sẽ xuất hiện các mode Kaluza–Klein khối lượng $\hbar c/L\approx2$ MeV (ở $L=100$ fm). Chúng đã bị máy gia tốc loại trừ, vì giới hạn là $1/L\gtrsim$ cỡ TeV, tức $L\lesssim10^{-19}$ m |
| $k\kappa/R=\chi(\hbar\omega_0/E_{\text{vt}})^2$ | **Đúng đại số.** Mô phỏng xác nhận: $k\kappa/R$ chỉ phụ thuộc $E_{\text{vt}}/(\sqrt\chi\hbar\omega_0)$ |
| "$k\kappa<R$ là biên tồn tại pha đồng bộ" | **Chưa chính xác.** Ở $k\kappa/R=4$, hệ vẫn giữ tĩnh được tới 21.8°–62.1° và có $\dot\theta_{\text{gãy}}>0$. Biên thật do độ cứng $\propto(Ik-m\Delta E)$ quyết định, nghĩa là phụ thuộc cả $I$ và $\Delta E$, không chỉ $\chi$ (mục 2.2) |
| $\dot\theta_{\text{gãy}}/\omega_0\approx E_{\text{vt}}/(\sqrt\chi\hbar\omega_0)$ | **Đúng tiệm cận.** Sai −0.14% tới −1.20% khi tỉ số = 20. Sai −18% tới −30% khi tỉ số = 1. Sai tới +125% khi tỉ số = 0.5 |
| $\chi\gg1$ dễ gãy; $\chi\ll1$ cần $\dot\theta$ khổng lồ | **Đúng** (ví dụ $L=1$ fm: $\dot\theta_{\text{gãy}}/\omega_0=624$ với $\chi=0.1$, còn $62.4$ với $\chi=10$) |
| $e=0$: "khóa ngay sau cú va chạm đầu tiên ($t=5.54$)" | **Sai.** $e=0$ trải qua 4 chu kỳ gãy/khóa ($\dot\theta(0)=5$) hoặc 13 chu kỳ ($\dot\theta(0)=6$). $t=5.54$ là lần khóa **cuối**, lần đầu là $t=1.41$ |
| $0<e\le0.3$: tiêu tán **cao nhất** (79.4%) | **Sai.** Cao nhất nằm ở $e=0.6$: 87.4% ($\dot\theta_0=5$) và 91.1% ($\dot\theta_0=6$). $e=0.9$ ở $\dot\theta_0=6$ còn mất 97.6% |
| Dải khóa lại là $0\le e\le0.3$ | **Đúng trong mô phỏng**, với 2 điều kiện đầu, $t\le200$: $e\le0.3$ luôn khóa, $e=0.4$ luôn không. Nhưng ngưỡng này phụ thuộc thông số, cửa sổ thời gian và ngưỡng khóa $\lvert u\rvert<10^{-3}$. Nó không phải hằng số phổ quát |
| $e\ge0.6$: "không thể tái đồng bộ", "plasma phi liên kết" | **Bị bác bỏ.** $e=0.9$, $\dot\theta_0=6$ **khóa lại** tại $t=175.65$, sau khi năng lượng tiêu tán tới 97.6%. Với mọi $e<1$, năng lượng giảm dần nên hệ có xu hướng khóa lại sau đủ lâu. $e$ lớn chỉ làm chậm quá trình này |
| Các diễn giải "đo đạc / nhiệt hóa / plasma" | **Ngoài cơ học cổ điển**, mô phỏng không kiểm chứng được. Đây là lớp diễn giải, không phải kết quả |

## 1. Thang $L$ cỡ hạt nhân (proton, $\hbar\omega_0=1$ MeV)

Các đại lượng trong bảng:
- $E_{\text{vt}}=\hbar c/L$ với $\hbar c=197.327$ MeV·fm
- $E_{\text{kép}}=\chi\hbar\omega_0$
- Biên $E_{\text{vt}}>\sqrt\chi\,\hbar\omega_0$

| $L$ (fm) | $E_{\text{vt}}$ (MeV) | $\chi$ | $\kappa/\ell_0$ | $k\kappa/R$ | Biên $\sqrt\chi\hbar\omega_0$ | $k\kappa<R$? | $\dot\theta_{\text{gãy}}/\omega_0$ |
|---|---|---|---|---|---|---|---|
| 0.2 | 986.63 | 0.1 / 1 / 10 | 0.0001 / 0.0010 / 0.0101 | $1.0\times10^{-7}$ / $10^{-6}$ / $10^{-5}$ | 0.32 / 1 / 3.16 | có | 3120 / 987 / 312 |
| 1 | 197.33 | 0.1 / 1 / 10 | 0.0005 / 0.0051 / 0.0507 | $2.6\times10^{-6}$ / $2.6\times10^{-5}$ / $2.6\times10^{-4}$ | 0.32 / 1 / 3.16 | có | 624 / 197 / 62.4 |
| 10 | 19.733 | 0.1 / 1 / 10 | 0.0051 / 0.0507 / 0.5068 | $2.6\times10^{-4}$ / $2.6\times10^{-3}$ / $0.026$ | 0.32 / 1 / 3.16 | có | 62.4 / 19.7 / 6.24 |
| 100 | 1.9733 | 0.1 / 1 / 10 | 0.0507 / 0.5068 / **5.068** | 0.026 / 0.257 / **2.57** | 0.32 / 1 / 3.16 | có / có / **không** | 6.24 / 1.97 / 0.62 |
| 197.33 | 1.0000 | 0.1 / 1 / 10 | 0.1 / **1.0** / **10** | 0.1 / **1.00** / **10** | 0.32 / 1 / 3.16 | có / **biên** / **không** | 3.16 / 1.00 / 0.32 |
| 200 | 0.9866 | 0.1 / 1 / 10 | 0.101 / **1.014** / **10.1** | 0.103 / **1.03** / **10.3** | 0.32 / 1 / 3.16 | có / **không** / **không** | 3.12 / 0.99 / 0.31 |

**Đọc bảng:**
- Khoảng an toàn là $L\lesssim10$ fm với mọi $\chi\le10$.
- Đầu trên "200 fm" của đề xuất chỉ dùng được khi $\chi<1$.
- Điều kiện $\kappa\le\ell_0$ (tức $E_{\text{kép}}\le E_{\text{vt}}$) còn chặt hơn: với $\chi=10$, nó hỏng ngay từ $L\ge\approx20$ fm.

## 2. Kiểm tra công thức $\chi$ trên động lực học đầy đủ

**Đơn vị không thứ nguyên:** $m=\omega_0=\hbar=\ell_0=1$, do đó $k=1$.

**Thông số:** $I=1$, $g=0$, $\Delta E=E_{\text{kép}}=\chi$, $E_{\text{vt}}=\varepsilon$, $\kappa=\chi/\varepsilon$, $R=\varepsilon$.

### 2.1. Số liệu

| $\chi$ | $\varepsilon$ | $\varepsilon/\sqrt\chi$ | $k\kappa/R$ | Góc giữ tĩnh max | $\dot\theta_{\text{gãy}}$ chính xác | Tiệm cận $\varepsilon/\sqrt\chi$ | Sai lệch |
|---|---|---|---|---|---|---|---|
| 0.1 | 0.1581 | 0.5 | 4.000 | 21.75° | 0.22256 | 0.5 | +124.66% |
| 0.1 | 0.3162 | 1 | 1.000 | 89.99° | 1.43613 | 1 | −30.37% |
| 0.1 | 0.6325 | 2 | 0.250 | 180° | 2.25557 | 2 | −11.33% |
| 0.1 | 1.5811 | 5 | 0.040 | 180° | 5.10860 | 5 | −2.13% |
| 0.1 | 6.3246 | 20 | 0.0025 | 180° | 20.02748 | 20 | −0.14% |
| 1 | 0.5 | 0.5 | 4.000 | 62.05° | 0.91611 | 0.5 | −45.42% |
| 1 | 1 | 1 | 1.000 | 90.00° | 1.22539 | 1 | −18.39% |
| 1 | 2 | 2 | 0.250 | 180° | 2.31901 | 2 | −13.76% |
| 1 | 5 | 5 | 0.040 | 180° | 5.17817 | 5 | −3.44% |
| 1 | 20 | 20 | 0.0025 | 180° | 20.04960 | 20 | −0.25% |
| 10 | 1.5811 | 0.5 | 4.000 | 51.80° | 0.88858 | 0.5 | −43.73% |
| 10 | 3.1623 | 1 | 1.000 | 60.68° | 1.25264 | 1 | −20.17% |
| 10 | 6.3246 | 2 | 0.250 | 180° | 2.23856 | 2 | −10.66% |
| 10 | 15.8114 | 5 | 0.040 | 180° | 5.34633 | 5 | −6.48% |
| 10 | 63.2456 | 20 | 0.0025 | 180° | 20.24257 | 20 | −1.20% |

### 2.2. Biên thật của pha đồng bộ

Lực liên kết khi đứng yên ($\dot\theta=0$), gần $\theta=0$, có độ dốc (sympy):

$$\frac{d\lambda}{d\theta}\Big|_{\theta=0,\dot\theta=0}=\frac{\kappa\,(I k-m\,\Delta E)-gI+gm\kappa^2}{I+m\kappa^2}$$

Góc giữ tĩnh vì vậy xấp xỉ:

$$\theta_{\text{giữ}}\approx\frac{R\,(I+m\kappa^2)}{\kappa\,\lvert Ik-m\Delta E\rvert}$$

**Hệ quả:**

1. Công thức $k\kappa/R$ chỉ là trường hợp riêng khi $m\Delta E\ll Ik$ và $m\kappa^2\ll I$. Khi đó $\theta_{\text{giữ}}\approx R/(k\kappa)$.
2. Khi $Ik=m\Delta E$, tức $\omega_0^2=\omega_\theta^2$ (cộng hưởng giữa hạt và động cơ), độ cứng triệt tiêu. Liên kết gần như không cần lực để giữ. Đây là lý do hàng $\chi=1$, $I=1$ giữ được tới 90° dù $k\kappa/R=1$.
3. Do đó biên tồn tại phụ thuộc **cả $I$ và $\Delta E$**, không chỉ $\chi$. Công thức $E_{\text{vt}}>\sqrt\chi\hbar\omega_0$ chỉ là tiêu chuẩn thô, đúng khi $\omega_\theta\ll\omega_0$.
4. **Chính xác hóa tài liệu trước:** với thang vũ trụ ($k\kappa/R\sim10^{17}$–$10^{77}$), liên kết không gãy "ngay khi đứng yên" theo nghĩa đen. Nó giữ được trong một góc $\theta_{\text{giữ}}\sim10^{-77}$–$10^{-17}$ rad, tức gần như không thể đồng bộ. Tài liệu trước đã được sửa câu chữ tương ứng.

## 3. Quét mịn hệ số hồi phục $e$

**Thông số:** $I=1,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.5,\ \kappa=0.2,\ R=2.5$, sao cho $R\kappa=\Delta E$. Thời gian $t\in[0,200]$. Coi là khóa khi $\lvert u\rvert<10^{-3}$.

### 3.1. $\dot\theta(0)=5$ ($E_0=12.5$)

| $e$ | Gãy | Nảy | Khóa | % $E_0$ mất | $t$ khóa đầu | Cuối | $\dot\theta$ cuối |
|---|---|---|---|---|---|---|---|
| 0.0 | 4 | 0 | 4 | 38.60% | 1.4088 | **Đồng bộ** | 3.86598 |
| 0.1 | 1 | 11 | 1 | 70.33% | 13.8735 | **Đồng bộ** | 2.70189 |
| 0.2 | 1 | 21 | 1 | 76.08% | 24.3292 | **Đồng bộ** | 2.58520 |
| 0.3 | 1 | 39 | 1 | 79.38% | 46.6186 | **Đồng bộ** | 2.33272 |
| 0.4 | 1 | 135 | 0 | 85.72% | — | Tự do | 1.63503 |
| 0.5 | 1 | 126 | 0 | 86.65% | — | Tự do | 1.85670 |
| 0.6 | 1 | 122 | 0 | **87.38%** | — | Tự do | 1.85944 |
| 0.7 | 1 | 119 | 0 | 86.75% | — | Tự do | 1.90404 |
| 0.8 | 1 | 120 | 0 | 85.36% | — | Tự do | 1.62379 |
| 0.9 | 1 | 122 | 0 | 77.67% | — | Tự do | 2.08168 |
| 1.0 | 1 | 130 | 0 | 0.00% | — | Tự do | −4.29655 |

### 3.2. $\dot\theta(0)=6$ ($E_0=18.22$)

| $e$ | Gãy | Nảy | Khóa | % $E_0$ mất | $t$ khóa đầu | Cuối | $\dot\theta$ cuối |
|---|---|---|---|---|---|---|---|
| 0.0 | 13 | 0 | 13 | 58.93% | 1.4079 | **Đồng bộ** | 3.91818 |
| 0.1 | 2 | 25 | 2 | 80.08% | 2.3081 | **Đồng bộ** | 2.80658 |
| 0.2 | 1 | 24 | 1 | 83.14% | 30.0287 | **Đồng bộ** | 2.41315 |
| 0.3 | 1 | 46 | 1 | 85.86% | 54.5853 | **Đồng bộ** | 2.09919 |
| 0.4 | 1 | 141 | 0 | 89.96% | — | Tự do | 1.62440 |
| 0.5 | 1 | 132 | 0 | 90.56% | — | Tự do | 1.53928 |
| 0.6 | 1 | 126 | 0 | 91.07% | — | Tự do | 1.50759 |
| 0.7 | 1 | 132 | 0 | 90.58% | — | Tự do | 2.00843 |
| 0.8 | 1 | 115 | 0 | 87.06% | — | Tự do | 2.09141 |
| 0.9 | 1 | 162 | 1 | **97.56%** | 175.6515 | **Đồng bộ** | −0.05990 |
| 1.0 | 1 | 132 | 0 | 0.00% | — | Tự do | −4.94411 |

### 3.3. Đọc kết quả

1. **Thời gian khóa tăng nhanh theo $e$.** Với $\dot\theta_0=5$: $t$ khóa $=1.41$ (e=0), $13.9$ (e=0.1), $24.3$ (e=0.2), $46.6$ (e=0.3), rồi vượt quá 200 ở $e\ge0.4$. Ngưỡng "0.3 / 0.4" là nơi thời gian khóa **vượt cửa sổ quan sát**. Đó không phải một chuyển pha sắc nét.
2. **$e$ cao vẫn khóa được** khi năng lượng đã cạn. Trường hợp $e=0.9$, $\dot\theta_0=6$ khóa sau 162 lần nảy, với $\dot\theta\approx-0.06$, tức pha gần như dừng hẳn.
3. **Tiêu tán lớn nhất nằm ở $e$ trung gian hoặc cao**, không phải ở $e\le0.3$. Chỉ $e=1$ không tiêu tán.
4. **Phát biểu đúng:** với mọi $0\le e<1$, hệ tiêu tán năng lượng và có xu hướng về lại trạng thái đồng bộ. $e$ càng nhỏ thì về càng nhanh. Nếu học thuyết cần "tái đồng bộ trong một thang thời gian $T$", thì dải $e$ chấp nhận được là $e<e^*(T)$, với $e^*$ tăng theo $T$. Trong mô phỏng, $e^*(200)\approx0.3$–$0.4$.

## Câu hỏi còn mở

- Giá trị cụ thể của $I$ (quán tính pha) là bao nhiêu? Biên tồn tại pha đồng bộ phụ thuộc vào $Ik-m\Delta E$, tức vào tỉ số $\omega_\theta/\omega_0$, mà hiện chưa được xác định.
- Hướng A (thang giam hãm QCD) hay hướng B (chiều phụ)? Hướng B chỉ hợp thực nghiệm nếu chỉ có hấp dẫn đi vào chiều phụ.
- Thang thời gian $T$ nào định nghĩa "tái đồng bộ"? Ngưỡng $e^*$ phụ thuộc trực tiếp vào $T$.
