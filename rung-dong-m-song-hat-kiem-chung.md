# Rung động của khối lượng $m$ sinh ra các mức năng lượng → lúc sóng, lúc hạt, lúc đồng thời: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I: $\hbar\omega_0=10$ MeV, $\kappa=\ell_0/3$, $k_BT=1$ MeV, $\gamma_x=0.70$, $f=1.05$.
> Số liệu tính bằng [giai_so_rung_dong_song_hat.py](giai_so_rung_dong_song_hat.py). Lệnh `python3 giai_so_rung_dong_song_hat.py` chạy ≈ 37 s trên 8 lõi.
> Phần B dùng lại `run()` của module Langevin chung qua một bản sao trong bộ nhớ, có ghi thêm $(x,v)$. File gốc **không bị sửa**.
> Bài này giải ý tưởng gốc của học thuyết theo hướng xây dựng. Bài trước, [song-hat-ba-pha-kiem-chung.md](song-hat-ba-pha-kiem-chung.md), đã kiểm một ánh xạ cụ thể hơn và bác bỏ nó.

## Thước đo sóng–hạt cho khối lượng $m$ đang rung

Khối lượng $m$ dao động với chu kỳ $T_0$. Mình so trạng thái của nó tại $t$ với trạng thái sau nửa chu kỳ, khi hạt đã sang phía bên kia:

| Đại lượng | Định nghĩa | Ý nghĩa |
|---|---|---|
| **V (sóng)** | $\lvert\langle\psi(t)\vert\psi(t+T_0/2)\rangle\rvert$ | V = 1: trạng thái không đổi, là sóng dừng |
| **D (hạt)** | Khoảng cách vết giữa hai trạng thái | D = 1: phân biệt được hạt ở bên này hay bên kia, tức là một hạt định xứ đang chạy |

Hai đại lượng tuân theo hệ thức lưỡng tính Englert $D^2+V^2\le1$, với dấu bằng khi trạng thái thuần.

## A. Năng lượng rung quyết định sóng hay hạt

Giải số trong cơ sở Fock ($n\le60$). So hai loại trạng thái của $m$:
- **trạng thái kết hợp:** rung trải trên **nhiều mức** năng lượng;
- **trạng thái riêng $\lvert n\rangle$:** đúng **một mức**.

| $E_{\text{rung}}$ | Biên độ / $\ell_0$ | V (kết hợp) | D (kết hợp) | $D^2+V^2$ | Độ trải năng lượng ΔE | V, D của $\lvert n\rangle$ | Nhìn thấy |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 1.000 | 0.000 | 1 | 0 | 1, 0 | **sóng** |
| 0.53 MeV | 0.33 | 0.899 | 0.437 | 1 | 2.3 MeV | 1, 0 | ranh giới sóng |
| **1.73 MeV** | 0.59 | **0.707** | **0.707** | 1 | 4.2 MeV | 1, 0 | **đồng thời, cân bằng** |
| 3.8 MeV | 0.87 | 0.468 | 0.884 | 1 | 6.2 MeV | 1, 0 | đồng thời |
| 4.15 MeV | 0.91 | 0.436 | 0.900 | 1 | 6.4 MeV | 1, 0 | ranh giới hạt |
| 10 MeV | 1.41 | 0.135 | 0.991 | 1 | 10 MeV | 1, 0 | **hạt** |
| 20.5 MeV | 2.03 | 0.017 | 1.000 | 1 | 14 MeV | 1, 0 | hạt |
| 30 MeV ($E_{\text{vt}}$) | 2.45 | 0.003 | 1.000 | 1 | 17 MeV | 1, 0 | hạt |

**Công thức giải tích**, khớp số tới 4 chữ số: $V=e^{-2E_{\text{rung}}/\hbar\omega_0}$ và $D=\sqrt{1-V^2}$.

**Kết quả:**
1. **Ý tưởng của học thuyết đứng được ở dạng chính xác.** Năng lượng rung của $m$ quyết định thứ mình nhìn thấy:
   - $E_{\text{rung}}<0.53$ MeV: thấy **sóng**;
   - 0.53–4.15 MeV: thấy **đồng thời cả hai**;
   - $>4.15$ MeV: thấy **hạt**.
   $D^2+V^2=1$ chính xác ở mọi mức. Có bao nhiêu "sóng" thì bớt đi bấy nhiêu "hạt".
2. **Điều kiện cần là "năng lượng khác nhau" phải có mặt cùng lúc.** Trạng thái $\lvert n\rangle$ có đúng một năng lượng nên **luôn là sóng dừng thuần** (V = 1, D = 0) ở mọi mức. Tính hạt chỉ xuất hiện khi $m$ rung đồng thời trên nhiều mức, với độ trải ΔE vài MeV. Nói cách khác: *năng lượng xác định → sóng; nhiều năng lượng chồng chập → hạt.*
3. Một va chạm Carnot (4.0–4.6 MeV) nằm **đúng ranh giới hạt** D = 0.9 (4.15 MeV).

## B. Dọc theo động cơ $\mathcal E_{\text{dual}}$ (Bảng I, 1 MeV)

Mô phỏng 8 quỹ đạo × 60 $T_0$, lấy mẫu mỗi 0.01. Khối lượng $m$ bị lái theo quỹ đạo cổ điển, nên trạng thái lượng tử của nó là trạng thái kết hợp với $\alpha=(x+iv)/\sqrt2$ (định lý dao động tử cưỡng bức). Năng lượng rung là $E_{\text{rung}}=\tfrac12(x^2+v^2)\,\hbar\omega_0$.

| Pha | Tỉ lệ thời gian | $E_{\text{rung}}$ TB | Max | V TB | D TB | Thời gian thấy **sóng** | **Đồng thời** | **Hạt** |
|---|---|---|---|---|---|---|---|---|
| Quá khứ (khóa) | 60% | 1.80 MeV | 5.8 | 0.71 | 0.67 | 1.0% | **94.8%** | 4.3% |
| Ngưỡng ($\pm0.05\,T_0$ quanh lần gãy) | 6% | 1.00 MeV | 5.2 | 0.83 | 0.55 | 7.3% | **92.4%** | 0.3% |
| Hiện tại (tự do, va chạm) | 12% | **0.97 MeV** | 6.9 | **0.83** | 0.52 | **31.7%** | 67.7% | 0.6% |
| Tương lai ($T_0$ đầu sau tái khóa) | 22% | **2.18 MeV** | 5.9 | 0.67 | **0.72** | 0.0% | 91.6% | **8.4%** |

**Kết quả:**
1. **Hệ thực sự lúc thấy sóng, lúc thấy hạt, nhưng chủ yếu là đồng thời.** Ở mọi pha, 68–95% thời gian $m$ nằm trong vùng "đồng thời". Sóng thuần và hạt thuần chỉ là những khoảnh khắc.
2. **Ánh xạ pha ngược với đoạn diễn giải trước:**
   - **Hiện tại là pha "sóng" nhất:** 32% thời gian là sóng thuần, V = 0.83. Khi liên kết gãy, rô-to thôi lái hạt. Ma sát $\gamma_x=0.7$ hãm hạt về gần điểm không trong ~0.2 $T_0$, nên năng lượng rung giảm còn ~1 MeV.
   - **Quá khứ và Tương lai "hạt" hơn:** D = 0.67–0.72, có 4–8% thời gian là hạt thuần. Lý do là rô-to lái $m$ theo $x=\kappa\sin\theta$ với vận tốc $\kappa\dot\theta$, đẩy năng lượng rung lên 1.8–2.2 MeV.
3. **Ở ngưỡng gãy:** V = 0.83, D = 0.55, và 92% thời gian là "đồng thời". Đây là chỗ đoạn diễn giải nói "đồng thời cả hai". Số liệu xác nhận điều này, nhưng ngưỡng không khác biệt nhiều so với các pha khác.

## Kết luận

| Ý của học thuyết | Kết quả số |
|---|---|
| Rung động của $m$ sinh ra các năng lượng khác nhau | ✓ Các mức $(n+\tfrac12)\hbar\omega_0$; rung kiểu cổ điển là chồng chập nhiều mức với ΔE = √E·ħω₀ |
| Từ đó lúc thấy sóng, lúc thấy hạt, lúc đồng thời | ✓ $V=e^{-2E/\hbar\omega_0}$, $D=\sqrt{1-V^2}$, $D^2+V^2=1$. Ranh giới: 0.53 MeV (sóng) và 4.15 MeV (hạt) |
| Sóng ↔ Quá khứ, hạt ↔ Hiện tại | ✗ Ngược lại: Hiện tại "sóng" nhất (V = 0.83), Tương lai "hạt" nhất (D = 0.72) |
| Đồng thời ở ngưỡng | ✓ một phần: 92% "đồng thời", nhưng các pha khác cũng vậy (68–95%) |

**Câu có thể đưa vào học thuyết:**

> The vibrating mass $m$ displays wave- or particle-like behaviour according to its vibrational energy above zero point: the half-period overlap $V=e^{-2E_{\text{vib}}/\hbar\omega_0}$ and distinguishability $D=\sqrt{1-V^2}$ saturate the Englert duality relation $D^2+V^2=1$. Along the engine cycle at $k_BT=1$ MeV the mass is in the intermediate regime 68–95% of the time; the free Present phase is the most wave-like ($\langle V\rangle=0.83$), because the particle, released from the rotor drive, is damped towards its zero-point state.

## Câu hỏi còn mở
- Muốn hạt thuần (D ≥ 0.9) chiếm phần lớn thời gian, năng lượng rung phải vượt 4.15 MeV, tức biên độ ≥ 0.91 $\ell_0$ ≈ 2.7κ. Có thể quét $\kappa$ (tức tỉ số $E_{\text{vt}}/\Delta E$) để tìm điểm mà pha Quá khứ chuyển hẳn sang "hạt". Có muốn chạy không?
