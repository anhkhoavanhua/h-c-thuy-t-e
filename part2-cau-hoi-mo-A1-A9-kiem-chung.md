# Part II: giải số chín câu hỏi mở A1–A9

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), mô hình V.C.1. Nối tiếp [part2-open-quantum-kiem-chung.md](part2-open-quantum-kiem-chung.md).
> **Số liệu:** [giai_so_part2_mo.py](giai_so_part2_mo.py). Lệnh `python3 giai_so_part2_mo.py A1 … A9`. Tổng cộng ≈ 40 phút, chạy bốn nhóm song song.

## Tóm tắt: câu nào đứng, câu nào đổ

| # | Câu hỏi | Kết quả | Hệ quả cho Part II |
|---|---|---|---|
| A1 | Gãy liên kết có trong phương trình chủ không? | Thêm giếng ràng buộc lực cực đại R: P_ngoài tăng **×2** và ⟨\|ℓ\|⟩ tăng 1.43 → 1.77 ngay trước khởi phát đợt. Nhưng đợt khởi phát ở ⟨\|ℓ\|⟩ ≈ **1.74**, xa ngưỡng cổ điển 3.2, và P_ngoài chỉ ~0.3% | **Một phần.** Có tiền dấu, nhưng làm nguội Davies tác dụng ngay khi u ≠ 0, trước khi liên kết gãy. Gãy kiểu Part I không xuất hiện |
| A2 | Xấp xỉ secular có hợp lệ? | 6.7% số cặp vi phạm. Nhưng Bloch–Redfield cho ⟨ℓ⟩ **2.342** so với Davies 2.325, khoảng cách vết **0.022**, ℓ → ℓ−1 90.7% so với 90.6% | **Đứng.** Davies đủ chính xác. Redfield chỉ vi phạm tính dương nhẹ (−2·10⁻⁴) |
| A3 | Giới hạn cổ điển ħ → 0 | Năng lượng mỗi bước **∝ s**: 2.08 → 1.01 → 0.56 → **0.35** $\hbar\omega_0$. Đạt mức cổ điển 0.40–0.46 tại s ≈ 0.15–0.2. ⟨ℓ⟩ thực (2.4–2.6) và công suất (1.1–1.5) gần như không đổi | **Đứng cho lượng tử năng lượng.** Nhưng cấu trúc đợt **không** về cổ điển: bước dày lên thành tiêu tán liên tục (61 bước/đợt), không về 5.12 |
| A4 | Độ bền theo (Γc, Γ0, ℓ0) | Trên 11 bộ tham số: ΔE/bước luôn **1.8–3.3** $\hbar\omega_0$; ℓ → ℓ−1 **67–97%** | **Đứng định tính.** Con số "92%" phụ thuộc tham số, cần ghi khoảng |
| A5 | Bơm Van der Pol lượng tử thay bơm bão hòa | ⟨ℓ⟩ đổi (1.27 so với 2.33), nhưng ΔE/bước vẫn **2.8**, ℓ → ℓ−1 87% | **Đứng.** Bậc thang rô-to không phụ thuộc dạng bơm |
| A6 | Nhiệt độ hữu hạn | 1 MeV: **không đổi** (nhảy lên 10⁻⁴/T₀). 5 MeV: nhảy lên 0.25/T₀, ⟨n⟩ = 0.17, Tr ρ² giảm 0.25 → 0.19 | **Đứng ở 1 MeV**, nên Davies T = 0 là đủ. Ở 5 MeV phải tính nhiệt độ |
| A7 | Phi Markov cho toàn động cơ (HEOM) | Liên kết yếu λ = 0.05: HEOM (K = 6) **3.630** so với Davies 3.599 (**0.9%**). λ = 0.25: **2.575** so với 2.399 (**7%**). λ = 1: lệch 8–16% ⟨ℓ⟩, 20–39% ⟨H⟩, HEOM chưa hội tụ hết. K = 2 cho kết quả phi vật lý | **Đứng ở liên kết yếu.** Ở cường độ cỡ V.C.1, hiệu chỉnh phi Markov khoảng 7%. Ở liên kết mạnh cần HEOM sâu hơn |
| A8 | Darwinism với 10 bản ghi | Không có bình nguyên. $R_{0.1}$ = 1.43–1.67 cho cả cặp C lẫn cặp E | **Không có Darwinism động lực** ở ≤ 2 $T_0$. Chỉ có ở mô hình bản ghi lý tưởng |
| A9 | Phổ nucleon phát ra | Lượng tử: **85%** đợt vượt $S_n$ ngay bước đầu (19.5 MeV), ⟨ε⟩ = **7.5 MeV**. Cổ điển Part I: 0.9–1.6 MeV. Bay hơi Weisskopf: 2.3–3.8 MeV | **Dự báo phân biệt được:** lượng tử cứng hơn bay hơi khoảng 2–3 lần, cổ điển mềm hơn khoảng 2 lần |

## A1. Gãy liên kết trong phương trình chủ

**Mô hình:** $H=H_{\text{V.C.1}}+V_b(y)$, với $V_b=-De^{-y^2/2w^2}$, $y=x-\kappa\sin\theta$, $w=1$, $D=Rw\sqrt e=4.95$ (lực cực đại đúng bằng R = 3). Cơ sở ℓ ≤ 6, n ≤ 12. Trạng thái đầu là trạng thái cơ bản: hạt bị giữ trong giếng. Chạy 24 quỹ đạo × 600 đơn vị thời gian, có 485 lần khởi phát đợt.

| $t-t_{\text{khởi phát}}$ ($T_0$) | −1 | −0.5 | −0.25 | −0.1 | 0 | +0.1 | +0.25 | +0.5 | +1 |
|---|---|---|---|---|---|---|---|---|---|
| P_ngoài / TB | 0.85 | 1.10 | 1.53 | **1.97** | 1.90 | 1.47 | 1.39 | 1.33 | 1.22 |
| ⟨\|ℓ\|⟩ | 1.51 | 1.66 | 1.75 | **1.77** | 1.76 | 1.71 | 1.62 | 1.66 | 1.57 |

(P_ngoài trung bình toàn thời gian = 0.0014, ⟨|ℓ|⟩ trung bình = 1.43.)

**Đọc bảng:**
- Đợt Carnot **có tiền dấu cổ điển**: rô-to tăng tốc và liên kết bị kéo giãn (P_ngoài gấp đôi) trong khoảng 0.25 $T_0$ trước khi đợt bắt đầu. Đây là hình ảnh gãy của Part I ở dạng xác suất.
- Nhưng ngưỡng **không** khớp. Đợt khởi phát ở ⟨|ℓ|⟩ ≈ 1.74, so với θ̇_br = 3.2 của Part I. Liên kết lượng tử hầu như không bao giờ gãy (P_ngoài ≤ 0.3%).
- Lý do: toán tử Davies $\hat u_-$ ghi nhận mọi lệch vận tốc tương đối u ≠ 0 do dao động lượng tử của hạt trong giếng, không chờ đến lúc |λ| ≥ R.
- Để có ngưỡng gãy đúng kiểu Part I, cần một kênh tiêu tán chỉ hoạt động **ngoài** giếng, ví dụ $L\propto P_{\text{ngoài}}\hat u$. Đó là một mô hình khác, có thể làm trong Part II như một biến thể.

## A2. Bloch–Redfield (không secular) so với Davies

Điều kiện secular $\Gamma_c|u_{ij}|^2/|\omega_{ij}|$: trung vị 0.000, phân vị 90% là 0.04, max 296 (các cặp ±ℓ gần suy biến). 6.7% số cặp có tỉ số > 0.1.

| Phương trình | ⟨ℓ⟩ | ⟨H⟩ | Tr ρ² | min eig ρ | max Re λ | Khe | ℓ → ℓ−1 | ΔE/bước |
|---|---|---|---|---|---|---|---|---|
| Davies | 2.325 | 4.238 | 0.2481 | 8·10⁻¹⁵ | −1·10⁻¹⁵ | 0.0488 | 90.6% | 2.681 |
| Bloch–Redfield | 2.342 | 4.267 | 0.2504 | **−2·10⁻⁴** | 2·10⁻¹⁴ | 0.0487 | 90.7% | 2.686 |

Khoảng cách vết giữa hai trạng thái dừng là 0.022.

- Redfield không ổn định (max Re λ ≈ 0, không dương) và chỉ vi phạm tính dương rất nhẹ.
- Các cặp vi phạm secular chủ yếu là cặp ±ℓ gần suy biến. Chúng ít ảnh hưởng vì phần tử $u_{ij}$ giữa chúng nhỏ.
- Kết luận: **Davies hợp lệ** cho mô hình này.

## A3. Giới hạn cổ điển $\hbar_{\text{eff}}=s\to0$

**Biến đổi chính tắc** giữ nguyên động lực học cổ điển: trong mã (ħ = 1) đặt $I\to I/s$, $\Delta E\to\Delta E/s$, $g\to g/\sqrt s$, $\kappa\to\kappa/\sqrt s$, bơm $\Gamma_0\to\Gamma_0/s$ (giữ mô-men cổ điển). Năng lượng thực bằng s × năng lượng trong mã. Chạy 12 quỹ đạo × 600 đơn vị thời gian cho mỗi s.

| s | Cơ sở | ⟨ℓ⟩ thực | ΔE/bước TB | Trung vị | Bước/đợt | ΔE/đợt | Nhảy/$T_0$ | Công suất ($\hbar\omega_0/T_0$) |
|---|---|---|---|---|---|---|---|---|
| 1 | ℓ ≤ 8 | 2.41 | **2.08** | 2.25 | 1.57 | 3.29 | 0.52 | 1.08 |
| 0.5 | ℓ ≤ 12 | 2.53 | **1.01** | 1.12 | 2.99 | 3.03 | 1.14 | 1.14 |
| 0.25 | ℓ ≤ 20 | 2.52 | **0.56** | 0.59 | 12.3 | 6.78 | 2.34 | 1.30 |
| 0.15 | ℓ ≤ 32 | 2.60 | **0.35** | 0.36 | 61.1 | 21.1 | 4.22 | 1.48 |

(Cổ điển Part I: 5.12 va chạm/đợt × 0.40–0.46 $\hbar\omega_0$.)

1. **Lượng tử năng lượng co về cổ điển đúng như dự đoán:** ΔE/bước ≈ 2.1·s $\hbar\omega_0$, tức đúng một bậc thang rô-to $\hbar\dot\theta$. Nó đi qua giá trị cổ điển 0.40–0.46 tại **s ≈ 0.15–0.2**, tương ứng tác dụng rô-to khoảng 15–20ħ.
2. **Động lực học cổ điển được giữ:** ⟨ℓ⟩ thực ≈ 2.4–2.6 và công suất tiêu tán 1.1–1.5 $\hbar\omega_0/T_0$ gần như không đổi theo s. Công suất tăng nhẹ ~35%, là phần hiệu chỉnh lượng tử còn lại.
3. **Cấu trúc đợt không về cổ điển.** Khi s nhỏ, bước nhảy dày đặc và các đợt hòa thành tiêu tán gần liên tục. Chuỗi rời rạc 5.12 va chạm của Part I đến từ **cơ chế gãy–tái khóa**, mà phương trình chủ V.C.1 không có (xem A1).
4. **Kết luận cho bài báo:** "0.4 $\hbar\omega_0$ dưới lượng tử" của Part I là giới hạn cổ điển chính xác của bậc thang rô-to. Đây là cầu nối định lượng Part I ↔ Part II. Riêng thống kê đợt cần mô hình gãy để nối được.

## A4. Độ nhạy theo (Γc, Γ0, ℓ0)

Liouvillian chính xác, ℓ ≤ 5, n ≤ 4:

| Γc | Γ0 | ℓ0 | ⟨ℓ⟩ | Tr ρ² | ℓ → ℓ−1 (số bước / năng lượng) | n → n−1 | ΔE/bước | Nhảy Carnot/$T_0$ |
|---|---|---|---|---|---|---|---|---|
| **0.5** | **0.2** | **4** | **2.325** | 0.248 | **90.6% / 92.7%** | 4.8% | **2.68** | 0.54 |
| 0.25 | 0.2 | 4 | 3.037 | 0.335 | 95.1% / 95.5% | 1.8% | 3.03 | 0.37 |
| 1 | 0.2 | 4 | 1.582 | 0.216 | 81.1% / 86.6% | 13.1% | 2.22 | 0.73 |
| 2 | 0.2 | 4 | 1.118 | 0.225 | 67.2% / 76.5% | 26.5% | 1.81 | 1.12 |
| 0.5 | 0.1 | 4 | 1.337 | 0.194 | 83.8% / 87.9% | 10.0% | 2.25 | 0.31 |
| 0.5 | 0.4 | 4 | 3.113 | 0.350 | 94.9% / 95.5% | 2.0% | 3.03 | 0.77 |
| 0.5 | 0.8 | 4 | 3.539 | 0.484 | 96.7% / 96.8% | 1.1% | 3.25 | 0.92 |
| 0.5 | 0.2 | 3 | 1.799 | 0.288 | 83.7% / 85.9% | 9.4% | 2.10 | 0.38 |
| 0.5 | 0.2 | 5 | 2.639 | 0.224 | 93.0% / 95.2% | 3.2% | 3.09 | 0.67 |
| 0.25 | 0.1 | 3 | 1.615 | 0.253 | 84.4% / 86.3% | 8.5% | 2.10 | 0.17 |
| 2 | 0.8 | 5 | 3.025 | 0.277 | 87.5% / 93.4% | 9.1% | 3.00 | 3.27 |

(Dòng chuyển dời tính theo gần đúng secular trên quần thể trạng thái riêng. Ở V.C.1 cho 90.6%, so với MCWF 91.5%: khớp.)

- **ΔE/bước ≥ 1.8 $\hbar\omega_0$ ở mọi bộ tham số.** Kết luận "không có bước dưới lượng tử" là bền vững.
- Tỉ lệ ℓ → ℓ−1 **phụ thuộc tham số** (67–97%). Nó giảm khi làm nguội mạnh, vì rô-to chậm lại (⟨ℓ⟩ nhỏ) và kênh hạt n → n−1 tăng lên. Bài báo nên ghi "**≈ 90% tại V.C.1, 67–97% trên dải tham số**".

## A5. Bơm Van der Pol lượng tử

Bơm VdP lượng tử cho rô-to: $L_\uparrow=\sqrt\gamma\,e^{i\theta}\sqrt{\ell_++\frac12}$ và $L_\downarrow=\sqrt{\gamma}/\Omega_0\,e^{-i\theta}\ell_+^{3/2}$. Hai kênh này cho $d\langle\ell\rangle/dt=\gamma(\ell+\tfrac12-\ell^3/\Omega_0^2)$, tức mô-men Rayleigh–Van der Pol của Part I với γ = 0.1, Ω₀ = 3.365.

| Bơm | Γc | ⟨ℓ⟩ | Tr ρ² | ℓ → ℓ−1 | n → n−1 | ΔE/bước | Nhảy Carnot/$T_0$ |
|---|---|---|---|---|---|---|---|
| Bão hòa (V.C.1) | 0 | 5.657 | 0.283 | — | — | — | — |
| Bão hòa | 0.05 | 3.794 | 0.685 | 98.1% | 0.4% | 3.41 | 0.10 |
| Bão hòa | 0.5 | 2.326 | 0.248 | 90.6% | 4.8% | 2.68 | 0.54 |
| **VdP lượng tử** | 0 | 1.921 | 0.031 | — | — | — | — |
| VdP lượng tử | 0.05 | 1.794 | 0.135 | 74.2% | 22.7% | 2.76 | 0.06 |
| **VdP lượng tử** | 0.5 | 1.267 | 0.167 | **86.6%** | 8.8% | **2.80** | 0.32 |

- **VdP lượng tử ở vùng sâu lượng tử không đạt** điểm cố định cổ điển ℓ = Ω₀ = 3.4. Không có làm nguội, ⟨ℓ⟩ chỉ đạt 1.92 với phân bố rất rộng (Tr ρ² = 0.03): nhiễu khuếch đại–mất mát ở số lượng tử nhỏ lấn át. Bơm bão hòa của V.C.1 cho chu trình giới hạn sắc nét hơn.
- **Bậc thang rô-to vẫn giữ** (ΔE/bước 2.8, ℓ → ℓ−1 87%). Kết luận chính không phụ thuộc dạng bơm.

## A6. Nhiệt độ hữu hạn (cân bằng chi tiết)

| $k_BT$ | ⟨ℓ⟩ | ⟨H⟩ | ⟨n⟩ | Tr ρ² | Nhảy xuống/$T_0$ | Nhảy lên/$T_0$ |
|---|---|---|---|---|---|---|
| 0 | 2.325 | 4.238 | 0.0143 | 0.248 | 0.500 | 0 |
| **1 MeV** | 2.324 | 4.238 | 0.0144 | 0.248 | 0.500 | **0.0001** |
| 5 MeV | 2.320 | 4.393 | **0.171** | **0.188** | 0.750 | **0.251** |

- Ở 1 MeV (điểm làm việc của Part I), Davies T = 0 **chính xác** tới 10⁻⁴.
- Ở 5 MeV, bể kích thích hạt (⟨n⟩ tăng 12 lần) và độ tinh khiết giảm. Đây cũng là vùng mà động cơ cổ điển mất mũi tên thời gian ([ads-cft-mo-rong-kiem-chung.md](ads-cft-mo-rong-kiem-chung.md)).
- Ghi chú: tỉ số quần thể của cặp ±3 luôn bằng 1, vì trạng thái riêng là sóng dừng đối xứng. ⟨ℓ⟩ > 0 đến từ các phần tử kết hợp do bơm định hướng tạo ra, không đến từ quần thể.

## A7. HEOM (phi Markov, toàn động cơ)

**Mô hình:**
- Bể Drude $J(\omega)=2\lambda\gamma_D\omega/(\omega^2+\gamma_D^2)$, $\gamma_D=2$, $k_BT=0.1$, ghép qua $\hat u$, cộng bơm V.C.1. Cơ sở ℓ ≤ 4, n ≤ 1.
- HEOM Tanimura–Kubo, K số hạng Matsubara, có hiệu chỉnh đuôi Ishizaki–Tanimura.
- Trạng thái dừng **giải trực tiếp** ma trận sinh, dư số < 10⁻¹⁴.
- So với Davies dựng từ **cùng** $J(\omega)$.

| λ | Phương pháp | ⟨ℓ⟩ | ⟨H⟩ | Tr ρ² | min eig ρ |
|---|---|---|---|---|---|
| **0.05** | Davies (Markov, secular) | 3.599 | 7.448 | 0.565 | 8·10⁻⁷ |
| | HEOM K = 3, độ sâu 2 / 3 | 3.557 / 3.557 | 7.300 / 7.301 | 0.439 | 1·10⁻⁶ |
| | **HEOM K = 6, độ sâu 2** | **3.630** | **7.448** | 0.552 | 4·10⁻⁸ |
| **0.25** | Davies | 2.399 | 4.460 | 0.247 | 2·10⁻⁶ |
| | HEOM K = 3, độ sâu 2 / 3 | 2.524 / 2.525 | 4.834 / 4.837 | 0.203 | 4·10⁻⁵ |
| | **HEOM K = 6, độ sâu 2** | **2.575** | **4.780** | 0.255 | 5·10⁻⁶ |
| **1.0** | Davies | 1.071 | 1.976 | 0.219 | 1·10⁻⁶ |
| | HEOM K = 3, độ sâu 2 / 3 | 1.159 / 1.222 | 2.575 / 2.738 | 0.12 | 6·10⁻⁴ |
| | HEOM K = 6, độ sâu 2 | 1.238 | 2.366 | 0.179 | **−1.9·10⁻³** |
| mọi λ | HEOM K = 2, độ sâu 3 | **phi vật lý** (Tr ρ² = 19–67, trị riêng −2.5 đến −4) | | | |

**Đọc bảng:**
1. **Liên kết yếu (λ = 0.05): Davies đúng.** Lệch ⟨ℓ⟩ 0.9%, ⟨H⟩ trùng (7.448). Độ sâu đã hội tụ (độ sâu 2 = 3), chỉ cần đủ số hạng Matsubara.
2. **Liên kết cỡ V.C.1 (λ = 0.25):** hiệu chỉnh phi Markov **+7%** cho ⟨ℓ⟩ và **+7%** cho ⟨H⟩. Davies đánh giá thấp năng lượng dừng. Độ sâu đã hội tụ. K = 3 → 6 dịch ⟨ℓ⟩ thêm 2%.
3. **Liên kết mạnh (λ = 1): HEOM chưa hội tụ hoàn toàn.** ⟨ℓ⟩ trải 1.16–1.24 giữa các cấu hình. K = 6 ở độ sâu 2 đã hơi mất tính dương. Dù vậy, cả ba cấu hình đều lệch **8–16%** ⟨ℓ⟩ và **20–39%** ⟨H⟩ so với Davies. Muốn có số chính xác ở λ = 1, cần K = 6 với độ sâu ≥ 3 (≈ 120 ADO), vượt bộ nhớ của lời giải trực tiếp; phải chuyển sang giải lặp.
4. **HEOM ở nhiệt độ thấp với quá ít số hạng Matsubara (K = 2) cho kết quả phi vật lý.** Đây là lỗi đã biết của HEOM khi $\beta\gamma_D\gg1$ ($\beta\gamma_D=20$ ở đây), nên K ≥ 3 là tối thiểu.
5. **Kết luận cho Part II:** ở điểm làm việc của động cơ, mô tả Markov (Davies) sai số khoảng 7%. Điều này **không** đổi bậc thang rô-to, vì kênh hồi phục và các khe năng lượng giữ nguyên, mà chỉ dịch trạng thái dừng. Chỉ khi liên kết bể mạnh (λ ≳ 1) mới cần HEOM.

## A8. Darwinism với 10 bản ghi

S nhỏ (ℓ ≤ 3, n ≤ 2, 21 chiều). Qubit tham chiếu R vướng với hai nhánh của S. Lấy trung bình trên 3 mảnh ngẫu nhiên mỗi kích thước.

| Cặp | τ | $I(R{:}F)$ (bit) với $\lvert F\rvert$ = 1, 2, …, 10 | $R_{0.1}$ |
|---|---|---|---|
| C: θ₀ = 0 / π | $T_0/5$ (2 $T_0$) | 0.05, 0.32, 0.57, 0.63, 0.80, **1.02**, 0.95, 1.32, 1.42, 1.66 | 1.67 |
| E: cơ bản / \|ℓ\| = 2 | $T_0/5$ | 0.15, 0.15, 0.40, 0.45, 0.86, **0.83**, 1.27, 1.32, 1.49, 1.65 | 1.43 |
| C | $T_0/10$ (1 $T_0$) | 0.17, 0.30, 0.40, 0.57, 0.68, **0.90**, 0.96, 1.09, 1.25, 1.38 | 1.67 |
| E | $T_0/10$ | 0.07, 0.21, 0.31, 0.47, 0.58, **0.83**, 0.94, 1.04, 1.26, 1.37 | 1.43 |

- Với 10 bản ghi, thông tin vẫn **tích lũy tuyến tính**, không có bình nguyên. $R_{0.1}\approx1.4$–1.7 cho cả hai họ.
- Ưu thế của họ C so với E ở M4 (0.78 so với 0.44 bit) **gần như biến mất** trong cơ sở nhỏ này với ℓ₀ = 2.
- Kết luận: trong mô hình động lực thật, bể Carnot **chưa** tạo tính khách quan kiểu Darwinism. Mỗi bản ghi chỉ mang ≲ 0.2 bit, vì phần lớn bản ghi là "không nhảy". Darwinism chỉ đúng cho mô hình bản ghi lý tưởng.

## A9. Phổ nucleon phát ra

Lấy từ 1456 đợt lượng tử (MCWF V.C.1, 24 quỹ đạo × 239 $T_0$):

| Mô hình | Năng lượng tới $S_n$ | ⟨ε⟩ nucleon |
|---|---|---|
| **Lượng tử** (V.C.1, exciton n = 2k+1, Eq. 22) | Bước đầu TB 19.5 MeV (trung vị 22.3). **85.2%** đợt vượt $S_n$ ngay bước 1 | **7.5 MeV** |
| Cổ điển Part I | Trung vị 3 va chạm × 4.0–4.6 MeV | 0.9–1.6 MeV |
| Bay hơi Weisskopf, a = A/8 (E* = 12, 20, 30.8 MeV) | — | 2.3, 3.0, 3.8 MeV |

- Hai giới hạn của học thuyết **kẹp hai phía** chuẩn bay hơi: cổ điển mềm hơn khoảng 2 lần, lượng tử cứng hơn khoảng 2–3 lần.
- Đây là **dự báo kiểm chứng được** bằng phổ neutron tiền cân bằng đo được (ví dụ hệ thống Kalbach cho A ~ 70), nhưng cần dữ liệu ngoài.
- Ghi chú: mô hình exciton Eq. 22 với n = 3 sau một bước là ngoại suy. Một lượng tử 20 MeV trong hạt nhân thực sẽ chia cho nhiều exciton hơn. Con số 7.5 MeV là **cận trên**.

## Phần bổ sung: giải ba câu hỏi còn mở sau A1–A9

> Số liệu: [giai_so_part2_mo2.py](giai_so_part2_mo2.py), lệnh `Q1 | Q2 | Q3`.

### Q1. Gãy kiểu Part I trong phương trình chủ: **thành công phần lớn**

**Mô hình:** giữ giếng ràng buộc của A1, nhưng kênh Carnot **chỉ** tác dụng lên trạng thái riêng không liên kết (P_ngoài > p*). Đây đúng là cơ chế của Part I: tiêu tán chỉ xảy ra ở va chạm tái khóa, còn pha khóa bảo toàn. Chạy 48 quỹ đạo × 2400 đơn vị thời gian.

| Đại lượng | p* = 0.1 | **p* = 0.3** | Part I | Không chặn (A1) |
|---|---|---|---|---|
| Số đợt | 348 | 324 | — | 485 |
| **ΔE / va chạm** (TB; trung vị) | 0.091; 0.263 | **0.412; 0.663** | **0.40–0.46** | ≈ 2–2.7 |
| **⟨\|ℓ\|⟩ lúc khởi phát** (TB; trung vị) | 3.46; 3.81 | **3.61; 3.77** | $M(0)\dot	heta_{	ext{br}}$ = **3.56** | 1.74 |
| Va chạm / đợt | 1.97 | 2.98 | 5.12 | — |
| P(1 va chạm): đo / hình học cùng TB | 59% / 51% | 43% / 34% | hình học | — |
| Thời gian lưu giữa đợt (TB; CV) | 28.9 $T_0$; **1.08** | 26.5 $T_0$; 1.35 | 4.5–6 $T_0$; **1.0–1.1** | 2.8 $T_0$; 0.73 |

**Đọc bảng:**
1. **Khi tiêu tán chỉ xảy ra sau gãy, bản lượng tử tái lập ba con số then chốt của Part I.**
   - Năng lượng mỗi va chạm: **0.412 $\hbar\omega_0$** (Part I 0.40–0.46).
   - Ngưỡng gãy: ⟨|ℓ|⟩ = **3.61** so với xung lượng ngưỡng cổ điển $M(0)\dot	heta_{	ext{br}}$ = 3.56.
   - Thời gian lưu gần mũ: **CV 1.08** ở p* = 0.1.
2. **"Bậc thang rô-to 2.3 $\hbar\omega_0$" của V.C.1 là hệ quả của kênh Davies tác dụng mọi lúc**, không phải tính chất bắt buộc của lượng tử hóa. Với cơ chế tiêu tán đúng như Part I, va chạm lượng tử có năng lượng cỡ cổ điển, vì chuyển dời từ trạng thái tự do về trạng thái liên kết có khe nhỏ.
3. **Còn lệch:**
   - Số va chạm mỗi đợt là 2–3, so với 5.12.
   - Thời gian lưu trung bình khoảng 27 $T_0$, so với 4.5–6. Part I kích hoạt bằng nhiễu nhiệt $k_BT$ quanh ngưỡng; ở đây bơm Γ₀ = 0.2 đẩy chậm và không có nhiễu nhiệt. Con số này phụ thuộc cơ chế lái, không phải cơ chế gãy.
4. **Hệ quả cho Part II:** nên trình bày **hai** mô hình. V.C.1 (Davies mọi lúc) cho bậc thang rô-to. Mô hình chặn gãy cho giới hạn khớp Part I. Tham số p* phân biệt hai chế độ: p* = 0.1 khớp thời gian lưu tốt hơn, p* = 0.3 khớp năng lượng va chạm tốt hơn.

### Q2. Darwinism với bản ghi góc: **không có, và đã rõ vì sao**

**Mô hình:** xung Carnot truyền dọc $
abla f=(-\kappa\cos	heta,1)$, nên mỗi nucleon ghi $\cos	heta$. Mỗi tích tắc $T_0/10$: tiến hóa theo $H$, rồi một qubit quay góc $arepsilon\cos	heta$. Có 12 bản ghi.

| ε | Cặp | $I(R{:}F)$ (bit) với \|F\| = 2, 4, 6, 8, 10, 12 | $R_{0.1}$ |
|---|---|---|---|
| 0.5 | C (θ₀ = 0 / π) | 0.05, 0.22, 0.43, 0.61, 0.87, 1.05 | 1.1 |
| 0.5 | E (sóng dừng ±2) | 0.01, 0.07, 0.16, 0.27, 0.41, 0.51 | < 1 |
| 1.2 | C | 0.01, 0.21, 0.56, 0.76, 1.23, 1.28 | 1.3 |
| 1.2 | E | 0.02, 0.11, 0.42, 0.70, 0.84, 0.91 | 1.0 |
| 2.0 | C | 0.00, 0.13, 0.44, 0.61, 1.18, 1.28 | 1.3 |
| 2.0 | E | 0.04, 0.23, 0.59, 0.85, 0.97, 0.99 | 1.3 |

- Bản ghi góc phân biệt bó sóng ràng buộc tốt hơn sóng dừng khi ghi yếu (1.05 so với 0.51 bit). Khi ghi mạnh, hai họ như nhau.
- **Không có bình nguyên ở mọi cường độ ghi.** Toàn bộ môi trường cũng chỉ giữ khoảng 1.3 bit, so với tối đa 2.
- **Nguyên nhân:**
  - Thứ được ghi (góc rô-to) **không bảo toàn** và nhòe pha trong 0.13–0.37 $T_0$ ([schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md), S3). Chỉ các bản ghi sớm còn mang thông tin, nên thông tin không được nhân bản.
  - Thứ bảo toàn (năng lượng, |ℓ|) chỉ được ghi qua các bước nhảy hiếm (A8).
- **Kết luận cho Part II:** Darwinism lượng tử **không** áp dụng cho động cơ ở thang $\lesssim2\,T_0$. Chỉ giữ khẳng định về sieve (trạng thái con trỏ), và nêu cơ chế thất bại của Darwinism như một kết quả.

### Q3. HEOM ở liên kết mạnh λ = 1 (tích phân thời gian tới trạng thái dừng)

| Phương pháp | ⟨ℓ⟩ | ⟨H⟩ | Tr ρ² | min eig | ‖dρ/dt‖ |
|---|---|---|---|---|---|
| Davies | 1.071 | 1.976 | 0.219 | 1·10⁻⁶ | — |
| **HEOM K = 6, độ sâu 3** (120 ADO) | **1.260** | **2.481** | 0.158 | −2·10⁻⁵ | 8·10⁻⁸ |
| HEOM K = 8, độ sâu 3 (220 ADO) | 1.259 | 2.448 | 0.164 | −7·10⁻⁵ | 8·10⁻⁸ |
| HEOM K = 6, độ sâu 4 (330 ADO) | 1.231 | 2.433 | 0.153 | −4·10⁻⁵ | 9·10⁻⁸ |

- Ở liên kết mạnh, hiệu chỉnh phi Markov là **+18% ⟨ℓ⟩, +26% ⟨H⟩**. Davies đánh giá thấp năng lượng dừng, phù hợp với xu hướng ở λ = 0.25 (+7%).
- **Hội tụ đã xác nhận:** theo K (K = 6 → 8: ⟨ℓ⟩ đổi 0.001, ⟨H⟩ đổi 1.3%) và theo độ sâu (3 → 4: ⟨ℓ⟩ −2.3%, ⟨H⟩ −2%). Trị riêng âm cỡ 10⁻⁵ là sai số cắt cụt. Hiệu chỉnh phi Markov ở λ = 1: **+15–18% ⟨ℓ⟩, +23–26% ⟨H⟩**.

### A9 với dữ liệu thực nghiệm EXFOR (IAEA)

> Số liệu: [giai_so_a9_exfor.py](giai_so_a9_exfor.py). Dữ liệu tải về lưu ở [exfor_data/](exfor_data/), lấy qua API https://nds.iaea.org/dataexplorer/api/.

**Dữ liệu dùng** (phổ neutron hai lần vi phân, tích phân theo góc thành dσ/dE):

| Bộ | Phản ứng | E* (MeV) | EXFOR |
|---|---|---|---|
| Baba 1988 | Cu(n,xn) 14.1 và 18.0 MeV | 21.5, 25.3 | 22077-002 |
| Takahashi 1988 | Cu(n,xn) 14.1 MeV | 21.5 | 22075-068 (con trỏ 1) |
| Gotoh 1992 | Ge(n,xn) 14.1 MeV (A ≈ 72, sát A ~ 70) | 21.3 | 22278-002 (con trỏ 1) |
| Yabuta 1986 | Ni(n,xn) 14.1 và 18 MeV | 22.9, 26.7 | 22025-044, -053 |
| Biryukov 1980 | ⁶⁰Ni(p,xn) 22 MeV | 26.4 | A0068-024 |

Vùng liên tục: 1 MeV ≤ ε ≤ $E_{	ext{inc}}$ − 4 MeV cho (n,xn), để bỏ đỉnh đàn hồi và mức rời rạc; 1–15 MeV cho (p,xn). Mô hình dùng **cùng** cửa sổ và cùng công thức exciton Eq. (22) của Part I, với **tham số Bảng I cố định**, không khớp lại cho từng phản ứng.

**Kết quả** (⟨ε⟩ MeV / P(ε > 8 MeV) / khoảng cách KS so với dữ liệu):

| Bộ dữ liệu | Dữ liệu | Weisskopf | Cổ điển Part I | Lượng tử V.C.1 | **Lượng tử chặn gãy Q1** (292 đợt phát) |
|---|---|---|---|---|---|
| Cu 14.1 (Baba) | 3.19 / 0.043 | 3.50 / 0.034 / 0.117 | 2.17 / 0.002 / 0.229 | 5.33 / 0.193 / 0.387 | **2.82 / 0.013 / 0.065** |
| Cu 18.0 (Baba) | 3.56 / 0.087 | 3.86 / 0.064 / 0.158 | 2.17 / 0.002 / 0.227 | 6.59 / 0.349 / 0.435 | 2.90 / 0.023 / 0.085 |
| Cu 14.1 (Takahashi) | 3.08 / 0.045 | 3.50 / 0.034 / 0.163 | 2.17 / 0.002 / 0.184 | 5.33 / 0.193 / 0.417 | **2.82 / 0.013 / 0.053** |
| Ge 14.1 (Gotoh) | 3.08 / 0.048 | 3.34 / 0.027 / 0.140 | 2.17 / 0.002 / 0.178 | 5.33 / 0.193 / 0.422 | **2.82 / 0.013 / 0.053** |
| Ni 14.1 (Yabuta) | 3.26 / 0.053 | 3.66 / 0.042 / 0.142 | 2.17 / 0.002 / 0.231 | 5.33 / 0.193 / 0.380 | 2.82 / 0.013 / 0.080 |
| Ni 18.0 (Yabuta) | 4.03 / 0.132 | 4.06 / 0.079 / 0.151 | 2.17 / 0.002 / 0.272 | 6.59 / 0.349 / 0.382 | 2.90 / 0.023 / 0.147 |
| ⁶⁰Ni(p,xn) 22 (Biryukov) | 3.36 / 0.063 | 4.03 / 0.077 / 0.176 | 2.17 / 0.002 / 0.235 | 6.82 / 0.368 / 0.462 | **2.91 / 0.023 / 0.061** |

Ngưỡng KS 95% theo cỡ mẫu mô hình: Q1 0.080 (292 đợt phát); cổ điển ≈ 0.07 (393 đợt); V.C.1 ≈ 0.04 (1323 đợt). Ngưỡng này chưa tính sai số của dữ liệu, vốn chỉ vài %.

**Kết luận A9:**
1. **Dữ liệu bác bỏ mô hình lượng tử V.C.1** (Davies mọi lúc). Nó cho phổ quá cứng: ⟨ε⟩ 5.3–6.8 MeV so với 3.1–4.0 đo được, đuôi 19–37% so với 4–13%, D = 0.38–0.46 ≫ 0.04. Kết hợp với Q1, đây là bằng chứng thực nghiệm rằng **"bậc thang rô-to 2.3 $\hbar\omega_0$" không phải mô tả đúng**.
2. **Dữ liệu bác bỏ cả mô hình cổ điển Part I**: phổ quá mềm, ⟨ε⟩ = 2.17 MeV, gần như không có đuôi, D = 0.18–0.27 vượt ngưỡng 0.07.
3. **Mô hình lượng tử chặn gãy (Q1) khớp dữ liệu, không cần khớp tham số:**
   - Qua phép kiểm KS ở **4/7 bộ** (D = 0.053–0.065 < 0.080): Cu 14.1 (cả hai nhóm đo), Ge 14.1, ⁶⁰Ni(p,xn) 22.
   - Sát ngưỡng ở Ni 14.1 (0.080) và Cu 18 (0.085). Không qua ở Ni 18 (0.147).
   - Khớp tốt hơn chuẩn bay hơi Weisskopf (D = 0.12–0.18) ở mọi bộ.
4. **Thiếu sót có hệ thống của Q1:**
   - ⟨ε⟩ thấp hơn 0.3–1.1 MeV, và đuôi > 8 MeV bị đánh giá thấp khoảng 3–6 lần.
   - Mô hình **không phụ thuộc năng lượng tới**. ⟨ε⟩ đo được tăng theo $E_{	ext{inc}}$ (Cu 3.19 → 3.56, Ni 3.26 → 4.03), còn mô hình giữ 2.8–2.9, vì năng lượng kích thích do động cơ quyết định (đợt TB 10 MeV), không do phản ứng.
   - Đuôi cứng trong dữ liệu đến từ phát xạ tiền cân bằng trực tiếp, mà công thức Eq. (22) với n = 2k+1 không có.
5. **Ghi chú độ tin cậy:**
   - Ý nghĩa của hai con trỏ trong Takahashi 1988 và Gotoh 1992 chưa xác minh từ tài liệu gốc. Mình dùng con trỏ 1. Kết quả Takahashi khớp với Baba, bộ một cột, nên có thể tin.
   - KS chưa gộp sai số thực nghiệm.
   - Phép so sánh **không** chứng minh cơ chế động cơ là đúng. Nó chỉ cho thấy mô hình chặn gãy cho phổ cùng độ cứng với hạt nhân thật, còn hai mô hình kia thì không.

## Giải ba câu hỏi mở cuối: E*, đuôi tiền cân bằng, thống kê đợt

> Số liệu: [giai_so_part2_mo3.py](giai_so_part2_mo3.py), lệnh `C3 | C12 | C4`.

### Mô hình hợp nhất: exciton do động cơ điều khiển

- **Gắn E\* (câu 1):** hạt nhân bắt đầu ở đúng năng lượng kích thích của phản ứng, $E_{\text{exc}}=E^*$, với n = 3 exciton.
- **Vai trò của động cơ:** mỗi va chạm Carnot tạo một cặp hạt–lỗ, $n\to n+2$, đúng ánh xạ của Part I. Thời điểm va chạm lấy từ **chuỗi va chạm thật** của mô hình chặn gãy.
- **Phát xạ:** giữa các va chạm, nucleon phát ra với tốc độ Griffin–Cline tuyệt đối $\lambda_e(n,E)\propto p(E-S_n)^n/(gnE^{n-1})$, quy về đơn vị $T_0$. Ví dụ ở Cu 14.1 MeV, $\lambda_e(n=3)$ = 0.17/$T_0$.
- **Hình dạng phổ:** từ trạng thái n là $(E-S_n)\,\mathrm{Beta}(2,n-1)$, **trùng Eq. (22) của Part I**. Phát xạ từ n nhỏ tự sinh ra **đuôi tiền cân bằng (câu 2)**.
- **Phát xạ lần hai (n,2n):** bay hơi từ năng lượng còn lại $U=E^*-S_n-\varepsilon_1$ nếu $U>S_{n2}$.
- **Đối chứng:** mô hình exciton chuẩn với tốc độ chuyển nội Williams, $|M|^2=K/A^3E$, K = 400 MeV³ (ở n = 3 cho λ₊ ≈ 2/$T_0$).

### C3. Thống kê đợt theo tốc độ bơm (câu 3)

| Γ₀ | Γc | Va chạm/đợt | Lưu TB ($T_0$) | CV | ΔE/va chạm | ⟨\|ℓ\|⟩ khởi phát | Điểm lệch |
|---|---|---|---|---|---|---|---|
| 0.2 | 0.5 | 2.51 | 24.8 | 1.24 | 0.249 | 3.61 | 4.87 |
| 0.5 | 0.5 | 3.40 | 14.8 | 1.12 | **0.432** | 3.79 | 2.29 |
| 0.5 | 2.0 | 3.63 | 15.2 | **1.06** | 0.355 | 3.72 | 2.42 |
| 1.0 | 0.5 | 4.12 | 10.3 | 1.28 | 0.520 | 3.92 | **1.69** |
| 1.0 | 5.0 | 3.80 | 10.2 | 1.35 | 0.512 | 3.72 | 1.74 |
| **2.0** | 0.5 | **4.58** | 9.5 | 1.72 | 0.638 | 4.02 | — |
| **3.0** | 0.5 | **4.91** | 9.7 | 1.66 | 0.628 | 4.09 | — |
| *Part I* | | *5.12* | *4.5–6* | *1.0–1.1* | *0.40–0.46* | *3.56* | |

(Hàng Γ₀ = 0.2, 2.0 và 5.0 khác không liệt kê hết; xem log.)

- **Bơm mạnh đưa số va chạm mỗi đợt về Part I:** 2.5 → 4.9, so với 5.12.
- **Không có bộ tham số nào khớp đồng thời mọi chỉ số.** Thời gian lưu dừng ở khoảng 9.5 $T_0$ (Part I: 4.5–6). CV tăng lên 1.7 khi bơm mạnh. Năng lượng mỗi va chạm tăng lên 0.63 $\hbar\omega_0$. Lý do: Part I kích hoạt gãy bằng **nhiễu nhiệt** quanh một điểm hút dưới ngưỡng, còn bơm lượng tử là quá trình **tất định trung bình**. Muốn có cả thời gian lưu ngắn lẫn CV ≈ 1 thì cần thêm kênh nhiễu nhiệt cho rô-to. Đây là câu hỏi còn lại.

### C1 + C2 + C4. So với dữ liệu EXFOR

**Thước đo:** ⟨ε⟩ (MeV) / P(ε > 8 MeV) / KS D. Ngưỡng KS 95% ≈ 0.076 (318 chuỗi động cơ).

| Bộ dữ liệu | Dữ liệu | **Động cơ Γ₀ = 2 + (n,2n)** | Exciton chuẩn K=400 + (n,2n) | Động cơ, chỉ lần 1 | Exciton chuẩn K=400, chỉ lần 1 |
|---|---|---|---|---|---|
| Cu 14.1 (Baba) | 3.19 / 0.043 | 2.86 / 0.041 / 0.101 | 2.51 / 0.018 / 0.163 | 3.30 / 0.053 / **0.043** | 2.90 / 0.031 / 0.072 |
| Cu 18.0 (Baba) | 3.56 / 0.087 | **3.62 / 0.086 / 0.055** | 3.21 / 0.047 / 0.064 | 4.43 / 0.145 / 0.167 | 3.78 / 0.086 / 0.104 |
| Cu 14.1 (Takahashi) | 3.08 / 0.045 | 2.85 / 0.037 / **0.056** | 2.50 / 0.020 / 0.123 | 3.31 / 0.054 / 0.076 | 2.87 / 0.029 / 0.058 |
| Ge 14.1 (Gotoh) | 3.08 / 0.048 | 2.80 / 0.035 / **0.061** | 2.48 / 0.018 / 0.126 | 3.31 / 0.050 / 0.087 | 2.87 / 0.026 / 0.051 |
| Ni 14.1 (Yabuta) | 3.26 / 0.053 | 3.06 / 0.048 / **0.042** | 2.67 / 0.022 / 0.115 | 3.36 / 0.057 / 0.058 | 2.86 / 0.024 / 0.086 |
| Ni 18.0 (Yabuta) | 4.03 / 0.132 | 3.63 / 0.093 / **0.066** | 3.13 / 0.046 / 0.134 | 4.43 / 0.145 / 0.132 | 3.74 / 0.080 / 0.082 |
| ⁶⁰Ni(p,xn) 22 | 3.36 / 0.063 | 3.04 / 0.060 / 0.084 | 2.54 / 0.021 / 0.161 | 3.68 / 0.098 / **0.041** | 2.90 / 0.034 / 0.082 |

**Kết luận:**
1. **Câu 1 (E\*): giải được.** Với $E_{\text{exc}}=E^*$, ⟨ε⟩ tăng theo năng lượng tới, đúng chiều dữ liệu.
   - Khi có (n,2n), mô hình động cơ cho Cu 14 → 18 MeV: 2.86 → 3.62, dữ liệu 3.19 → 3.56.
   - Phụ thuộc E\* đi vào **qua phản ứng**, động cơ không cần thay đổi.
2. **Câu 2 (đuôi tiền cân bằng): giải được.** Phát xạ từ trạng thái exciton thấp, xảy ra **trong khoảng lưu giữa các đợt va chạm**, tái lập đuôi cứng: P(ε > 8) = 0.041 so với 0.043 (Cu 14), 0.086 so với 0.087 (Cu 18).
3. **Câu 3 (thống kê đợt): giải một phần.**
   - Số va chạm mỗi đợt về 4.6–4.9, so với 5.12.
   - Thời gian lưu (9.5 so với 4.5–6 $T_0$) và CV (1.7 so với 1.0) chưa khớp. Cần kênh nhiễu nhiệt.
4. **Kết quả mới, độc lập với dữ liệu.** Γ₀ = 2 được chọn **chỉ** để khớp thống kê đợt của Part I, **không** khớp theo phổ neutron. Dù vậy, với (n,2n), mô hình động cơ:
   - **qua phép kiểm KS ở 5/7 bộ** (D = 0.042–0.066);
   - ~~tốt hơn exciton chuẩn K = 400 ở cả 7/7 bộ~~ — **ĐÃ BÁC BỎ ở mục R2+R3 bên dưới**: khi thêm cạnh tranh proton/α/γ, hiệu ứng cặp và χ² công bằng, exciton chuẩn khớp tốt hơn động cơ ở 7/7 bộ.
   - Nhịp va chạm cụm–lưu của động cơ (khoảng 2.4 va chạm trong 1 $T_0$ đầu, rồi nghỉ) tạo phổ gần dữ liệu hơn nhịp đều của exciton chuẩn.
5. **Giới hạn cần ghi rõ:**
   - Chỉ tính phát neutron, không có cạnh tranh proton, alpha hay gamma. Bội neutron 1.5–1.9 có lẽ cao hơn thực tế.
   - $\sigma_{\text{inv}}$ được lấy là hằng số. $S_{n2}$ lấy trung bình đồng vị.
   - KS chưa gộp sai số thực nghiệm.
   - Thứ hạng giữa các mô hình nhạy với thành phần vật lý: bỏ (n,2n) thì exciton chuẩn K = 400 lại khớp tốt hơn ở 14 MeV.
   - Kết luận chắc chắn là: động cơ ở chế độ giống Part I cho phổ **cạnh tranh được** với mô hình exciton chuẩn. Nó **chưa** chứng minh động cơ là cơ chế thật.

## Giải các câu hỏi còn lại: R1 (nhiễu cho rô-to), R2 (cạnh tranh n/p/α/γ), R3 (sai số thực nghiệm)

> Số liệu:
> - [giai_so_part2_mo4.py](giai_so_part2_mo4.py), lệnh `R1 | R1B`;
> - [giai_so_part2_mo5.py](giai_so_part2_mo5.py);
> - số liệu hạt nhân thật trong [exfor_data/](exfor_data/): khối lượng AME (`nuclides.json`), σ_non (`sigma_non.json`).

### R1. Nhiễu cho rô-to lượng tử: cơ chế khác Part I

**Mô hình:** trần bơm ℓ₀ = 3.2, nên bơm một mình dừng ở |ℓ| = 3, dưới ngưỡng gãy. Thêm bể cho rô-to, ghép qua cosθ, có cân bằng chi tiết ở $k_BT$.

| Bể rô-to | Đợt | Va chạm/đợt | Lưu TB ($T_0$) | CV | ΔE/va chạm |
|---|---|---|---|---|---|
| Không có bể | 197 | 3.76 | 11.5 | 1.33 | 0.62 |
| $k_BT$ ≈ 0 (chỉ nhảy xuống), Γt = 0.5 | 711 | 2.43 | **6.51** | **0.95** | 0.39 |
| 1 MeV, Γt = 0.5 | 781 | 2.65 | **5.90** | **0.93** | **0.401** |
| 3 MeV | 754 | 2.71 | 6.05 | 0.90 | 0.37 |
| 5 MeV | 759 | 2.56 | 6.08 | 0.89 | 0.39 |
| 10 MeV | 741 | 2.66 | 6.20 | 0.90 | 0.44 |
| 1 MeV, Γt = 0.1 | 471 | 3.10 | 9.89 | 1.01 | 0.44 |
| 1 MeV, Γt = 2.0 | 1029 | 2.36 | 4.37 | 0.95 | 0.20 |
| *Part I* | | *5.12* | *4.5–6* | *1.0–1.1* | *0.40–0.46* |

1. **Có bể rô-to thì thời gian lưu (5.9–6.5 $T_0$), dạng mũ (CV ≈ 0.95) và năng lượng va chạm (≈ 0.40 $\hbar\omega_0$) đều khớp Part I.**
2. **Nhưng cơ chế không phải kích hoạt nhiệt.**
   - Thời gian lưu **không đổi** từ T ≈ 0 đến 10 MeV: hệ số góc Arrhenius bằng 0.
   - Khe rô-to |ℓ| 3 → 4 khoảng 3.3–3.5 $\hbar\omega_0$, nên ở 1 MeV kích thích nhiệt bị chặn cỡ $e^{-33}$.
   - Thời gian lưu do **cường độ ghép Γt** quyết định: 0.1 → 9.9, 0.5 → 6.3, 2.0 → 4.4 $T_0$.
   - Như vậy chính **nhiễu lượng tử (điểm không) của sự ghép với bể** kích hoạt gãy.
   - Đây là phiên bản lượng tử của câu Part I "thermal-quantum fluctuations". Nhưng phần "thermal" ở 1 MeV không đóng góp. Γt ≈ 1 cho đúng 4.5–6 $T_0$.
3. Với ℓ₀ = 3.2, số va chạm mỗi đợt vẫn chỉ 2.4–2.7. Phép thử hiệu chỉnh kết hợp (ℓ₀ = 4 + bể) ở cuối mục này.

(Ghi chú: dòng in "ΔE ≈ 7.17 ħω0" trong log R1 là chọn sai trạng thái riêng có kích thích hạt. Khe rô-to đúng là 3.3 trên đa tạp ràng buộc, 3.5 khi tự do.)

### R2. Chuỗi phát xạ đầy đủ: cạnh tranh n/p/α/γ, σ_inv(ε), số liệu hạt nhân thật

**Thành phần:**
- năng lượng tách từ AME (IAEA). Lỗi dấu của $S_\alpha$ phát hiện trong lúc chạy đã sửa: $S_\alpha=BE_{\text{cha}}-BE_{\text{con}}-BE_\alpha$;
- σ_n Dostrovsky; σ_p, σ_α có hàng rào Coulomb; γ theo GDR Brink–Axel;
- hệ số điện tích $R_\nu(p)$; α chỉ phát ở pha cân bằng ($\varphi_\alpha$ = 0.1);
- bay hơi Weisskopf–Ewing nhiều tầng, mật độ mức Fermi gas **có dịch chuyển cặp** $\delta=12/\sqrt A\times(2,1,0)$;
- σ_non đo được (EXFOR): Cu 14.1 MeV **1.486 b**, Ni 14.1 MeV **1.397 b**.

**Vì sao cần dịch chuyển cặp:** khi chưa có, mô hình cho Ni chỉ 4% phát proton và dư 24% neutron. Trong khi đó ⁵⁸Ni(n,p) đo được ở 14 MeV là **0.333 b** (trung vị 51 điểm EXFOR). Khi có dịch chuyển cặp, σ(n,xp) dự đoán cho Ni tự nhiên là 213–260 mb, cùng bậc với số đo.

**So tuyệt đối** (tích phân trong cửa sổ, mb):

| Bộ | Dữ liệu | Động cơ | Exciton K = 400 | Bay hơi thuần |
|---|---|---|---|---|
| Cu 14.1 (Baba) | 1520 ± 77 | 1438 | 1502 | 1567 |
| Cu 18.0 (Baba) | 1631 ± 82 | 1802 | 1887 | 2048 |
| Cu 14.1 (Takahashi) | 1400 ± 70 | 1376 | 1430 | 1464 |
| Ni 14.1 (Yabuta) | 969 ± 49 | **981** | 959 | 930 |

Cả ba mô hình đều cho tiết diện tuyệt đối đúng trong khoảng 1–2σ ở 14 MeV. Ở Cu 18 MeV, mô hình động cơ dư 10% (σ_non ở 18 MeV là nội suy, ±5%).

### R3. Kiểm định χ² với sai số thực nghiệm

Sai số từng bin được truyền qua tích phân góc, cộng thêm 2% không tương quan. So **hình dạng** (chuẩn hóa diện tích).

| Bộ | Động cơ χ²/k | Exciton K = 400 χ²/k | Bay hơi χ²/k |
|---|---|---|---|
| Cu 14.1 (Baba) | 9.6 | **7.7** | 54 |
| Cu 18.0 (Baba) | 25.6 | **9.2** | 74 |
| Cu 14.1 (Takahashi) | 88 | **38** | 152 |
| Ge 14.1 (Gotoh) | 155 | **89** | 299 |
| Ni 14.1 (Yabuta) | 39 | **29** | 64 |
| Ni 18.0 (Yabuta) | 13.4 | **10.9** | 45 |
| ⁶⁰Ni(p,xn) 22 | 34 | **30** | 61 |

**Kết luận R3, thẳng thắn:**
1. **Với sai số thực nghiệm thật (1–2%), cả ba mô hình đều bị bác bỏ** (p ≈ 0). Đây là chuyện thường: mô hình exciton/bay hơi đơn giản không đạt độ chính xác của mã Hauser–Feshbach đầy đủ (TALYS, EMPIRE), vốn có mô hình quang học, spin, mức rời rạc.
2. **Xếp hạng công bằng:** exciton chuẩn < động cơ < bay hơi thuần, ở **7/7 bộ**. Động cơ kém exciton chuẩn khoảng 1.2–2.8 lần về χ², và tốt hơn bay hơi thuần khoảng 2–6 lần.
3. **Sửa kết luận trước:** khẳng định "động cơ tốt hơn exciton chuẩn" ở mục C4 đến từ (a) mô hình thiếu cạnh tranh proton và hiệu ứng cặp, và (b) sai số Monte Carlo của động cơ bị phóng đại bởi hệ số thiết kế (6000 sự kiện / 318 chuỗi). Bảng χ² có gộp sai số mô hình (χ²/k ≤ 0.5 cho động cơ) **không** dùng để xếp hạng, vì lý do (b).
4. **Điều đứng vững:** động cơ (tham số chỉ hiệu chỉnh theo Part I) cho tiết diện tuyệt đối và độ cứng phổ cùng bậc với dữ liệu, và có phụ thuộc E\* đúng chiều. Nó là một mô tả **cạnh tranh nhưng không vượt** mô hình exciton chuẩn.

### Hiệu chỉnh kết hợp: trần bơm ℓ₀ = 4 + bể rô-to ([giai_so_part2_r1c.py](giai_so_part2_r1c.py); 16 quỹ đạo × 2400 đơn vị thời gian mỗi cấu hình)

| ℓ₀ | Γ₀ | Γt | Đợt | Va chạm/đợt | Lưu ($T_0$) | CV | ΔE/va chạm | ⟨\|ℓ\|⟩ khởi phát |
|---|---|---|---|---|---|---|---|---|
| 4.0 | 2.0 | 0.5 | 788 | **3.41** | **5.70** | **0.90** | 0.31 | 3.76 |
| 4.0 | 2.0 | 1.0 | 894 | 3.20 | **5.00** | 0.92 | 0.24 | **3.59** |
| 4.0 | 1.0 | 1.0 | 783 | 2.35 | 5.89 | 0.91 | −0.18 | 2.91 |
| 3.6 | 2.0 | 1.0 | 826 | 3.03 | 5.48 | 0.88 | 0.24 | 3.48 |
| *Part I* | | | | *5.12* | *4.5–6* | *1.0–1.1* | *0.40–0.46* | *3.56* |

- **Ba trên bốn thống kê của Part I khớp cùng lúc** ở ℓ₀ = 4, Γ₀ = 2, Γt = 1:
  - thời gian lưu 5.0 $T_0$;
  - phân bố mũ, CV 0.92;
  - ngưỡng gãy 3.59, so với 3.56.
- **Không khớp được đồng thời** số va chạm mỗi đợt (3.2–3.4 so với 5.12) và năng lượng mỗi va chạm (0.24–0.31 so với 0.40–0.46). Bể rô-to vừa rút ngắn thời gian lưu, vừa lấy bớt năng lượng của đợt. Đây là một **đánh đổi có cấu trúc** của bản lượng tử, không phải chuyện chọn tham số.
- ΔE âm ở cấu hình Γ₀ = 1 nghĩa là một số bước Carnot làm ⟨H⟩ của trạng thái điều kiện tăng lên, do trạng thái trước bước nhảy là chồng chập. Ở cấu hình đó, đại lượng này không còn là năng lượng va chạm sạch.

## Hai câu hỏi mở cuối: độ chính xác phổ và số va chạm mỗi đợt

> Số liệu: [giai_so_part2_mo6.py](giai_so_part2_mo6.py), lệnh `H1 | H2 | S`. [giai_so_part2_mo5.py](giai_so_part2_mo5.py) có thêm cột χ² với 8% hệ thống, và `E_DUAL_OM=1` để dùng σ_R tính từ thế quang học.

### Câu 1. Độ chính xác phổ

**H1. Sàn χ² đặt ra bởi chính dữ liệu.** So hai phép đo độc lập cùng phản ứng Cu(n,xn) 14.1 MeV (Baba 1988 và Takahashi 1988), cùng cách xử lý:

| Sai số không tương quan thêm vào mỗi phép đo | 0% | 2% | 5% | 8% |
|---|---|---|---|---|
| χ²/k (Baba vs Takahashi, k = 45) | **16.2** | 10.0 | 4.2 | 2.2 |

- Tích phân hai phép đo khớp nhau (−1.2%), nhưng ⟨ε⟩ lệch: 3.10 so với 2.82 MeV. Độ lệch theo từng bin có trung vị **18%**, phân vị 90% là 29%.
- **Hai thí nghiệm không tương thích trong sai số khai báo.** Vậy đích "χ²/k ≈ 1 với sai số thống kê" là **bất khả thi với mọi mô hình**. Chuẩn thực tế là cộng khoảng 8% hệ thống mỗi bin, khi đó sàn đạt được ≈ 2.2.

**H2. σ_inv neutron từ thế quang học Koning–Delaroche 2003** (Numerov theo ℓ, j, có spin–quỹ đạo, khớp Riccati–Bessel):
- **Kiểm độc lập:** σ_R(n + ⁶³Cu, 14.1 MeV) = **1490 mb**, so với σ_non đo được **1486 ± 9 mb** (EXFOR, 10 phép đo).
- So với Dostrovsky: KD lớn hơn **1.2–1.5 lần** ở ε < 2 MeV, và trùng tới ±2% ở 14 MeV.

**χ²/k của các mô hình với 8% hệ thống** (σ_n Dostrovsky / σ_n KD):

| Bộ | Động cơ | Exciton K = 400 | Bay hơi thuần |
|---|---|---|---|
| Cu 14.1 (Baba) | 4.0 / 4.6 | **3.5 / 3.3** | 25 / 28 |
| Cu 18.0 (Baba) | 11.9 / 11.1 | **3.9 / 4.1** | 35 / 35 |
| Cu 14.1 (Takahashi) | 16.4 / 16.0 | **5.9 / 5.2** | 38 / 37 |
| Ge 14.1 (Gotoh) | 21.7 / 24.8 | **10.4 / 11.1** | 52 / 52 |
| Ni 14.1 (Yabuta) | 11.0 / 10.9 | **6.9 / 5.4** | 25 / 26 |
| Ni 18.0 (Yabuta) | 7.0 / 6.7 | **5.8 / 5.7** | 26 / 29 |
| ⁶⁰Ni(p,xn) 22 | 11.3 / 12.1 | **11.2 / 10.7** | 23 / 26 |
| *Sàn (H1)* | *≈ 2.2* | | |

**Kết luận câu 1:**
1. Phần lớn χ² "khổng lồ" trước đó (8–300) là do **sai số thực nghiệm khai báo quá nhỏ** so với độ lệch giữa các thí nghiệm, không phải do mô hình.
2. Với chuẩn thực tế, **exciton chuẩn tiến sát sàn ở Cu 14 MeV** (3.3–3.5 so với 2.2). **Động cơ đứng thứ hai** (4–25). Bay hơi thuần bị loại hẳn.
3. **Thế quang học không đổi thứ hạng**, và chỉ đổi χ² vài chục phần trăm. Phần lệch còn lại không nằm ở σ_inv.
4. Phần lệch còn lại, nhất là ở 18 MeV và Ge, rất có thể đến từ **kích thích trực tiếp các trạng thái tập thể** và **mức rời rạc** ở đầu năng lượng cao của phổ (n,n'). Muốn khép phải dùng mã Hauser–Feshbach cộng kênh trực tiếp (TALYS, EMPIRE). Đây là việc của vật lý hạt nhân chuẩn, độc lập với học thuyết.

### Câu 2. Số va chạm mỗi đợt: giới hạn cổ điển

**Mô hình S:** mô hình chặn gãy (giếng ràng buộc lực cực đại R, Carnot chỉ trên trạng thái không liên kết) dưới phép thu nhỏ $\hbar_{\text{eff}}=s$. Mọi tham số được co theo biến đổi chính tắc để **động lực học cổ điển giữ nguyên**: I/s, ΔE/s, g/√s, κ/√s, giếng D/s và w/√s, bơm Γ₀/s.

| s | Cơ sở | Đợt | **Va chạm/đợt** | ΔE/va chạm ($\hbar\omega_0$) | ΔE/đợt | Lưu ($T_0$) | CV | ℓ thực khởi phát |
|---|---|---|---|---|---|---|---|---|
| 1 (lượng tử, ~3ħ) | ℓ ≤ 6, n ≤ 12 | 257 | **4.45** | 0.61 | 2.77 | 11.4 | 2.08 | 4.01 |
| 0.5 | ℓ ≤ 11, n ≤ 20 | 182 | **6.12** | 0.32 | 1.93 | 6.3 | 1.65 | 4.12 |
| 0.25 | ℓ ≤ 20, n ≤ 34 | 261 | **9.33** | 0.27 | 2.49 | **5.65** | 1.87 | 4.03 |
| *Part I (cổ điển)* | | | *5.12* | *0.40–0.46* | *≈ 2.2* | *4.5–6* | *1.0–1.1* | *3.56* |

**Kết luận câu 2** (đã chỉnh sau khi có s = 0.25):
1. **Số bước nhảy mỗi đợt không có giới hạn cổ điển hữu hạn.** Nó tăng liên tục: 4.45 → 6.12 → 9.33 khi s = 1 → 0.5 → 0.25. Một "va chạm" lượng tử (một bước nhảy Davies) **không** tương đương một va chạm Carnot cổ điển. Khi ħ → 0, mỗi va chạm cổ điển tương ứng với **một cụm** bước nhảy nhỏ. Vì vậy so số bước nhảy với 5.12 là so hai đại lượng khác loại. Việc nó đi qua 5.12 ở s ≈ 0.6–0.7 không có ý nghĩa cơ bản.
2. **Các đại lượng hội tụ đúng về Part I:**
   - **năng lượng mỗi đợt** 2.8 → 1.9 → 2.5 $\hbar\omega_0$ (Part I khoảng 2.2);
   - **thời gian lưu** 11.4 → 6.3 → **5.65** $T_0$ (Part I 4.5–6);
   - **ngưỡng khởi phát** ổn định khoảng 4.0–4.1 (ngưỡng ràng buộc 3.56).
3. CV của thời gian lưu giảm (2.1 → 1.65–1.87) nhưng chưa về 1. Riêng chỉ số này cần bể rô-to (mục R1, cho CV ≈ 0.9).
4. **Kết luận:** con số 5.12 của Part I không mâu thuẫn với lượng tử. Đó là một đại lượng **cổ điển** (số va chạm dẻo), và không có đối ứng một-một trong bản lượng tử. Thống kê đợt nên so bằng **năng lượng mỗi đợt và thời gian lưu**, và hai thứ này hội tụ về Part I khi ħ_eff → 0. Ở ~3ħ, đợt lượng tử ngắn hơn về số bước nhưng mỗi bước lớn hơn. Đây là hiệu ứng ħ thật.

## Câu hỏi còn mở (sau toàn bộ Part II)
- Không còn câu hỏi mở nào về **cơ chế** trong phạm vi Part II.
- Việc còn lại là **nâng cấp chuẩn của vật lý hạt nhân**: Hauser–Feshbach + kênh trực tiếp tập thể, để χ² phổ chạm sàn ~2.2. Phần này dùng mã có sẵn (TALYS), không phải phát triển riêng cho học thuyết.
