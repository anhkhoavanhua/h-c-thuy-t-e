# Đồng hồ hạt nhân làm đồng hồ Page–Wootters trong vũ trụ $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I ($T_0=4.136\cdot10^{-22}$ s). Nối tiếp [schrodinger-pw-wdw-e7-tong-hop.md](schrodinger-pw-wdw-e7-tong-hop.md) và [paper_part3.tex](paper_part3.tex).
> **Đồng hồ hạt nhân thật:** chuyển dời đồng phân ²²⁹ᵐTh, $E_c$ = **8.355733 eV**, tần số ν ≈ 2.0204·10¹⁵ Hz, chu kỳ 1/ν = 4.95·10⁻¹⁶ s. Kích thích laser trực tiếp năm 2024: Tiedau et al., PRL 132, 182501; Zhang et al., Nature 633, 63.
> **Số liệu:** [giai_so_dong_ho_hat_nhan.py](giai_so_dong_ho_hat_nhan.py), chạy ≈ 2 phút trên 8 lõi.

## Câu hỏi

Trong Page–Wootters, đồng hồ là một hệ lượng tử **nằm trong** vũ trụ tĩnh. Thời gian của người quan sát là số đọc của nó. Muốn dùng đồng hồ hạt nhân cho học thuyết, nó phải **phân giải** được ba pha Quá khứ → Hiện tại → Tương lai của động cơ. Ba pha đó diễn ra ở thang $T_0$, tức 10⁻²² s.

## N1. Đồng hồ phải chính xác đến đâu?

Mô phỏng Langevin Part I (Bảng I, 1 MeV, 193 đợt): thời gian lưu Quá khứ có trung vị 2.98 $T_0$; đợt Hiện tại có trung vị 0.88 $T_0$, phân vị 10% là 0.22 $T_0$. Mỗi lần đọc đồng hồ có sai số Gauss $\sigma_t$:

| $\sigma_t$ ($T_0$) | $\sigma_t$ (s) | Xếp đúng cả 3 pha | Hiện tại sau Quá khứ | Tương lai sau Hiện tại |
|---|---|---|---|---|
| 0.03 | 1.3·10⁻²³ | **99.8%** | 100% | 99.8% |
| 0.18 | 7.4·10⁻²³ | 93.9% | 100% | 93.9% |
| 1 | 4.1·10⁻²² | 69.6% | 94.3% | 75.0% |
| 5.6 | 2.3·10⁻²¹ | 32.4% | 68.7% | 56.6% |
| 1000 | 4.1·10⁻¹⁹ | **16.6%** (đoán bừa: 1/6) | 50% | 50% |

- **Ngưỡng:** xếp đúng 90% cần $\sigma_t\le$ **0.27 $T_0$ = 1.1·10⁻²² s**; xếp đúng 99% cần $\sigma_t\le$ **0.071 $T_0$ = 2.9·10⁻²³ s**.
- Điểm nghẽn là pha **Hiện tại → Tương lai**, vì đợt gãy ngắn (khoảng 0.9 $T_0$).

## N2. Đồng hồ hạt nhân ²²⁹Th có đủ không? Giới hạn Mandelstam–Tamm

Một đồng hồ lượng tử có trải năng lượng ΔE thì không đọc thời gian chính xác hơn $\sigma_t\ge\hbar/(2\Delta E)$.

| Đại lượng | Giá trị |
|---|---|
| Một hạt nhân ²²⁹Th: $\sigma_t^{\min}=\hbar/E_c$ | 7.9·10⁻¹⁷ s = **1.9·10⁵ $T_0$** |
| Một chu kỳ dao động của Th | 4.95·10⁻¹⁶ s = **1.2·10⁶ $T_0$** |
| ΔE cần cho 90% / 99% | **3.0 MeV / 11.2 MeV** |
| Số hạt nhân Th **độc lập** cần ($\Delta E=\sqrt N E_c/2$) | 5.0·10¹¹ / 7.2·10¹² |
| Số hạt nhân Th **vướng víu GHZ** cần ($\Delta E=NE_c/2$) | 7.1·10⁵ / 2.7·10⁶ |

1. **Một hạt nhân Th chậm hơn động cơ khoảng một triệu lần.** Trong một tích tắc của đồng hồ Th, động cơ đã quay khoảng 10⁶ chu kỳ. Nó không thể phân biệt Quá khứ, Hiện tại hay Tương lai.
2. **Muốn dùng Th, phải gộp rất nhiều hạt nhân:** khoảng 7·10⁵ hạt nhân ở trạng thái GHZ, hoặc khoảng 5·10¹¹ hạt nhân độc lập. GHZ cỡ 10⁶ là ngoài khả năng hiện nay (GHZ lớn nhất cỡ 10²–10³ qubit). Đồng hồ khi đó còn phải tương tác với động cơ trong 10⁻²² s.
3. **Trải năng lượng cần có là 3–11 MeV, tức khoảng 0.3–1.1 ħω₀.** Đây đúng là **thang của chính dao động hạt nhân**, vùng chuyển dời gamma MeV, lớn gấp khoảng 10⁶ lần chuyển dời đồng phân Th (eV).
4. **Kết luận:** đồng hồ hạt nhân tự nhiên cho học thuyết là một **đồng hồ gamma** (chuyển dời MeV), không phải đồng hồ đồng phân Th. Nói cách khác, đồng hồ cùng thang với $\hbar\omega_0$ của động cơ. Điều này khớp với V.C.4 (cần ≥ 25 xung mỗi $T_0$) và với Part III, nơi đồng hồ chuỗi xung có $\tau=T_0/25$.

## N3. Page–Wootters: động cơ trông thế nào qua đồng hồ độ phân giải Δ?

**Thiết lập:** đồng hồ thang M = 40 mức, độ phân giải Δ. Đọc số chỉ thời gian của nó tương đương làm trơn bằng nhân Fejér $|K(\tau)|^2$. Hệ là động cơ lượng tử V.C.1, 4 quỹ đạo MCWF (một bản ghi môi trường mỗi quỹ đạo, mỗi quỹ đạo 600 đơn vị thời gian). Tín hiệu là ⟨|ℓ|⟩(t).

| Δ ($T_0$) | Δ (s) | Độ thấy được của đợt V(Δ) | Tương quan với tín hiệu gốc | ΔE đồng hồ (MeV) |
|---|---|---|---|---|
| 0.02 | 8.3·10⁻²⁴ | **0.997** | 0.999 | 40 |
| 0.1 | 4.1·10⁻²³ | 0.985 | 0.995 | 8.0 |
| 0.5 | 2.1·10⁻²² | 0.935 | 0.973 | 1.6 |
| 2 | 8.3·10⁻²² | 0.816 | 0.909 | 0.40 |
| 10 | 4.1·10⁻²¹ | 0.393 | 0.564 | 0.08 |
| 50 | 2.1·10⁻²⁰ | **0.018** | (nhiễu cửa sổ) | 0.02 |
| *Th, 1 hạt nhân* | *7.9·10⁻¹⁷* | *→ 0* | | *8·10⁻⁶* |

(Ở Δ = 50 $T_0$ nhân làm trơn dài gần bằng chuỗi dữ liệu, nên tương quan −0.50 là sai số cửa sổ, không có ý nghĩa vật lý.)

1. **Đồng hồ MeV thấy được động lực học**: V ≥ 0.93 khi Δ ≤ 0.5 $T_0$, tức ΔE ≳ 1.6 MeV.
2. **Đồng hồ eV (Th) chỉ thấy vũ trụ tĩnh.** Mọi dao động bị lấy trung bình, người quan sát chỉ thấy **trạng thái dừng** của động cơ. Đây là minh họa cụ thể của Wheeler–DeWitt: cùng một vũ trụ, "có thời gian" hay "tĩnh" tùy độ phân giải của đồng hồ bên trong nó.

## Cách dùng đúng đồng hồ hạt nhân để phát triển học thuyết

| Vai trò | Đồng hồ | Có làm được không |
|---|---|---|
| Đồng hồ Page–Wootters **phân giải ba pha** | Chuyển dời gamma MeV (ΔE ≈ 3–11 MeV), hoặc đồng hồ chuỗi xung τ ≤ $T_0$/25 | ✓ về nguyên lý (N1–N3) |
| Đồng hồ Page–Wootters bằng ²²⁹Th | Th: σ_t ≈ 1.9·10⁵ $T_0$ | ✗ chỉ thấy trạng thái dừng |
| **Phép đo thống kê dài hạn** bằng ²²⁹Th | Th đo tần số dừng cực chính xác (cỡ 10⁻¹⁸–10⁻¹⁹) | Khả thi, nhưng chỉ kiểm được **đại lượng dừng**: nếu học thuyết dự đoán chân không hạt nhân ($E_{\text{vt}}$) làm dịch tần số đồng phân Th, đồng hồ Th đo được độ dịch đó. Học thuyết hiện **chưa** có dự đoán định lượng này |

## Câu hỏi còn mở
- Học thuyết có dự đoán nào cho **độ dịch tần số đồng phân ²²⁹Th** do $E_{\text{vt}}$ (chân không hạt nhân) không? Nếu có, đó là cách duy nhất để đồng hồ Th (eV, rất chính xác nhưng rất chậm) kiểm học thuyết bằng thực nghiệm thật. Cần một mô hình ghép giữa động cơ và các nucleon của ²²⁹Th.
