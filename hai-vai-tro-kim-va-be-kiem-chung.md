# Hai vai trò: kim $\hat K$ cho thời gian, bể Carnot cho mũi tên — giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Part III. Nối tiếp [kim-chi-gio-page-wootters-kiem-chung.md](kim-chi-gio-page-wootters-kiem-chung.md).
> Cùng vũ trụ C ⊗ S ⊗ E: S có 36 chiều, τ = T0/25, mỗi tích tắc một ancilla {· không, C Carnot, B bơm}.
> **Số liệu:** [giai_so_hai_vai_tro.py](giai_so_hai_vai_tro.py): 4000 nhánh bản ghi × 1500 tích tắc (60 T0), chạy khoảng 1.8 phút.

## Tóm tắt

| Mệnh đề | Kiểm số | Kết quả |
|---|---|---|
| **V1a.** $\hat K$ tạo ra Lindblad với **mọi** bể | Tỉ số khi chia đôi τ là 2.00–2.04 ở cả ba mô hình: Carnot + bơm, chỉ Carnot, chỉ bơm | ✓ |
| **V1b.** Người đọc kim **không cần** đọc bể | Trung bình 4000 nhánh khớp Λᵏρ₀ tới 0.011–0.017, bằng mức nhiễu thống kê 0.016 | ✓ |
| **V2a.** Entropy bản ghi thô mang mũi tên **trên từng nhánh** | s_k giảm **0 / 6 000 000** bước. Mỗi tích tắc cộng ít nhất 0.027 nat. H tăng tuyến tính **5.23 nat/T0**, còn S(E) bão hòa ở 1.667 | ✓ |
| **V2a'.** Riêng bộ đếm N_C không đủ | H(N_C) chỉ tăng kiểu logarit: 3.00 nat ở 60 T0, so với 316 nat của toàn chuỗi | mũi tên nằm ở **toàn chuỗi bản ghi** |
| **V2b.** Hết tiêu tán thì hết mũi tên | Chỉ Carnot ở T = 0: h giảm từ 5.3 về **0.0000 nat/T0** sau 20 T0 | ✓ mũi tên cần động cơ còn chạy |
| **Thử chéo.** Bể và $\hat K$ cùng chiều? | Cùng nhánh: luôn cùng chiều. Hai nhánh độc lập: cần Δt ≈ **16 T0** để xếp đúng 99% | bể cho **chiều**, $\hat K$ cho **số đo** |

## Vai trò 1. Kim lý tưởng $\hat K$: thời gian và phương trình chuyển động

**1a. Lindblad không phụ thuộc bể là gì** (sai số vết so với e^{tℒ} tại t = 1 T0):

| Mô hình bể | τ = T0/25 | T0/50 | T0/100 | T0/200 | Tỉ số khi chia đôi τ |
|---|---|---|---|---|---|
| Carnot + bơm (E_dual) | 3.96e-2 | 1.94e-2 | 9.64e-3 | 4.82e-3 | 2.04, 2.01, 2.00 |
| Chỉ Carnot (T = 0) | 3.98e-2 | 1.98e-2 | 9.94e-3 | 4.98e-3 | 2.01, 1.99, 2.00 |
| Chỉ bơm | 2.75e-2 | 1.19e-2 | 5.65e-3 | 2.77e-3 | 2.30, 2.11, 2.04 |

Bể chỉ quyết định **ℒ là gì**. Việc có một tham số t để ℒ tác dụng theo thì đến từ $\hat K$.

**1b. Đọc kim mà không đọc bể.** Mỗi nhánh bản ghi là một quỹ đạo lượng tử thuần (Tr ρ² = 1). Trung bình trên các nhánh:

| k | ‖trung bình 4000 nhánh − Λᵏρ₀‖ | Nhiễu thống kê 1/√M | Tr(Λᵏρ₀)² |
|---|---|---|---|
| 25 (1 T0) | 0.011 | 0.016 | 0.377 |
| 125 (5 T0) | 0.017 | 0.016 | 0.247 |

Trạng thái ρ_S(t) mà kim $\hat K$ cho ra là **tổng trên mọi bản ghi**. Người đọc kim không cần biết bể ghi gì. Hai vai trò vì vậy tách được về mặt toán học.

## Vai trò 2. Bể Carnot: mũi tên qua entropy bản ghi thô

**Đại lượng dùng:** surprisal của chuỗi bản ghi trên từng nhánh, $s_k=-\ln p(r_1\dots r_k)$, đo trong cơ sở con trỏ {·, C, B}.
- Vì $p(r_1\dots r_{k+1})\le p(r_1\dots r_k)$, s_k **không thể giảm trên bất kỳ nhánh nào**.
- Trung bình ⟨s_k⟩ = H(r_1…r_k), tức entropy Shannon của bản ghi.

**2a. Số đo** (E_dual: Carnot + bơm):

| k | t/T0 | S(E) = S(ρ_S) | Σᵢ S(Eᵢ) | **H(chuỗi) = ⟨s_k⟩** | sd(s_k) | H(N_C) |
|---|---|---|---|---|---|---|
| 1 | 0.04 | 0.063 | 0.06 | 0.16 | 0.72 | 0.06 |
| 7 | 0.28 | 0.440 | 0.71 | 1.98 | 2.17 | 0.67 |
| 25 | 1 | 1.269 | 4.37 | 6.59 | 4.02 | 1.22 |
| 125 | 5 | 1.623 | 25.6 | 28.2 | 8.65 | 1.79 |
| 500 | 20 | 1.667 | 106.1 | 107.0 | 17.2 | 2.43 |
| 1500 | 60 | **1.667** (bão hòa) | 320.6 | **316.0** (tuyến tính) | 29.5 | **3.00** (kiểu log) |

- **Trên từng nhánh:** s_k giảm **0 lần** trong 6·10⁶ bước. Bước tăng nhỏ nhất là 0.027 nat, nên mỗi tích tắc đều để lại thông tin, kể cả tích tắc "không có gì xảy ra".
- **Tốc độ entropy** h = **5.23 nat/T0**, khớp tốc độ của Σᵢ S(Eᵢ) là 5.36. Hệ số khuếch tán D = 14.8 nat²/T0.
- **S(E) toàn cục** dừng ở 1.667 từ khoảng 20 T0 nên không thể làm mũi tên.
- **Riêng bộ đếm N_C:** H(N_C) chỉ tăng thêm 0.21 nat từ 40 đến 60 T0, kiểu ½ ln(phương sai). Phần lớn entropy nằm ở **thời điểm và loại** của từng bản ghi, không nằm ở số đếm.

**2b. Mũi tên tắt khi tiêu tán tắt.** Chỉ giữ kênh Carnot ở T = 0, xuất phát từ trạng thái riêng thứ 6:

| t | S(E) | H(chuỗi) | Tăng trong 1 T0 vừa qua |
|---|---|---|---|
| 1 T0 | 0.764 | 5.295 | 5.30 nat |
| 5 T0 | 0.021 | 6.714 | 0.008 nat |
| 20 T0 | 0.016 | 6.717 | 0.0000 nat |
| 100 T0 | 0.014 | 6.718 | 0.0000 nat |

- Khi hệ đã về trạng thái cơ bản, không còn bản ghi mới nên h → 0. Bản ghi cũ vẫn còn (H giữ ở 6.72), nhưng mũi tên **dừng**.
- Trong E_dual, bơm giữ động cơ ở trạng thái dừng không cân bằng, nên h = 5.23 nat/T0 được duy trì mãi.
- **Mũi tên thời gian của học thuyết là chữ ký của một động cơ đang chạy.**

## Phép thử chéo: mũi tên của bể có cùng chiều với $\hat K$?

**Cùng một nhánh:** $s_{k+\Delta k}\ge s_k$ với mọi Δk ≥ 0, nên bể **không bao giờ** chỉ ngược chiều $\hat K$.

**Hai nhánh độc lập** (hai người quan sát, mỗi người đọc bản ghi của mình), từ k₀ = 20 T0:

| Δt | 0.04 T0 | 1 T0 | 2 T0 | 5 T0 | 10 T0 | 12 T0 | **16 T0** | 24 T0 |
|---|---|---|---|---|---|---|---|---|
| P(s_b(k₀+Δk) > s_a(k₀)) | 0.503 | 0.585 | 0.665 | 0.852 | 0.974 | 0.985 | **0.998** | 1.000 |

**Dùng bể làm đồng hồ** (đọc $\hat t=s_k/h$ trên một nhánh): σ_t = sd(s_k)/h, so với yêu cầu 0.071 T0.

| t | σ_t |
|---|---|
| 1 T0 | 0.77 T0 |
| 5 T0 | 1.65 T0 |
| 20 T0 | 3.28 T0 |
| 60 T0 | 5.64 T0 |

- Bể cho **chiều** tuyệt đối trên mỗi nhánh. Nhưng **số đo** thời gian từ bể thô 10–80 lần so với yêu cầu, và càng tệ theo √t.
- Đây cùng bậc với kết quả TUR trong [nam-dong-ho-luong-tu-kiem-chung.md](nam-dong-ho-luong-tu-kiem-chung.md) và σ_t = 2–7.6 T0 của bộ đếm N_C.
- $\hat K$ cho số đo với σ_t = 0 theo cách dựng.

## Kết luận: chứng minh hai vai trò

1. **$\hat K$ → thời gian + động lực.** Chiếu lên ⟨k\| cho $\dot\rho=\mathcal L\rho$ với sai số O(τ), dù bể là gì (ba mô hình, tỉ số 2.00). Phép chiếu không cần đọc bản ghi: trung bình các nhánh bằng Λᵏρ₀ trong phạm vi nhiễu thống kê.
2. **Bể Carnot → mũi tên.** Entropy bản ghi thô $s_k=-\ln p(r_1\dots r_k)$ không giảm trên **từng** nhánh (0/6·10⁶). Trung bình của nó tăng tuyến tính h = 5.2 nat/T0. Trong khi đó S(E) toàn cục bão hòa ở 1.667, và riêng N_C chỉ tăng kiểu log.
3. **Mũi tên cần động cơ chạy.** Tắt bơm thì hệ về trạng thái cơ bản và h → 0.
4. **Hai vai trò nhất quán.** Bể không bao giờ ngược chiều $\hat K$. Nhưng bể chỉ đo được thời gian với σ_t ≈ 0.8–5.6 T0, và hai người quan sát cần cách nhau khoảng 16 T0 mới xếp đúng 99%. **Chiều** đến từ bể, **số đo** đến từ $\hat K$.

**Câu tiếng Anh cho bản thảo:**

> The two roles separate cleanly. The ideal pulse-train clock K̂ supplies the time parameter: conditioning on it yields the Lindblad equation with O(τ) error for every bath model tested (convergence ratio 2.00), and it does so without reading the environment records (the branch average reproduces Λᵏρ₀ within statistical noise). The Carnot bath supplies the arrow: the record surprisal −ln p(r₁…r_k) never decreases on any branch (0 of 6×10⁶ steps) and its mean grows linearly at h = 5.2 nat/T₀, whereas the global environment entropy saturates at 1.67 nat and the Carnot count alone grows only logarithmically. The arrow vanishes when dissipation stops (h → 0 for T = 0 relaxation), and as a clock the bath resolves time only to σ_t ≈ 0.8–5.6 T₀: the bath fixes the direction, K̂ the measure.

## Câu hỏi còn mở
- $\hat K$ vẫn là đồng hồ lý tưởng hình thức. Một đồng hồ vật lý đạt σ_t ≤ 0.071 T0 cần năng lượng vượt ngân sách hạt nhân (kết quả TUR, Salecker–Wigner). Học thuyết cần nói rõ $\hat K$ là đồng hồ **ngoài** hạt nhân, ví dụ chuỗi xung laser/gamma của phòng thí nghiệm, hay chỉ là tham số toán học.
