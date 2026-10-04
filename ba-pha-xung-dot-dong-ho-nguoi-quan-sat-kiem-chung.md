# Ghép Quá khứ + Hiện tại + Tương lai: các xung đột, và người quan sát làm đồng hồ Page–Wootters

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf).
> **Mô hình:** bản lượng tử hiệu chỉnh tốt nhất của Part II:
> - giếng ràng buộc lực R; Carnot chỉ tác dụng lên trạng thái không liên kết (p\* = 0.3);
> - bơm Γ₀ = 2 (ℓ₀ = 4); bể rô-to Γt = 1 ở 1 MeV;
> - khớp Part I ở thời gian lưu 5.0 $T_0$, CV 0.92, ngưỡng 3.59 ([part2-cau-hoi-mo-A1-A9-kiem-chung.md](part2-cau-hoi-mo-A1-A9-kiem-chung.md)).
>
> **"Người quan sát"** là bể Carnot ([nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md)).
> **Số liệu:** [giai_so_ba_pha_xung_dot.py](giai_so_ba_pha_xung_dot.py): 16 nhánh × 318 $T_0$, chạy ≈ 1.5 phút trên 8 lõi.

Mỗi quỹ đạo MCWF là **một nhánh** của trạng thái lịch sử Page–Wootters, ứng với một bản ghi môi trường cụ thể. Pha được gán theo bản ghi Carnot:
- **Hiện tại:** khoảng chứa một đợt bản ghi (các bản ghi cách nhau < $T_0$);
- **Tương lai:** $T_0$ ngay sau đợt;
- **Quá khứ:** thời gian còn lại.

## A. Dòng thời gian ba pha ghép lại

| Pha | Tỉ lệ thời gian | Độ dài một lần (TB / trung vị) |
|---|---|---|
| Quá khứ | **71.8%** | 3.88 / 2.37 $T_0$ |
| Hiện tại | **9.8%** | 0.89 / 0.69 $T_0$ |
| Tương lai | **18.4%** | 1.00 / 1.00 $T_0$ |

**Ma trận chuyển pha** (tổng trên 16 nhánh):

| Từ \ Sang | Quá khứ | Hiện tại | Tương lai |
|---|---|---|---|
| Quá khứ | — | **559** | 374 |
| Hiện tại | 0 | — | **562** |
| Tương lai | **930** | 1 | — |

- **Chu trình Quá khứ → Hiện tại → Tương lai → Quá khứ luôn đúng chiều.** Không có lần nào từ Hiện tại hay Tương lai quay về Quá khứ mà bỏ qua bước. Chỉ **1** lần gãy lại ngay trong pha Tương lai.
- **374 lần Quá khứ → Tương lai trực tiếp (40% số đợt):** đó là các đợt chỉ có **một** bản ghi Carnot. Ở các đợt này, **"Hiện tại" co lại thành một khoảnh khắc**, tức một bước nhảy lượng tử, không có độ dài. Đây là đặc điểm lượng tử, khác Part I nơi đợt luôn kéo dài.

## B. Các xung đột có thể xảy ra

### B1. Chồng chập pha: trước khi đọc bản ghi, vũ trụ ở cả ba pha cùng lúc

Tại một thời điểm, tỉ lệ nhánh ở Quá khứ / Hiện tại / Tương lai là **0.72 / 0.10 / 0.18**.

- **99%** thời điểm có ít nhất 2 pha cùng hiện diện trong các nhánh; **78%** có đủ cả 3 pha.
- Khi chưa đọc bản ghi môi trường, câu hỏi "bây giờ là pha nào?" **không có một đáp án duy nhất**. Trạng thái của động cơ là hỗn hợp cả ba pha. Xung đột này chỉ được giải khi **đọc bản ghi**, tức điều kiện hóa trên một nhánh.

### B2. Xung đột giữa hai người quan sát

Bản ghi Carnot không được nhân bản dư thừa ($R_\delta\approx1.4$–1.7, mục A8). Hai người quan sát A, B đọc hai mảnh môi trường khác nhau: mỗi người nhận một bản ghi Carnot với xác suất p, độc lập nhau.

| p (phần bản ghi mỗi người thấy) | Thời gian A, B **gán pha khác nhau** | A: QK – B: HT | A: HT – B: TL | A: QK – B: TL | Số đợt A thấy / thật |
|---|---|---|---|---|---|
| 1.0 | **0%** | 0% | 0% | 0% | 1.00 |
| 0.9 | 3.6% | 0.5% | 1.0% | 2.1% | 0.97 |
| 0.7 | 9.1% | 1.4% | 2.0% | 5.7% | 0.89 |
| 0.5 | 12.6% | 2.1% | 2.1% | **8.4%** | 0.78 |
| 0.3 | **13.2%** | 2.0% | 1.4% | **9.8%** | 0.58 |

- Khi bản ghi không đầy đủ, hai người quan sát **bất đồng tới 13% thời gian** về việc đang ở pha nào.
- Loại xung đột chủ yếu là **"A: Quá khứ – B: Tương lai"**: một người đã thấy một đợt Carnot, người kia bỏ lỡ nó (chỉ thấy 58–78% số đợt). Họ **không mâu thuẫn về thứ tự** mà mâu thuẫn về việc **một sự kiện đã xảy ra hay chưa**.
- Đây là hệ quả số của việc thiếu Darwinism: tính khách quan của "Hiện tại" chỉ có khi bản ghi được nhân bản đủ.

### B3. Đồng hồ không phân biệt được hai mốc

Nếu dùng **mọi bản ghi** (Carnot, bơm, bể) làm tích tắc, thì **2.4%** cặp mốc pha liên tiếp (59/2410) **không có tích tắc nào ở giữa**. Trong các trường hợp này, đồng hồ đọc cùng một số cho hai pha khác nhau.

### B4. Đảo thứ tự (đã đo ở các mục trước)

- Đồng hồ có sai số σ_t: tỉ lệ xếp đúng ba pha giảm từ 99.8% (σ ≤ 0.03 $T_0$) xuống 16.7% (đoán bừa) khi σ ≫ $T_0$ ([dong-ho-hat-nhan-page-wootters-kiem-chung.md](dong-ho-hat-nhan-page-wootters-kiem-chung.md), N1).
- Ở nhiệt độ cao (5 MeV), động cơ cổ điển đọc "Tương lai trước Hiện tại" ở 54.5% số đợt ([dong-ho-photon-kiem-chung.md](dong-ho-photon-kiem-chung.md)).

## C. Người quan sát làm đồng hồ Page–Wootters

Đồng hồ là **bộ đếm bản ghi của bể Carnot**.

| Tiêu chí của đồng hồ PW | Chỉ bản ghi Carnot | Mọi bản ghi (Carnot + bơm + bể) | Đòi hỏi để phân giải ba pha |
|---|---|---|---|
| Đơn điệu (không đếm lùi) | **✓ có** | **✓ có** | ✓ |
| Tuần hoàn / vấn đề đếm vòng | **Không tuần hoàn** | **Không tuần hoàn** | phải không tuần hoàn |
| Tốc độ tích tắc | 0.61 / $T_0$ | 4.31 / $T_0$ | ≥ 25 / $T_0$ |
| Độ đều (CV khoảng giữa tích tắc) | 2.05 | 1.52 | ≈ 0 (đều) |
| Độ phân giải ở tích tắc thứ 5 / 20 / 80 | 62% / 44% / 21% | 40% / 21% / 15% | σ ≤ 0.07 $T_0$ |
| Độc lập với hệ: tương quan tốc độ tích tắc ↔ ⟨\|ℓ\|⟩ | −0.21 | **−0.63** | 0 (đồng hồ lý tưởng không tương tác) |
| ⟨\|ℓ\|⟩ ngay sau tích tắc | 3.04 ± 0.73 | 3.08 ± 0.65 | — |

**Kết luận phần C:**
1. **Người quan sát làm được đồng hồ Page–Wootters theo nghĩa quan hệ.** Nó **đơn điệu** và **không tuần hoàn**, tức tự giải được vấn đề đếm vòng mà đồng hồ rô-to, photon phân cực hay OAM đều mắc. Nó cũng **không thuận nghịch**, nên chính nó là mũi tên thời gian.
2. **Nó không phải đồng hồ PW lý tưởng:**
   - **chậm**: 0.6–4.3 tích tắc mỗi $T_0$, so với ≥ 25 cần có;
   - **không đều**: CV 1.5–2.1;
   - **phụ thuộc mạnh vào hệ**: tương quan tới −0.63, vì tốc độ tích tắc đổi theo trạng thái động cơ. Đây đúng là trường hợp "đồng hồ tương tác" trong Page–Wootters mở rộng, khi đồng hồ và hệ không tách nhau. Thời gian đọc được bị co giãn theo trạng thái hệ.
3. **Tích tắc rơi vào cùng một trạng thái**: ⟨|ℓ|⟩ ≈ 3.0–3.1 ngay sau mỗi tích tắc, tức sát ngưỡng gãy. Người quan sát là một **đồng hồ sự kiện**: nó tích tắc khi động cơ chạm ngưỡng, chứ không tích tắc đều theo thời gian. Thời gian của người quan sát vì vậy **được định nghĩa bởi chính các lần chuyển pha**.
4. **Độ phân giải cải thiện theo √k** (40% → 15% với mọi bản ghi). Đồng hồ này dùng tốt cho **khoảng thời gian dài** (thống kê nhiều chu trình), nhưng không phân giải được từng đợt.

## Bức tranh hợp nhất

| Tầng | Câu trả lời của số liệu |
|---|---|
| **Vũ trụ toàn cục** (Wheeler–DeWitt) | Tĩnh, thuần, không có pha (Part III: ràng buộc thỏa tới 4·10⁻¹⁵) |
| **Trước khi đọc bản ghi** | Hỗn hợp cả ba pha (B1: 99% thời điểm có ≥ 2 pha) |
| **Một người quan sát đọc bản ghi** | Một dòng thời gian Quá khứ → Hiện tại → Tương lai luôn đúng chiều (A); 40% đợt có "Hiện tại" tức thời |
| **Hai người quan sát** | Bất đồng 0–13% thời gian, tùy phần bản ghi mỗi người có (B2) |
| **Người quan sát làm đồng hồ** | Đơn điệu, không tuần hoàn, là mũi tên thời gian; nhưng chậm, không đều, phụ thuộc hệ (C) |

**Câu cho học thuyết:** trong $\mathcal E_{\text{dual}}$, người quan sát và đồng hồ là **cùng một thứ**: bể nhận năng lượng Carnot. Thời gian là số bản ghi nó đã nhận. Quá khứ là phần đã ghi, Hiện tại là lúc đang ghi, Tương lai là phần chưa ghi. Xung đột giữa các người quan sát đến từ việc mỗi người chỉ đọc được một phần bản ghi.

## Câu hỏi còn mở
- Muốn hai người quan sát **luôn đồng ý** (p → 1), bể phải nhân bản mỗi bản ghi Carnot vào nhiều mảnh môi trường. Học thuyết có cơ chế nào (ví dụ nhiều nucleon cùng nhận xung) để tăng độ dư thừa không? Mục A8 cho thấy với cơ chế hiện tại R ≈ 1.5.
