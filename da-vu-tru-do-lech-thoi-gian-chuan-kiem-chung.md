# Độ lệch thời gian chuẩn của đa vũ trụ E_dual

> **Câu hỏi:** đa vũ trụ gồm nhiều vũ trụ, mỗi vũ trụ có một động cơ E_dual với E năng lượng riêng. Tính một **con số chuẩn** cho độ lệch thời gian giữa chúng, mà vẫn tuân theo Wheeler–DeWitt.
> **Nối tiếp:** [hai-vu-tru-e-dual-xung-dot-kiem-chung.md](hai-vu-tru-e-dual-xung-dot-kiem-chung.md). Ở đó, năng lượng khác nhau làm nhịp thời gian khác nhau (×3.1 / ×1.7), và ghép ở tầng chân không kéo các nhịp lại gần nhau.
> **Số liệu:** [giai_so_da_vu_tru_do_lech_thoi_gian.py](giai_so_da_vu_tru_do_lech_thoi_gian.py) `[A|B|C]`, kết quả thô [da_vu_tru_ket_qua.txt](da_vu_tru_ket_qua.txt). B chạy khoảng 4 giây; A và C cùng khoảng 8 phút.

## Cách đo giữ đúng Wheeler–DeWitt

- Trạng thái chung của đa vũ trụ là **dừng**: toàn cục tĩnh. Không có đồng hồ ngoài nào.
- **Nhịp thời gian của vũ trụ i** là $r_i$, số bản ghi của chính nó trong một đơn vị chung. Chỉ các **tỉ số** $r_i/r_j$ mới có nghĩa; đó là tương quan nội tại, không phải số đọc của một đồng hồ ngoài.
- **Độ lệch thời gian chuẩn** được định nghĩa là
  $$\sigma_t=\operatorname{std}_i\ln\frac{r_i}{\langle r\rangle}$$
  tức độ lệch logarit so với nhịp trung bình của chính đa vũ trụ. Đây là thước đo **tương đối thuần túy**, nên tương thích với WdW.
- **Phân bố E năng lượng** giữa các vũ trụ thì học thuyết **không** cho sẵn. Con số chuẩn dưới đây dùng **cửa sổ vận hành** mà chính bài báo cho động cơ: tỉ số lái f ∈ (0.90, 1.139), với f phân bố đều trong cửa sổ.

## A. Con số chuẩn (cổ điển, Part I, cửa sổ vận hành của bài báo)

Mô phỏng Langevin đầy đủ ở k_BT = 1 MeV, mỗi f dùng 4 quỹ đạo × 100–120 T0.

| f | 0.900 | 0.950 | 1.000 | 1.050 | 1.100 | 1.130 | 1.139 |
|---|---|---|---|---|---|---|---|
| Nhịp r (bản ghi Carnot / T0) | 0.033 | 0.092 | 0.290 | 0.692 | 1.360 | 1.546 | 1.845 |
| Số lần gãy (đếm được) | 13 | 45 | 117 | 335 | 544 | 743 | 739 |

Đa vũ trụ có f đều trong cửa sổ:

| Đại lượng | Giá trị |
|---|---|
| **Độ lệch thời gian chuẩn σ_t** | **1.22 ± 0.06** (bootstrap Poisson) |
| Tức hệ số lệch điển hình $e^{\sigma_t}$ | **× 3.4** so với nhịp trung bình |
| Hệ số biến thiên CV = std(r)/⟨r⟩ | 0.89 |
| Tỉ số nhịp cực đại r_max / r_min | **57** |
| Đồng hồ theo lần gãy hay theo bản ghi Carnot | như nhau (1.214 so với 1.218) |

Trong cửa sổ vận hành của học thuyết, nhịp thời gian tăng gần **hàm mũ** theo f (gấp khoảng 3 lần mỗi khi f tăng 0.05). Vì vậy chỉ một dải lái hẹp (±12%) cũng đã đủ để các vũ trụ lệch thời gian cỡ **×3.4 điển hình** và **×57 cực đại**.

## B. Hằng số nội tại của một động cơ (lượng tử, Part III, chính xác)

| a (bơm) | 0.001 | 0.01 | 0.1 | 0.3 | **1 (E_dual)** | 3 | 10 | 100 | 1000 |
|---|---|---|---|---|---|---|---|---|---|
| ΔE (ħω0) | 0.20 | 0.57 | 0.96 | 1.73 | 4.26 | 6.66 | 7.92 | 8.51 | 8.53 |
| r (/T0) | 0.0007 | 0.005 | 0.044 | 0.150 | 0.494 | 0.824 | 0.997 | 1.057 | **1.079** |
| r / r(E_dual) | 0.001 | 0.010 | 0.089 | 0.304 | 1 | 1.67 | 2.02 | 2.14 | **2.18** |

- **Nhịp thời gian bão hòa.** r_max ≈ 1.08/T0, nên một vũ trụ **không thể chạy nhanh hơn E_dual quá ×2.18**, dù E năng lượng lớn đến đâu. Nửa bão hòa nằm ở a ≈ 1, tức E_dual nằm đúng ở "đầu gối" của đường cong.
- **Phía chậm thì không có chặn:** khi a → 0 thì r → 0 (đóng băng). Độ lệch thời gian giữa các vũ trụ vì thế **bất đối xứng**: nhanh hơn thì có trần (×2.18), chậm hơn thì không có sàn.

## C. Đa vũ trụ ghép ở tầng chân không: độ lệch co lại

Dùng 12 động cơ lượng tử với a log-đều trong [0.3, 3]. Mỗi động cơ tiêu tán thêm sang một "biển chung" với cường độ κ, và nhận lại bơm tỉ lệ với κ·J̄, trong đó J̄ là nhịp trung bình được giải tự hợp.

| κ | σ_t | r_max / r_min | ⟨r⟩ (/T0) |
|---|---|---|---|
| 0 (tách rời) | **0.548** | 5.5 | 0.48 |
| 0.3 | 0.451 | 4.0 | 0.60 |
| 1 | 0.216 | 1.95 | 1.05 |
| 3 | 0.021 | 1.07 | 2.96 |
| 10 | 0.002 | 1.01 | 8.46 |

- **Ghép đồng bộ hóa thời gian.** Khi κ ≳ 1, σ_t giảm gần như **κ⁻²** (0.216 → 0.021 → 0.002). Đa vũ trụ ghép mạnh có **một nhịp thời gian chung**.
- **Cảnh báo:** ⟨r⟩ tăng theo κ vì biển chung trong mô hình trường trung bình **bơm thêm năng lượng** (κ·J̄), tức kênh trao đổi không bảo toàn năng lượng. Xu hướng σ_t giảm là đáng tin; giá trị tuyệt đối của ⟨r⟩ ở κ lớn thì không.

## Kết luận: con số chuẩn

| Đa vũ trụ | Độ lệch thời gian chuẩn σ_t | Nghĩa |
|---|---|---|
| **Các vũ trụ tách rời, E năng lượng trải trong cửa sổ vận hành của học thuyết** | **1.22 ± 0.06** | lệch điển hình ×3.4, cực đại ×57 |
| Các vũ trụ lượng tử tách rời, E năng lượng ±nửa thập phân quanh E_dual | 0.55 | lệch điển hình ×1.7 |
| Ghép ở tầng chân không với κ ≳ 3 | ≤ 0.02 | gần như một nhịp chung |
| Giới hạn nội tại | nhanh hơn tối đa ×2.18; chậm hơn không có chặn | |

**Con số chuẩn đề xuất: σ_t = 1.22 ± 0.06 (e-fold)** cho đa vũ trụ không tương tác, với E năng lượng phân bố đều trong cửa sổ vận hành f ∈ (0.90, 1.139) của bài báo. Con số này:
- là **tỉ số nội tại** giữa các đồng hồ bản ghi, không cần đồng hồ ngoài, nên **tương thích với Wheeler–DeWitt**: trạng thái chung vẫn tĩnh;
- **phụ thuộc phân bố E năng lượng** (đều trong cửa sổ là một giả định) và phụ thuộc cường độ ghép. Ghép càng mạnh thì σ_t càng nhỏ, gần như κ⁻².

**Câu tiếng Anh cho bản thảo:**

> Define the multiverse time dispersion relationally, σ_t = std ln(r_i/⟨r⟩), where r_i is the record rate of universe i; no external clock is involved and the joint state remains Wheeler–DeWitt stationary. For engines whose drive is spread uniformly over the paper's operating window f ∈ (0.90, 1.139), the full Langevin dynamics gives σ_t = 1.22 ± 0.06: clocks typically differ by a factor 3.4 from the multiverse mean and by up to 57 between the window edges. A single engine's clock rate saturates at 2.18 times the E_dual rate, but has no lower bound (r → 0 as the drive vanishes). Vacuum-level coupling synchronizes the multiverse: in a self-consistent mean-field model σ_t falls from 0.55 to 0.02 between κ = 0 and κ = 3, roughly as κ⁻² beyond κ ≈ 1.

## Giới hạn

- **Phân bố E năng lượng là giả định.** "Đều trong cửa sổ vận hành" là lựa chọn tự nhiên nhất từ bài báo, nhưng không phải định luật. Phân bố khác thì σ_t khác.
- Ở f = 0.90 chỉ có 13 lần gãy (sai số khoảng ±28%). Bootstrap Poisson đã tính sai số này vào ±0.06.
- Phần C là trường trung bình (mỗi động cơ chỉ thấy trung bình của các động cơ khác) và kênh trao đổi **không bảo toàn năng lượng**. Một mô hình bảo toàn năng lượng giữa các vũ trụ có thể cho tốc độ đồng bộ khác.
- "Vũ trụ" ở đây là một động cơ E_dual. Con số là của mô hình, không phải phép đo vũ trụ học.

## Câu hỏi còn mở

- Thay phân bố đều bằng một phân bố "tự nhiên" hơn, ví dụ phân bố cực đại entropy với E năng lượng trung bình cố định, rồi xem σ_t đổi bao nhiêu.
- Làm lại phần C với kênh trao đổi **bảo toàn năng lượng** (mỗi lượng tử một vũ trụ nhả ra thì một vũ trụ khác nhận đúng lượng đó), để có σ_t(κ) đáng tin ở κ lớn.
