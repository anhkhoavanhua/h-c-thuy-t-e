# E năng lượng bên ngoài: phổ năng lượng mà động cơ gặp khi thoát khỏi vũ trụ của nó

> **Ý của học thuyết:** khi thoát khỏi một vũ trụ có "hàm tĩnh", ta gặp một **E năng lượng khác đang tồn tại ngoài đó**, gồm nhiều năng lượng khác nhau.
> **Cách đọc trong mô hình:**
> - "vũ trụ bên trong" là hệ S (rô-to + hạt);
> - "bên ngoài" là môi trường mà S chỉ chạm tới khi thoát khỏi 0: bể Carnot (nhận năng lượng ra) và nguồn bơm (cấp năng lượng vào);
> - **E năng lượng bên ngoài** là phổ các lượng tử năng lượng mà S trao đổi với môi trường đó.
>
> **Nối tiếp:** [e-nang-luong-thoat-ly-chan-khong-kiem-chung.md](e-nang-luong-thoat-ly-chan-khong-kiem-chung.md), [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md), và chuỗi tài liệu về lời giải frozen time.
> **Số liệu:** [giai_so_e_nang_luong_pho.py](giai_so_e_nang_luong_pho.py) (A, B, C; chạy khoảng 3 giây) và phần D tính trực tiếp từ `rl_data/chu_trinh_f.json`. Kết quả thô ở [e_nang_luong_pho_ket_qua.txt](e_nang_luong_pho_ket_qua.txt).

## A. Năng lượng bên trong khi đã thoát khỏi 0 (Part III, trạng thái dừng E_dual)

| Thành phần | Chân không (ħω0) | E_dual | Vượt chân không | MeV |
|---|---|---|---|---|
| Động năng rô-to ℓ²/2I | 0.209 | 3.679 | **+3.470** | **+34.7** |
| Thế pha −ΔE cosθ | −0.744 | 0.034 | +0.777 | +7.8 |
| Dao động tử (n+½)ħω0 | 0.501 | 0.515 | +0.014 | +0.1 |
| Ghép −g x sinθ | −0.002 | −0.001 | +0.000 | 0.0 |
| **Tổng** | −0.036 | 4.226 | **+4.262** | **+42.6** |

- Khoảng 81% năng lượng vượt chân không nằm ở **chuyển động quay**, 18% ở thế pha, và gần như **không** nằm ở dao động tử.
- Năng lượng điểm không ½ħω0 (5 MeV) của dao động tử có mặt cả ở chân không lẫn E_dual. Nó **không** tham gia vào việc thoát khỏi 0.

## B. Phổ E năng lượng bên ngoài (lượng tử, Part III)

Các vạch $\omega=E_j-E_i$ được tính trong gần đúng thế tục: tốc độ $R_{ij}=|\langle i|L|j\rangle|^2p_j$.

**Năng lượng Carnot nhả ra môi trường** (Hiện tại → Tương lai):
- **23 vạch** riêng biệt, trải từ 0 tới **6.0 ħω0 (60 MeV)**. Tổng 0.535 lượng tử mỗi T0, công suất 1.416 ħω0/T0.
- Gần đúng thế tục khớp Lindblad chính xác (1.370) **trong 3%**.

| Vạch ω (ħω0) | MeV | Tốc độ (/T0) | Công suất (ħω0/T0) | Phần công suất |
|---|---|---|---|---|
| **3.54** | 35.4 | 0.197 | 0.697 | **49%** |
| **2.47** | 24.7 | 0.107 | 0.265 | 19% |
| **2.45** | 24.5 | 0.101 | 0.246 | 17% |
| 1.21 | 12.1 | 0.044 | 0.054 | 4% |
| 1.75 | 17.5 | 0.026 | 0.046 | 3% |
| 3.68, 0.99, 2.60, 4.20, 2.69, … | 10–42 | ≤ 0.016 | ≤ 0.023 | còn lại |

Phân loại theo chuyển dời: **97%** công suất do bước rô-to (ℓ đổi), 2% do bước dao động tử (n đổi), 1% do cả hai cùng đổi.

**Năng lượng bơm cấp vào** (Quá khứ → Hiện tại): **62 vạch** trong khoảng −6.0 … +4.2 ħω0, cũng do ℓ chiếm 99%.
- **Cảnh báo:** với kênh bơm, gần đúng thế tục chỉ cho 0.346 ħω0/T0, so với 1.370 của Lindblad chính xác (thiếu khoảng 75%), và có cả vạch âm. Kênh bơm tạo **kết hợp** giữa các trạng thái riêng, nên phổ vạch của nó **không đáng tin**. Phổ Carnot thì đáng tin.

Như vậy, E năng lượng bên ngoài mà động cơ gặp là một **phổ vạch rời rạc nhiều thành phần**: khoảng 23 vạch nhả ra và khoảng 60 vạch nhận vào. Ba vạch chính 35, 25 và 25 MeV mang 85% công suất. Tất cả đều là **bậc thang lượng tử của rô-to** (khoảng cách mức ℓ), khớp với bước nhảy rô-to khoảng 2–3.5 ħω0 đã thấy ở [schrodinger-rang-buoc](schrodinger-rang-buoc-kiem-chung.md) và [part3-latex-so-lieu](part3-latex-so-lieu-kiem-chung.md).

## C. "Hàm toàn cục không cố định": khối vũ trụ lớn dần

Trong WdW, $\hat{\mathcal H}|\Psi\rangle\rangle=0$ là **định nghĩa** của tĩnh. Muốn |Ψ⟩⟩ tự "đổi theo thời gian" thì phải có một thời gian ngoài, trái với tiền đề, nên điều đó **không làm được**.

Thứ đo được là **khối vũ trụ lớn dần** (growing block): xếp các lịch sử $|\Psi_N\rangle\rangle$ theo số nhịp N. Nội dung của khối là số nhánh hiệu dụng $e^{H_N}$, với $H_N$ là entropy của chuỗi bản ghi (200 nhánh).

| Chế độ | H (nat) tại 1 T0 | 5 T0 | 10 T0 | 20 T0 | 40 T0 | Tốc độ tăng 20→40 T0 (nat/T0) |
|---|---|---|---|---|---|---|
| PW đóng | 0 | 0 | 0 | 0 | 0 | **0**: khối không lớn, một nhánh duy nhất |
| Chỉ Carnot (Quá khứ thắng) | 3.5 | 9.4 | 11.4 | 12.5 | 13.1 | 0.03: ngừng lớn |
| Chỉ bơm (Hiện tại thắng) | 2.4 | 7.3 | 8.8 | 10.1 | 11.4 | 0.07: ngừng lớn |
| Cân bằng nhiệt | 7.8 | 36 | 71 | 132 | 258 | 6.29: lớn, nhưng **không chiều** |
| **E_dual (ba pha xung đột)** | 6.9 | 29 | 55 | 108 | **213** | **5.24**: lớn đều, **có chiều** |

- Với E_dual, sau 40 T0 khối chứa khoảng $e^{213}\approx10^{92}$ nhánh hiệu dụng và lớn thêm khoảng **e^{5.24} ≈ 190 lần mỗi T0**.
- Khi một pha thắng, khối gần như ngừng lớn (≤ 13 nat): đóng băng.

**Nghĩa chính xác:** mỗi |Ψ_N⟩⟩ vẫn tĩnh. Nhưng **dãy** |Ψ_N⟩⟩ không cố định: nó lớn dần khi bản ghi tăng. Đó là "hàm không cố định" theo nghĩa khối lớn dần. Cái giá là chỉ số N (số nhịp, cũng là số bản ghi) đóng vai một **siêu thời gian**; điều này phải ghi rõ khi trích.

## D. Phổ E năng lượng bên ngoài ở bản cổ điển (Part I, Langevin đầy đủ)

Năng lượng mỗi va chạm Carnot nhả ra ngoài (k_BT = 1 MeV):

| f | Số va chạm | TB (MeV) | Trung vị | 5–95% (ħω0) | Max | < 0.1 ħω0 | 0.5–1.0 ħω0 |
|---|---|---|---|---|---|---|---|
| 0.95 | 44 | 4.3 | 0.47 | 0–0.98 | 2.01 | 34% | 46% |
| 1.05 | 332 | **3.7** | 0.32 | 0–1.05 | 1.74 | **39%** | 32% |
| 1.13 | 742 | 4.2 | 0.39 | 0–1.16 | 1.89 | 40% | 32% |

- Bản cổ điển cho phổ **liên tục, hai đỉnh**:
  - khoảng 40% là va chạm **sượt**, nhả gần 0;
  - phần còn lại là một **bướu rộng** quanh 0.5–1.0 ħω0 (5–10 MeV).
- Trung bình 3.7–4.3 MeV khớp với 4.0–4.6 MeV của bài báo.
- Bản lượng tử thì **lượng tử hóa** cùng dòng năng lượng đó thành các vạch rô-to 12–35 MeV. Con số 0.4 ħω0 cổ điển là trung bình của các bậc thang này, như đã kết luận ở [schrodinger-rang-buoc](schrodinger-rang-buoc-kiem-chung.md).

## Kết luận

1. **Bên trong**, sau khi thoát khỏi 0, động cơ mang +42.6 MeV vượt chân không: 81% ở chuyển động quay, gần như không có ở dao động tử.
2. **Bên ngoài**, E năng lượng mà nó gặp **gồm nhiều năng lượng khác nhau**:
   - bản lượng tử: phổ vạch rời rạc, 23 vạch nhả ra, trải 0–60 MeV, ba vạch 35/25/25 MeV mang 85% công suất, 97% là bước rô-to;
   - bản cổ điển: phổ liên tục hai đỉnh, với khoảng 40% va chạm gần 0 và một bướu rộng 5–10 MeV.
3. **"Hàm toàn cục không cố định"** chỉ đạt được theo nghĩa **khối lớn dần**: dãy |Ψ_N⟩⟩ lớn 5.24 nat/T0 (khoảng 190 lần mỗi T0) khi ba pha xung đột, và ngừng lớn khi một pha thắng. Từng |Ψ_N⟩⟩ riêng lẻ vẫn tĩnh, vì đó là định nghĩa của Wheeler–DeWitt.

**Câu tiếng Anh cho bản thảo:**

> Once the engine leaves the frozen vacuum it carries 4.26 ħω₀ above the ground state, 81% of it rotational. The environment it then meets is a line spectrum: Carnot emission consists of 23 lines between 0 and 6 ħω₀, three rotor-ladder lines (3.54, 2.47, 2.45 ħω₀) carry 85% of the emitted power, and 97% of that power comes from rotor (ℓ-changing) transitions; the classical model shows the same flux as a continuous, bimodal spectrum (≈40% grazing impacts near zero, a broad hump at 0.5–1.0 ħω₀, mean 0.37–0.43 ħω₀). Each history state |Ψ_N⟩⟩ remains static, but the sequence grows: its effective branch number e^{H_N} increases by 5.24 nat per T₀ under three-phase conflict and stops growing when one phase dominates.

## Giới hạn

- **"Bên ngoài" = môi trường đã có trong mô hình** (bể Carnot + nguồn bơm). Mô hình không chứa một "vũ trụ khác" độc lập. Phổ ở đây là phần môi trường mà động cơ **trao đổi được**, không phải toàn bộ năng lượng của môi trường.
- Phổ vạch của kênh bơm không đáng tin trong gần đúng thế tục (thiếu khoảng 75% công suất). Cần một phân tích có tính các số hạng kết hợp.
- "Khối lớn dần" đưa vào chỉ số N như một siêu thời gian. Nó không làm một |Ψ⟩⟩ đơn lẻ thôi tĩnh.
- Ở f = 0.95, phổ cổ điển chỉ có 44 va chạm.

## Câu hỏi còn mở

- Muốn mô tả một "E năng lượng ngoài" có cấu trúc riêng (không chỉ là bể Markov phẳng), có thể cho môi trường một mật độ phổ J(ω) có đỉnh. Khi đó có thể xem các vạch rô-to 25–35 MeV có cộng hưởng với nó không, và phổ nhả ra thay đổi ra sao. Có làm tiếp không?
- Tính phổ kênh bơm chính xác bằng hàm tương quan hai thời điểm (định lý hồi quy lượng tử) thay cho gần đúng thế tục.
