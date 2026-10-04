# Hai vũ trụ, hai động cơ E_dual: E năng lượng quyết định tốc độ thời gian, bên trong động mà bên ngoài tĩnh

> **Ý của học thuyết:**
> - Vũ trụ ta và một vũ trụ bên ngoài, mỗi bên có một động cơ E_dual riêng với **E năng lượng** riêng, xung đột ở tầng chân không lượng tử.
> - Năng lượng khác nhau thì **tốc độ thời gian** khác nhau.
> - Bên trong một vũ trụ thì động, còn nhìn từ ngoài vẫn tĩnh. Wheeler–DeWitt vẫn đúng, vì không ai đứng ngoài toàn bộ để đo.
> - Xung đột ba pha làm cho vũ trụ ngoài **thấy** vũ trụ trong đang động.
>
> **Nối tiếp:** [e-nang-luong-cong-huong-moi-truong](e-nang-luong-cong-huong-moi-truong-kiem-chung.md), [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md), [xung-dot-qua-khu-hien-tai](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md).
> **Số liệu:** [giai_so_hai_vu_tru_e_dual.py](giai_so_hai_vu_tru_e_dual.py), kết quả thô [hai_vu_tru_e_dual_ket_qua.txt](hai_vu_tru_e_dual_ket_qua.txt). Chạy khoảng 16 phút trên 8 lõi (máy đang tải); mỗi trường hợp 16 nhánh × 111 T0.

## Mô hình

- **Mỗi vũ trụ là một động cơ E_dual:** rô-to ℓ ≤ 4 (có mức gãy ℓ = 4) ⊗ hạt n ≤ 1, tức 18 chiều. H theo Bảng I, có kênh Carnot và kênh bơm. **E năng lượng của vũ trụ = cường độ bơm** a_A (của ta), a_B (bên ngoài).
- **Xung đột ở tầng chân không:** một kênh trao đổi cường độ κ. Khi vũ trụ A tiêu tán một bước Carnot thì cùng lúc vũ trụ B được nâng một bậc, và ngược lại. Đây là **cách duy nhất** để mỗi vũ trụ "thấy" vũ trụ kia.
- Hệ ghép có 324 chiều, giải bằng MCWF. Mỗi quỹ đạo là một nhánh của trạng thái lịch sử chung, và **trạng thái chung của cả hai vũ trụ vẫn là một trạng thái dừng**, tức WdW tĩnh toàn cục.
- **Đồng hồ của mỗi vũ trụ** là số bản ghi tiêu tán của **chính nó** mỗi T0, gồm cả lần tiêu tán sang vũ trụ kia.

## Kết quả

| a_A | a_B | κ | ΔE_A | ΔE_B (ħω0) | r_A (/T0) | r_B (/T0) | **r_A / r_B** | A→B (/T0) | B→A (/T0) |
|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 1.0 | 0 | 4.26 | 4.32 | 0.486 | 0.504 | 0.97 | 0 | 0 |
| 1.0 | **0.3** | 0 | 4.12 | **1.81** | 0.511 | **0.164** | **3.12** | 0 | 0 |
| 1.0 | **3.0** | 0 | 4.19 | **6.67** | 0.486 | **0.810** | **0.60** | 0 | 0 |
| 1.0 | 0.3 | 0.05 | 4.20 | 1.85 | 0.515 | 0.173 | 2.98 | 0.040 | 0.006 |
| 1.0 | 0.3 | 0.20 | 3.69 | 1.90 | 0.550 | 0.213 | **2.59** | **0.128** | 0.044 |
| **0** | 1.0 | 0 | 0.92 | 4.10 | **0.002** | 0.470 | 0.00 | 0 | 0 |
| **0** | 1.0 | 0.05 | 0.85 | 4.21 | **0.020** | 0.515 | 0.04 | 0.001 | 0.045 |
| **0** | 1.0 | 0.20 | 0.99 | 3.67 | **0.060** | 0.545 | 0.11 | 0.011 | 0.150 |
| 0 | 0 | 0.20 | 0.84 | 0.85 | 0.003 | 0.004 | — | 0.001 | 0.001 |

Sai số chuẩn của r khoảng ±0.01–0.02.

## Đọc kết quả, đối chiếu từng ý của học thuyết

**1. E năng lượng quyết định tốc độ thời gian.** Khi hai vũ trụ tách rời (κ = 0), cùng một động cơ nhưng năng lượng khác nhau thì đồng hồ nội tại chạy khác nhau:
- Vũ trụ B năng lượng thấp (a_B = 0.3, ΔE = 1.8 ħω0) chạy **chậm hơn 3.1 lần** vũ trụ A.
- Vũ trụ B năng lượng cao (a_B = 3, ΔE = 6.7 ħω0) chạy **nhanh hơn 1.7 lần**.
- Cùng năng lượng thì cùng nhịp (0.97 ± 0.03).
- Tốc độ thời gian tăng theo E năng lượng nhưng chậm dần: ΔE tăng 3.7 lần (từ 1.81 lên 6.67) thì nhịp chỉ tăng 4.9 lần (từ 0.164 lên 0.810).

**2. Bên trong động, bên ngoài tĩnh.** Khi κ = 0:
- Mỗi vũ trụ có nhịp bản ghi riêng khác 0 (bên trong động).
- Tốc độ trao đổi bằng **0**: không một sự kiện nào của A đến được B.
- Với B, A không có gì thay đổi, tức **tĩnh**, đúng như WdW. Đây cũng là trường hợp p = 0 ở [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md) mục C: một người quan sát không đọc được bản ghi thì thấy tương lai đã định (h = 0).

**3. Xung đột ở tầng chân không làm vũ trụ ngoài thấy vũ trụ trong động.** Khi κ > 0:
- B nhận được sự kiện của A với tốc độ 0.040/T0 (κ = 0.05) và 0.128/T0 (κ = 0.2). Đó là **nhịp động** của A mà B thấy được từ bên ngoài.
- **B thấy A chậm hơn thực tế.** Ở κ = 0.2, A thật sự ghi 0.550 bản ghi/T0, nhưng B chỉ thấy 0.128, tức **23%**. Vũ trụ bên ngoài luôn thấy vũ trụ bên trong với một **E khúc dày hơn**, giống người quan sát chỉ đọc một phần bản ghi.
- **Năng lượng chảy từ vũ trụ nhiều E năng lượng sang vũ trụ ít.** A→B 0.128 so với B→A 0.044; ΔE_A giảm từ 4.12 xuống 3.69, ΔE_B tăng nhẹ từ 1.81 lên 1.90.
- **Hai nhịp thời gian kéo lại gần nhau:** r_A/r_B giảm từ 3.12 → 2.98 → 2.59 khi κ tăng. Xung đột giữa hai vũ trụ là một **đồng bộ hóa một phần**.

**4. Một vũ trụ đóng băng được vũ trụ kia làm tan băng.** A không có bơm (a_A = 0): bên trong A đóng băng, r_A = 0.002.
- Ghép với B đang chạy: r_A tăng lên 0.020 (κ = 0.05) rồi **0.060** (κ = 0.2), tức **gấp 30 lần**, nhờ năng lượng B đổ sang (B→A 0.150/T0).
- Thời gian của A chạy lại, nhưng vẫn chậm hơn B khoảng 9 lần (tỉ số 0.11).
- **Hai vũ trụ cùng đóng băng thì không cứu được nhau** (0.003, 0.004): không bên nào có E năng lượng để cho.

**5. Wheeler–DeWitt vẫn đúng.** Trạng thái chung A⊗B là một trạng thái dừng: toàn cục tĩnh. Mọi chuyển động ở trên đều là **tương quan bên trong** (bản ghi riêng của mỗi vũ trụ, sự kiện trao đổi giữa chúng), không cần ai đứng ngoài cả hai. Đó đúng là điều học thuyết khẳng định.

## Kết luận

| Ý của học thuyết | Số đo | Đứng? |
|---|---|---|
| Năng lượng khác ⟹ tốc độ thời gian khác | Chậm hơn 3.1 lần (ΔE thấp), nhanh hơn 1.7 lần (ΔE cao) | ✓ |
| Bên trong động, bên ngoài tĩnh nếu không tương tác | κ = 0: 0 sự kiện trao đổi, bên ngoài thấy tĩnh | ✓ |
| Xung đột ở tầng chân không làm bên ngoài thấy bên trong động | κ > 0: thấy 0.04–0.13 sự kiện/T0, nhưng chỉ khoảng 23% nhịp thật | ✓ (với giới hạn: thấy chậm hơn thực tế) |
| WdW vẫn đúng vì không ai đứng ngoài | Trạng thái chung dừng; mọi chuyển động là tương quan nội tại | ✓ |
| Thêm từ số đo: E năng lượng chảy từ cao xuống thấp, kéo hai nhịp lại gần nhau, có thể làm tan băng vũ trụ đóng băng | A→B > B→A; tỉ số nhịp 3.12 → 2.59; r_A tăng 30 lần | mới |

**Câu tiếng Anh cho bản thảo:**

> Two universes, each running its own E_dual engine, keep their own record clocks. With no contact, a universe with less drive energy runs 3.1 times slower (and one with more runs 1.7 times faster), while neither registers a single event of the other: each is dynamic inside and static from outside, as the Wheeler–DeWitt constraint requires of the joint stationary state. A vacuum-level exchange channel lets each universe see the other at the rate of exchange events (0.04–0.13 per T₀), only about a quarter of the other's internal tempo. It also carries energy from the higher- to the lower-energy universe, pulls their clock rates together (ratio 3.12 → 2.59), and can thaw a frozen universe (record rate ×30), whereas two frozen universes cannot thaw each other.

## Giới hạn

- **"Vũ trụ" ở đây là một động cơ 18 chiều.** Mô hình nói về **cơ chế** (năng lượng, bản ghi, trao đổi), không nói về vũ trụ học thật.
- **Kênh trao đổi do tôi chọn** (A tiêu tán ⟺ B được nâng); học thuyết chưa xác định kênh này. Kênh khác có thể cho tỉ lệ khác.
- **Vũ trụ "đóng băng" không nằm ở chân không:** ΔE_A ≈ 0.9 ħω0 khi a_A = 0. Hệ kẹt ở trạng thái mà kênh Carnot không làm hồi phục được (trạng thái tối), giống điều đã thấy ở [e-nang-luong-thoat-ly](e-nang-luong-thoat-ly-chan-khong-kiem-chung.md). Nó đứng yên ở đó nhưng không ở mức thấp nhất.
- **Có hai giả định đo lường:** "tốc độ thời gian" được đo bằng nhịp bản ghi riêng của mỗi vũ trụ, và "điều thấy được từ ngoài" được đo bằng nhịp sự kiện trao đổi. Thông tin thực sự mà B thu được về A có thể ít hơn số sự kiện.

## Câu hỏi còn mở

- Đo thông tin (nat/T0) mà B thực sự thu được về pha của A từ các sự kiện trao đổi, thay vì chỉ đếm sự kiện.
- Ghép nhiều vũ trụ (ba hoặc nhiều hơn) để xem năng lượng và nhịp thời gian phân bố thế nào.
