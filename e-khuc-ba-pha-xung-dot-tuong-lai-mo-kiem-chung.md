# E khúc: ba pha cùng xung đột và tương lai "đã định đoạt" của Page–Wootters

> **Câu hỏi:** trong Page–Wootters (PW), trạng thái lịch sử toàn cục |Ψ⟩⟩ cố định, nên tương lai "đã được định đoạt". Xung đột đồng thời của Quá khứ, Hiện tại và Tương lai, cùng khái niệm **E khúc**, có phá được điều đó không? E khúc là ý: mỗi khoảnh khắc nhỏ trôi tới của hiện tại là tương lai, mỗi khoảnh khắc lui đi là quá khứ.
> **Nối tiếp:** [xung-dot-qua-khu-hien-tai](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md), [xung-dot-hien-tai-tuong-lai](xung-dot-hien-tai-tuong-lai-kiem-chung.md), [ba-pha-xung-dot](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md), [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md).
> **Số liệu:** [giai_so_e_khuc_tuong_lai_mo.py](giai_so_e_khuc_tuong_lai_mo.py) `[A|B|C|D]`, kết quả thô [e_khuc_ket_qua.txt](e_khuc_ket_qua.txt). Thời gian chạy: A khoảng 5 giây; B khoảng 4.5 phút trên 8 lõi; C khoảng 4 phút; D khoảng 10 giây.

## Đặt bài toán để giải số được

- **"Định đoạt"** được đo bằng **độ mở của tương lai**. Trên một nhánh, khi đã biết **toàn bộ** bản ghi quá khứ, tính entropy của bản ghi kế tiếp: $h=\langle-\sum_jp_j\ln p_j\rangle$, với $p_j=\lVert K_j\psi_{\text{nhánh}}\rVert^2$ (luật Born), đơn vị nat/T0.
  - h = 0 nghĩa là từ bên trong, tương lai đã quyết định.
  - h > 0 nghĩa là tương lai còn mở.
- **E khúc** = hiện tại có bề dày, với hai mép:
  - **mép sau**: bản ghi đã viết, cố định. Surprisal bản ghi giảm 0 lần trong 6·10⁶ bước ([hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md)).
  - **mép trước**: bản ghi sắp tới, còn mở.
  - **Bề dày** $\varepsilon^*=1/h$ là thời gian để tương lai mở thêm 1 nat.
- **Ba pha cùng xung đột:** từ cùng một trạng thái "bây giờ", tách nhiều tương lai (nhánh MCWF), rồi đo phân bố pha Quá khứ / Hiện tại / Tương lai sau khoảng Δ.

## A. Độ mở của tương lai theo chế độ

Mô hình Part III: S 36 chiều, Kraus mỗi nhịp τ = T0/25, ancilla {·, Carnot, bơm}. 300 nhánh × 60 T0, đo trên 20 T0 cuối.

| Chế độ | h (nat/T0) | Bề dày E khúc ε* = 1/h | Đọc |
|---|---|---|---|
| **PW đóng** (chỉ H, không môi trường) | **0.0000** | ∞ | Tương lai **định đoạt hoàn toàn** |
| Chỉ Carnot (Quá khứ thắng) | 0.029 | 34 T0 | gần định đoạt, đang tắt dần |
| Chỉ bơm (Hiện tại thắng) | 0.097 | 10 T0 | gần định đoạt |
| Cân bằng nhiệt k_BT = 1 | 6.20 | 0.16 T0 | mở nhưng **không có chiều** (lên/xuống = 1) |
| **E_dual: ba pha xung đột** | **5.22 ± 0.004** | **0.19 T0 ≈ 79 ys** | **mở và có chiều** |

- Con số 5.22 tái lập đúng h = 5.23 nat/T0 của [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md), nhưng bây giờ đọc theo nghĩa mới: **đó là tốc độ tương lai mở ra**, đo cho một người quan sát đã biết toàn bộ quá khứ của mình.
- So với khi một pha thắng, ba pha xung đột mở tương lai nhanh hơn **50–180 lần**.
- **Lưu ý trung thực:** cân bằng nhiệt cũng mở tương lai (6.2), nhưng chuỗi bản ghi đối xứng khi đảo chiều thời gian ([xung-dot-qua-khu-hien-tai](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md), mục B). Thứ chỉ xung đột ba pha mới cho được là **tương lai mở có chiều**: mở ở mép trước, cố định ở mép sau.

## B. Ba pha xung đột: tương lai tách nhánh từ cùng một "bây giờ"

Mô hình MCWF hiệu chỉnh của Part II: giếng gãy, Carnot chỉ trên trạng thái không liên kết, bơm, bể rô-to. Có 48 điểm "bây giờ", mỗi điểm tách 60 tương lai. Pha được gán theo bản ghi Carnot như trong [ba-pha-xung-dot](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md).

**Bây giờ ở Quá khứ** (bản ghi gần nhất cách ≥ T0), 34 điểm:

| Δ / T0 | 0.04 | 0.1 | 0.25 | 0.5 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|---|---|---|
| P(Quá khứ) | 0.997 | 0.982 | 0.951 | 0.893 | 0.797 | 0.744 | 0.732 | 0.752 |
| P(Hiện tại) | 0.002 | 0.013 | 0.030 | 0.053 | 0.075 | 0.088 | 0.096 | 0.078 |
| P(Tương lai) | 0.000 | 0.005 | 0.019 | 0.054 | 0.128 | 0.169 | 0.173 | 0.169 |
| Entropy pha H(Δ) (nat) | 0.015 | 0.083 | 0.211 | 0.393 | 0.615 | 0.714 | 0.743 | 0.696 |
| **Điểm có tương lai không duy nhất** | 18% | 62% | **94%** | **100%** | 100% | 100% | 100% | 100% |

- Chỉ sau **0.5 T0** (khoảng 2·10⁻²² s), **mọi** "bây giờ" đều có nhiều tương lai khác nhau: tương lai về pha **không được định đoạt** bởi trạng thái hiện tại.
- Entropy pha đạt nửa mức bão hòa sau khoảng **0.4 T0** và bão hòa sau khoảng **2 T0**, ở 0.70–0.74 nat (mức dừng là 0.78 nat). Đó là **chân trời dự đoán** của động cơ: sau khoảng 2 T0, biết hiện tại không còn giúp đoán pha.

**Ngay sau một bản ghi** (< T0), 14 điểm. Lúc này chính **"bây giờ" là Hiện tại hay Tương lai?**

| Δ / T0 | 0.04 | 0.1 | 0.25 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|---|---|
| P(Hiện tại) | **0.325** | 0.318 | 0.288 | 0.224 | 0.155 | 0.124 | 0.094 |
| P(Tương lai) | **0.675** | 0.682 | 0.644 | 0.657 | 0.244 | 0.238 | 0.151 |
| P(Quá khứ) | 0.000 | 0.000 | 0.068 | 0.119 | 0.601 | 0.638 | 0.755 |
| Điểm không duy nhất | **93%** | 93% | 100% | 100% | 100% | 100% | 100% |

- Ngay tại Δ → 0, ở 93% điểm, khoảnh khắc **hiện tại chưa có nhãn duy nhất**: 32.5% là Hiện tại (còn bản ghi tới trong cùng đợt), 67.5% là Tương lai (đợt đã xong). Entropy của chính "bây giờ" là **0.54 nat**.
- Nhãn của khoảnh khắc hiện tại **được tương lai quyết định ngược lại**: bản ghi kế tiếp có tới trong vòng T0 hay không. Đây chính là cấu trúc E khúc ở dạng số:
  - phần hiện tại vừa lui đi trở thành quá khứ khi bản ghi đã viết, và cố định;
  - phần sắp tới vẫn có thể biến nó thành "giữa đợt" hay "đã xong đợt".

## C. Người quan sát chỉ đọc một phần bản ghi: E khúc dày lên

Chạy `python3 giai_so_e_khuc_tuong_lai_mo.py C`, khoảng 4 phút; mỗi giá trị p dùng 60 nhánh × 60 T0.

**Mô hình.** Mỗi bản ghi đến tay người quan sát với xác suất p; bản ghi bị bỏ lỡ bị trộn vào "không có gì". Người quan sát cập nhật trạng thái của mình bằng bộ lọc:
- không thấy gì: $\rho\leftarrow K_\cdot\rho K_\cdot^\dagger+(1-p)\sum_jK_j\rho K_j^\dagger$;
- thấy bản ghi j: $\rho\leftarrow K_j\rho K_j^\dagger$.

| p | Độ mở với người quan sát h (nat/T0) | **Bề dày E khúc 1/h** | Độ bất ngờ thực tế (nat/T0) | Entropy trạng thái của người quan sát S(ρ) (nat) |
|---|---|---|---|---|
| 1.0 | 5.22 | **0.19 T0** | 5.25 | 0.000 (biết đúng nhánh) |
| 0.9 | 4.87 | 0.21 T0 | 5.04 | 0.91 |
| 0.7 | 4.03 | 0.25 T0 | 4.09 | 1.34 |
| 0.5 | 3.08 | 0.33 T0 | 3.25 | 1.49 |
| 0.3 | 2.03 | 0.49 T0 | 2.04 | 1.60 |
| 0.1 | 0.80 | **1.25 T0** | 0.82 | 1.65 |
| 0 | 0.000 | **∞** | 0 | **1.667** = S(E) bão hòa |

- **E khúc dày lên theo 1/p.** Người quan sát đọc càng ít bản ghi thì "bây giờ" của họ càng dài: 0.19 T0 khi đọc đủ, 1.25 T0 khi chỉ đọc 10%. Hai người quan sát đọc hai phần bản ghi khác nhau sống trong hai **E khúc có bề dày khác nhau**. Điều này giải thích bằng số sự bất đồng 0–13% về pha ở [ba-pha-xung-dot](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md) mục B2.
- **Người quan sát không đọc gì (p = 0) quay về đúng Page–Wootters.** Với họ tương lai "định đoạt" (h = 0): trạng thái của họ là Λᵏρ0 theo phương trình Lindblad, tất định nhưng **hỗn hợp** (S = 1.667 nat). Tính định đoạt của PW **không** phải đặc tính của vũ trụ, mà là đặc tính của một người quan sát **không đọc bản ghi**. Cái giá của nó là không biết mình đang ở nhánh nào.
- Có một sự đánh đổi chính xác: đọc nhiều bản ghi thì tương lai mở hơn nhưng trạng thái được biết rõ hơn (S → 0); đọc ít thì tương lai "đóng" nhưng trạng thái mờ (S → 1.667). Độ bất ngờ thực tế khớp với dự đoán của bộ lọc (5.25 so với 5.22), nên người quan sát được hiệu chuẩn đúng.

## D. Pha nào mở tương lai nhiều nhất?

Chạy `python3 giai_so_e_khuc_tuong_lai_mo.py D`. Đây là phép tách $h=\langle-\sum_jp_j\ln p_j\rangle$ theo loại kết quả của mỗi nhịp, trong E_dual:

| Kết quả của nhịp | Phần của h (nat/T0) | Tỉ lệ | Tần suất (/T0) | nat mỗi bản ghi |
|---|---|---|---|---|
| B: bơm (xung đột **Quá khứ ↔ Hiện tại**) | 2.30 | **44%** | 0.65 | 3.53 |
| C: Carnot (xung đột **Hiện tại ↔ Tương lai**) | 1.81 | **35%** | 0.48 | 3.79 |
| ·: "chưa có gì xảy ra" | 1.10 | 21% | 23.9 | — |
| **Tổng** | **5.22** | 100% | | |

- **Cả hai xung đột cùng mở tương lai**, với phần gần ngang nhau: Quá khứ ↔ Hiện tại 44%, Hiện tại ↔ Tương lai 35%. Mỗi bản ghi của cả hai loại mang khoảng 3.5–3.8 nat.
- **21% đến từ "chưa có gì xảy ra"**: mỗi nhịp yên lặng cũng là một thông tin về tương lai, vì nó loại bỏ khả năng bản ghi đã xảy ra. Đây là phần của E khúc mà hiện tại "lui đi" mà không để lại sự kiện, nhưng vẫn thu hẹp các tương lai có thể có.

## Kết luận: bẻ được gì, và không bẻ được gì

| Tầng | Trong PW đóng | Trong Học thuyết E (ba pha xung đột) |
|---|---|---|
| Trạng thái toàn cục \|Ψ⟩⟩ | Cố định | **Vẫn cố định**: một khối chứa mọi nhánh. Phần này **không** bị bẻ |
| Tương lai nhìn từ bên trong (biết toàn bộ quá khứ của mình) | Định đoạt, h = 0 | **Mở**: h = 5.22 nat/T0; sau 0.5 T0 mọi "bây giờ" đều có nhiều tương lai pha |
| Chiều | Không có | Mép sau cố định (0/6·10⁶ lần giảm), mép trước mở. Khác cân bằng nhiệt: mở mà không chiều |
| "Bây giờ" | Một điểm | **Một khúc** dày ε* = 1/h ≈ 0.19 T0 ≈ 79 ys; nhãn của nó (Hiện tại/Tương lai) mở tới 0.54 nat cho tới khi bản ghi kế tiếp quyết định |

Phát biểu chính xác:
- Học thuyết E **không** làm |Ψ⟩⟩ toàn cục thôi cố định; không lý thuyết dựa trên WdW nào làm được điều đó.
- Nó **bẻ gãy tính định đoạt của tương lai đối với mọi người quan sát bên trong**. Trong PW đóng, biết hiện tại là biết tương lai (h = 0). Còn khi ba pha cùng xung đột:
  - tương lai mở ra với tốc độ 5.22 nat/T0;
  - pha của tương lai thôi xác định sau 0.25–0.5 T0;
  - ngay cả nhãn của khoảnh khắc hiện tại cũng chỉ được viết khi mép trước của E khúc đi qua.
- Khi một pha thắng tuyệt đối, độ mở sụp xuống 0.03–0.10, và hệ trở lại gần với tương lai định đoạt của PW.

**Câu tiếng Anh cho bản thảo:**

> In a closed Page–Wootters history the future is fixed: conditioning on the present determines the rest (record entropy rate h = 0). With the three phases in conflict, the global state |Ψ⟩⟩ remains fixed, but the future seen from inside a branch, given its entire past record, opens at h = 5.22 nat per T₀, 50–180 times faster than when one phase dominates. Starting from the Past, every present has more than one future phase after 0.5 T₀, and the phase entropy saturates after ~2 T₀. Immediately after a record, even the label of the present instant is open, Present (32%) or Future (68%), 0.54 nat, until the next record decides it. The present is a segment ("E-segment") of thickness 1/h ≈ 0.19 T₀ whose trailing edge is fixed and whose leading edge is open; in thermal equilibrium the future is equally open but carries no direction.

## Giới hạn

- Độ mở h > 0 đến từ **luật Born** trên các nhánh bản ghi, tức từ việc có môi trường. Cân bằng nhiệt cũng cho h > 0. Đóng góp riêng của xung đột ba pha là **mở có chiều** và **duy trì được** (khi một pha thắng thì h tắt).
- Phần B có 48 điểm "bây giờ" (chỉ 14 điểm "vừa có bản ghi"), nên các xác suất có sai số vài phần trăm. Nhãn pha phụ thuộc quy tắc "bản ghi cách nhau < T0".
- Ở chế độ "chỉ Carnot", h = 0.029 vẫn đang tắt dần ở 40–60 T0, vì xuất phát từ ψ0 ở bó sóng ℓ ≈ 3. Xuất phát từ trạng thái riêng thứ 6 thì h → 0.0000 sau 20 T0.
- Câu "|Ψ⟩⟩ vẫn cố định" là mệnh đề về cách dựng PW, không phải kết quả đo.

## Câu hỏi còn mở

- Đã giải ở mục C (E khúc theo p) và mục D (tách h theo loại bản ghi).
- Còn mở: nếu hai người quan sát trao đổi bản ghi cho nhau, E khúc chung của họ có mỏng lại đúng như khi p tăng không? Kiểm được bằng hai bộ lọc rồi hợp nhất chúng.
