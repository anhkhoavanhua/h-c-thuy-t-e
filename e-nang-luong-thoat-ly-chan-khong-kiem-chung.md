# E năng lượng: năng lượng để thoát khỏi "0" (chân không, thời gian đóng băng)

> **Câu hỏi:** cần bao nhiêu năng lượng để động cơ $\mathcal E_{\text{dual}}$ vượt khỏi môi trường chân không lượng tử, tức thoát khỏi "0": trạng thái cơ bản mà ở đó không có bản ghi và thời gian đóng băng (h = 0)?
> **Nối tiếp lời giải frozen time:** [xung-dot-qua-khu-hien-tai](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md), [xung-dot-hien-tai-tuong-lai](xung-dot-hien-tai-tuong-lai-kiem-chung.md), [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md).
> **Số liệu:** [giai_so_e_nang_luong.py](giai_so_e_nang_luong.py), kết quả thô [e_nang_luong_ket_qua.txt](e_nang_luong_ket_qua.txt), chạy khoảng 6 giây. Phần B dùng số lần gãy theo nhiệt độ đã đo ở [xung_dot_mo_bang_ket_qua.txt](xung_dot_mo_bang_ket_qua.txt).

## A. Lượng tử (Part III): năng lượng vượt chân không và giá của thời gian

Carnot cố định b = 1; bơm a được quét. Các đại lượng tính trên trạng thái dừng chính xác:
- năng lượng vượt chân không ΔE = ⟨H⟩ − E₀;
- công suất bơm vào P;
- dòng bản ghi Carnot J_C;
- độ mở tương lai h (tính trên nhánh như ở [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md)).

Chân không của S có E₀ = −0.036 ħω0 (đã gồm năng lượng điểm không). Khe tới mức kích thích đầu là 0.866 ħω0.

| a (bơm) | ΔE (ħω0) | P vào (ħω0/T0) | J_C (/T0) | **h (nat/T0)** | **P/h (ħω0 mỗi nat)** | P/J_C (ħω0 mỗi bản ghi) |
|---|---|---|---|---|---|---|
| **0 (chân không)** | **0** | 0 | 0 | **0.000** | — | — |
| 0.01 | 0.572 | 0.006 | 0.005 | 0.166 | **0.037** | 1.23 |
| 0.03 | 0.713 | 0.019 | 0.013 | 0.369 | 0.052 | 1.43 |
| 0.1 | 0.962 | 0.075 | 0.044 | 1.013 | 0.074 | 1.71 |
| 0.3 | 1.731 | 0.322 | 0.150 | 2.493 | 0.129 | 2.14 |
| **1 (E_dual)** | **4.262** | **1.370** | **0.494** | **5.234** | **0.262** | **2.77** |
| 3 | 6.658 | 2.606 | 0.824 | 6.786 | 0.384 | 3.16 |
| 10 | 7.922 | 3.340 | 0.997 | 7.231 | 0.462 | 3.35 |

1. **Năng lượng chân không không làm thời gian chạy.** Ở a = 0, hệ ở trạng thái cơ bản, có đủ năng lượng điểm không, nhưng h = 0.000 và J = 0. Năng lượng điểm không là năng lượng **bị khóa**: không có dòng qua đâu cả nên không sinh bản ghi.
2. **Không có ngưỡng sắc.** Bất kỳ dòng năng lượng nào đi qua hệ (a > 0) cũng làm h > 0. Ngay với a = 0.01, công suất 0.006 ħω0/T0 đã cho h = 0.17 nat/T0. Đại lượng quyết định việc "thoát khỏi 0" **không phải năng lượng chứa trong hệ** mà là **dòng năng lượng đi qua** nó (bơm vào, Carnot ra). Điều này khớp với lời giải frozen time: cần cả hai kênh.
3. **E năng lượng theo nghĩa giá:** mỗi nat tương lai mở ra tốn
   - 0.037 ħω0 (0.37 MeV) khi bơm yếu,
   - **0.26 ħω0 (2.6 MeV)** ở E_dual,
   - 0.46 ħω0 khi bơm rất mạnh.

   Thời gian chạy nhanh hơn thì đắt hơn trên mỗi nat: P/h tăng 12 lần trong khi h tăng 44 lần.
4. **Mỗi bản ghi Carnot** tốn 1.2–3.4 ħω0, và **2.8 ħω0 (28 MeV)** ở E_dual. Con số này cùng bậc với bước nhảy rô-to khoảng 2 ħω0 mỗi bản ghi ([part3-latex-so-lieu](part3-latex-so-lieu-kiem-chung.md)).
5. ΔE không giảm về 0 tuyến tính khi a → 0 (a = 0.01 vẫn cho ΔE = 0.57). Có thể một phần dân số nằm ở các trạng thái mà kênh Carnot không làm hồi phục được (trạng thái tối). Chưa kiểm.

**So với Landauer:** với k_BT = 1 MeV của học thuyết, giới hạn là k_BT = 0.1 ħω0 mỗi nat. E_dual tốn 0.26 ħω0/nat, tức **2.6 lần** giới hạn. Lưu ý: trong mô hình Part III này, bể Carnot ở T = 0, nên giới hạn Landauer đúng nghĩa bằng 0. Con số 0.037 < 0.1 ở bơm yếu **không vi phạm** gì; phép so với k_BT = 1 MeV chỉ để tham khảo.

## B. Cổ điển (Part I): thoát khỏi điểm hút Quá khứ

| k_BT (MeV) | 0 | 0.2 | 0.5 | 1.0 |
|---|---|---|---|---|
| Lần gãy / T0 | **0** | 0.188 | 0.525 | 1.025 |
| Gãy / T0 chia cho k_BT (MeV) | — | 0.94 | 1.05 | 1.03 |

- **Ở k_BT = 0, hệ không thoát được.** Điểm hút tất định nằm dưới ngưỡng 13.7%. **Khe tất định** là $\tfrac12\langle M\rangle\dot\theta_{\text{br}}^2(1-M_{\text{det}})=$ **0.74 ħω0 = 7.4 MeV**: đó là năng lượng phải thêm vào rô-to nếu không có nhiễu.
- **Khớp Arrhenius** cho Δ = **0.41 MeV**, trùng với Δ_eff ≈ 0.4 MeV của bài báo. Hàng rào hiệu dụng nhỏ hơn khe tất định 18 lần, vì nhiễu chỉ cần đẩy pha θ tới góc có |λ| lớn nhất, không cần nâng cả năng lượng quay.
- **Lưu ý trung thực:** ba điểm dữ liệu cũng khớp **tỉ lệ thuận với k_BT** gần như hoàn hảo (gãy/k_BT = 0.94, 1.05, 1.03). Khi k_BT ≥ Δ, hàm mũ không còn là yếu tố chi phối, nên con số Δ = 0.41 MeV **không phân biệt được** với quy luật tuyến tính. Bài báo cũng nói tốc độ thoát "bão hòa" khi vượt Δ_eff. Muốn chốt Δ cần đo ở k_BT ≪ 0.4 MeV (0.05–0.2 MeV), với quỹ đạo dài hơn vì sự kiện hiếm.

## Thang E năng lượng của học thuyết

| Đại lượng | Giá trị | Nghĩa |
|---|---|---|
| Năng lượng điểm không (chân không) | có, nhưng h = 0 | Không làm thời gian chạy |
| Hàng rào kích hoạt hiệu dụng (cổ điển) | ≈ 0.4 MeV (chưa chốt, xem B) | Nhiễu nhiệt cần cỡ này để Quá khứ bắt đầu gãy |
| **Giá của thời gian** (lượng tử, E_dual) | **0.26 ħω0 = 2.6 MeV mỗi nat** | Năng lượng tiêu tán cho mỗi nat tương lai mở ra |
| Khe tất định (cổ điển) | 0.74 ħω0 = 7.4 MeV | Năng lượng quay cần thêm nếu không có nhiễu |
| Năng lượng mỗi bản ghi Carnot (lượng tử) | 2.8 ħω0 = 28 MeV | Dòng vào trên mỗi bản ghi |
| Ngưỡng gãy E_br | 4.706 ħω0 = 47 MeV | Mức năng lượng của pha Hiện tại |

**Câu trả lời ngắn:**
- "Thoát khỏi 0" **không đòi một lượng năng lượng chứa cố định**. Chân không có năng lượng điểm không mà vẫn đóng băng.
- Nó đòi một **dòng năng lượng đi qua hệ**: có vào (Quá khứ bị đẩy lên ngưỡng) và có ra (Carnot ghi bản ghi). Giá ở E_dual là khoảng **2.6 MeV cho mỗi nat** tương lai, gấp 2.6 lần Landauer ở 1 MeV.
- Về phía cổ điển, ở k_BT = 0 không thoát được; khi có nhiễu, thang kích hoạt cỡ 0.4 MeV nhưng con số này chưa được chốt.

**Câu tiếng Anh cho bản thảo:**

> Escaping the frozen vacuum does not require a stored energy threshold: the ground state carries its zero-point energy yet produces no records (h = 0). What is required is an energy current through the engine, drive in and Carnot dissipation out. Any non-zero drive yields h > 0, at a cost of 0.04–0.46 ħω₀ per nat of opened future (0.26 ħω₀ ≈ 2.6 MeV per nat in the operating regime, 2.6 times the Landauer bound at k_BT = 1 MeV) and 2.8 ħω₀ per Carnot record. Classically, the deterministic attractor lies 0.74 ħω₀ below breakdown; thermally activated escape gives an Arrhenius barrier of 0.41 MeV, consistent with the effective 0.4 MeV, although at k_BT ≥ 0.2 MeV the breakdown rate is equally consistent with a linear dependence on k_BT.

## Câu hỏi còn mở

- Chạy lại phần cổ điển ở k_BT = 0.05–0.2 MeV với quỹ đạo dài (vài trăm T0) để tách Arrhenius khỏi quy luật tuyến tính và chốt Δ.
- Vì sao ΔE không về 0 khi bơm rất yếu? Kiểm phổ của Liouvillian khi chỉ có Carnot để tìm trạng thái tối.
