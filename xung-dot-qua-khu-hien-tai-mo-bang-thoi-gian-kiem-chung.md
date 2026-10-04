# Xung đột Quá khứ ↔ Hiện tại làm thời gian "tan băng": giải số

> **Luận điểm của học thuyết:** bài toán đóng băng thời gian của Wheeler–DeWitt được giải nhờ **xung đột giữa Quá khứ và Hiện tại**.
> **Cơ sở:** HocThuyetE.pdf … HocThuyetE (8).pdf (Part I: động cơ, gãy, Carnot; Part III: Lindblad lượng tử + bản ghi). Nối tiếp [ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md), [hai-vai-tro-kim-va-be-kiem-chung.md](hai-vai-tro-kim-va-be-kiem-chung.md), [wdw-bang-0-hoc-thuyet-e-tinh-ra-gi-kiem-chung.md](wdw-bang-0-hoc-thuyet-e-tinh-ra-gi-kiem-chung.md).
> **Số liệu:** [giai_so_xung_dot_mo_bang.py](giai_so_xung_dot_mo_bang.py), kết quả thô [xung_dot_mo_bang_ket_qua.txt](xung_dot_mo_bang_ket_qua.txt). Phần A, B chạy dưới 30 giây; phần C khoảng 30 giây trên 8 lõi.

## Phát biểu thành mệnh đề kiểm được

1. Trạng thái WdW là **dừng toàn cục**, và mọi trạng thái dừng đều "đóng băng" ở cấp toàn cục.
2. Thời gian **quan hệ** chỉ chạy khi hệ liên tục sinh **bản ghi mới, bất khả nghịch**: một dòng chu trình $J>0$ có **chiều**. Đây cũng là đại lượng đã cho mũi tên h = 5.23 nat/T0 ([hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md)).
3. **Luận điểm xung đột:**
   - Quá khứ, tức trạng thái khóa với lực đẩy lên ngưỡng gãy: kênh bơm/động cơ.
   - Hiện tại, tức gãy và tiêu tán Carnot kéo về: kênh Carnot.
   - Mệnh đề: J > 0 **chỉ khi cả hai cùng hoạt động**. Nếu một bên thắng tuyệt đối thì J = 0 và thời gian đóng băng.

## A. Lượng tử (Part III, S 36 chiều)

Bơm được nhân hệ số a, Carnot nhân hệ số b. Với mỗi (a, b), tìm trạng thái dừng chính xác (nhân của Liouvillian) và tính dòng bản ghi Carnot $J_C=\mathrm{Tr}(L_C\rho_{ss}L_C^\dagger)$, đơn vị bản ghi mỗi $T_0$. "Quá khứ" là các mức dưới ngưỡng gãy $E_{\text{br}}=4.706\,\hbar\omega_0$ (18/36 mức).

**Dòng bản ghi J_C (bản ghi / T0):**

| a \ b | 0 | 0.1 | 0.3 | 1 | 3 | 10 |
|---|---|---|---|---|---|---|
| **0** (không bơm) | — | **0** | **0** | **0** | **0** | **0** |
| 0.1 | **0** | 0.040 | 0.046 | 0.044 | 0.040 | 0.036 |
| 0.3 | **0** | 0.082 | 0.142 | 0.150 | 0.133 | 0.124 |
| 1 (E_dual) | **0** | 0.101 | 0.255 | **0.494** | 0.560 | 0.552 |
| 3 | **0** | 0.108 | 0.303 | 0.824 | 1.536 | 1.904 |
| 10 | **0** | 0.111 | 0.323 | 0.997 | 2.446 | 4.345 |

**Phần trạng thái ở Quá khứ, P(Quá khứ):**

| a \ b | 0.1 | 1 | 10 |
|---|---|---|---|
| 0 | **1.000** | **1.000** | **1.000** |
| 1 | 0.022 | 0.447 | 0.893 |
| 10 | 0.001 | 0.019 | 0.220 |
| chỉ bơm (b = 0) | **0.000** | | |

**Ba chế độ:**

| Chế độ | Trạng thái | J_C | Thời gian |
|---|---|---|---|
| **Quá khứ thắng** (không bơm) | ρ → cơ bản, P(QK) = 1 | 3.8·10⁻¹⁵ | **đóng băng** |
| **Hiện tại thắng** (không Carnot) | kẹt trên ngưỡng, P(QK) = 0; bơm cũng dừng (J_B = 0.011/T0) | 0 | **đóng băng** |
| **Xung đột** (cả hai) | Quá khứ và Hiện tại đồng hiện diện | 0.04–4.3 | **chạy** |

- **Hai bên xung đột cực đại.** Trạng thái đích của "chỉ bơm" và "chỉ Carnot" cách nhau ‖·‖_tr = **1.000**, tức hoàn toàn trực giao. Hệ bị kéo về hai đích loại trừ nhau, và dòng bản ghi chính là sự giằng co giữa chúng.
- **Bên yếu hơn quyết định nhịp.** Hàng a = 0.1 giữ J ≈ 0.04 với mọi b; cột b = 0.1 giữ J ≈ 0.1 với mọi a. Thời gian chạy nhanh hay chậm bị nút thắt bởi phía yếu hơn trong cuộc xung đột.
- **Tắt bơm thì thời gian đóng băng dần.** Từ trạng thái dừng a = b = 1, khi tắt bơm (giữ Carnot), J_C giảm: 0.494 (0) → 0.339 (1 T0) → 0.107 (5 T0) → **0.0048** (20 T0) → 0.0024 (100 T0).
- **Lưu ý trung thực:** J **không** đạt cực đại khi Quá khứ và Hiện tại cân bằng 50/50. Tương quan giữa J và P_QK·P_HT chỉ +0.34. Ví dụ a = b = 10 cho J lớn nhất (4.3) với P(QK) = 0.22. Điều kiện đúng là "**cả hai kênh khác 0**", và cường độ do bên yếu hơn quyết định, chứ không phải "cân bằng".

## B. Đối chứng: có tích tắc nhưng không có chiều (cân bằng nhiệt)

Không có bơm; bể Carnot ở nhiệt độ $k_BT$ có cân bằng chi tiết, tức có cả kênh lên. Bản ghi vẫn sinh ra liên tục:

| k_BT (ħω0) | J_xuống / T0 | J_lên / T0 | Dòng ròng / T0 | Lên / xuống |
|---|---|---|---|---|
| 0.1 | 0.0001 | 0.0001 | +3·10⁻⁸ | 0.9996 |
| 0.3 | 0.0599 | 0.0599 | +1·10⁻⁵ | 0.9998 |
| 1.0 | 0.8310 | 0.8310 | +2·10⁻⁵ | 1.0000 |
| 3.0 | 1.8666 | 1.8666 | +1·10⁻⁵ | 1.0000 |
| **E_dual** (bơm + Carnot T = 0) | **0.494** | **0** | **+0.494** | **0** |

Ở cân bằng, mỗi bản ghi "xuống" có một bản ghi "lên" tương ứng, nên chuỗi bản ghi chạy ngược vẫn hợp lệ như chạy xuôi. Có **tích tắc** nhưng **không có chiều**. Riêng việc có bản ghi chưa đủ để thời gian tan băng; phải có một **dòng ròng một chiều**. Trong học thuyết, dòng ròng đó là chu trình **Quá khứ → (gãy) → Hiện tại → (Carnot) → Tương lai → Quá khứ**, với ma trận chuyển 559 lần xuôi và 0 lần ngược ([ba-pha-xung-dot](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md)).

## C. Cổ điển (Part I, Langevin đầy đủ)

**Theo nhiệt độ** (f = 1.05; 2 quỹ đạo × 40 T0 mỗi giá trị):

| k_BT (MeV) | Lần gãy / T0 | Carnot tiêu tán / T0 (ħω0) |
|---|---|---|
| **0** | **0.000** | **0.000** |
| 0.2 | 0.19 | 0.10 |
| 0.5 | 0.53 | 0.24 |
| 1.0 | 1.03 | 0.38 |

**Theo tỉ số lái f** (k_BT = 1 MeV; 4 quỹ đạo × 120 T0, dữ liệu `rl_data`):

| f | 0.95 | 1.05 | 1.13 |
|---|---|---|---|
| Lần gãy / T0 | 0.094 | 0.698 | 1.548 |
| Carnot / T0 | 0.040 | 0.258 | 0.655 |

- **Ở k_BT = 0, Quá khứ thắng tuyệt đối.** Điểm hút tất định nằm dưới ngưỡng 13.7% ($M_{\text{det}}=0.863$), nên không gãy, không có bản ghi, cơ học cổ điển cũng "đóng băng" theo nghĩa này.
- Thăng giáng (nhiệt độ) và lực lái mạnh hơn (f) đẩy Quá khứ vào xung đột với ngưỡng Hiện tại, và nhịp bản ghi tăng theo.
- Ở f = 1.05, k_BT = 1 có hai ước lượng 1.03 và 0.70 lần gãy/T0. Chênh lệch là thăng giáng thống kê: quỹ đạo 40 T0 so với 120 T0.

## Kết luận

| Tầng | Kết quả số |
|---|---|
| Toàn cục (WdW) | Dừng: đóng băng ở cấp toàn cục là đúng, không phải nghịch lý |
| Chỉ Quá khứ | J = 0 (lượng tử: 4·10⁻¹⁵; cổ điển ở k_BT = 0: 0 lần gãy): **đóng băng** |
| Chỉ Hiện tại | J = 0, kẹt trên ngưỡng, bơm cũng dừng: **đóng băng** |
| Cân bằng nhiệt | Có tích tắc, dòng ròng 0, đối xứng đảo chiều: **đóng băng theo chiều** |
| **Xung đột Quá khứ ↔ Hiện tại** | **J > 0 một chiều** (0.04–4.3 bản ghi/T0), lên/xuống = 0: **thời gian chạy và có mũi tên** |

Luận điểm của học thuyết **đứng được** theo nghĩa chính xác sau:
- Trong một vũ trụ WdW tĩnh, dòng bản ghi một chiều (thời gian quan hệ có hướng) tồn tại **khi và chỉ khi** cả lực đẩy về phía gãy (Quá khứ) lẫn tiêu tán Carnot (Hiện tại) cùng khác 0.
- Hai lực này kéo hệ về hai trạng thái đích **trực giao** (khoảng cách vết bằng 1). Xung đột của chúng sinh ra chu trình ba pha, và chu trình ba pha là đồng hồ đồng thời là mũi tên.

**Câu tiếng Anh cho bản thảo:**

> The Wheeler–DeWitt state is globally stationary, so "frozen time" is correct at the global level. Relational time requires a non-vanishing, one-way record current. In the engine this current exists only when the drive toward breakdown (Past) and Carnot dissipation (Present) are both active: they pull the system toward orthogonal attractors (trace distance 1.000), and their conflict sustains the cycle Past → Present → Future at J_C = 0.49 records per T₀. If either side wins, J_C = 0 (ground state for pure dissipation, 4·10⁻¹⁵; stuck above threshold for pure drive), and in a detailed-balance bath records are produced but carry no direction (up/down ratio 1.0000). Time unfreezes through the conflict of Past and Present.

## Giới hạn

- "Trạng thái WdW" được thay bằng **trạng thái dừng** của Lindblad, tương đương phần hệ của trạng thái lịch sử khi đã ở chế độ dừng. Chưa dựng một |Ψ⟩⟩ có cả rô-to, hạt và thanh ghi bản ghi.
- J không cực đại ở cân bằng 50/50 (tương quan +0.34). "Xung đột" phải hiểu là **cả hai kênh khác 0**, không phải hai bên ngang sức.
- Sau khi tắt bơm 100 T0, J_C vẫn còn 0.0024/T0, giảm chậm (có thể do trạng thái gần tối). Tài liệu hai-vai-tro dùng h (entropy) và thấy 0.0000 sau 20 T0. Hai đại lượng này khác nhau.
- Phần cổ điển theo nhiệt độ dùng quỹ đạo ngắn (2 × 40 T0), sai số thống kê khoảng 20–30%.

## Câu hỏi còn mở

- Muốn chứng minh trọn vẹn thì phải dựng một |Ψ⟩⟩ hữu hạn có cả rô-to, hạt và thanh ghi bản ghi, rồi kiểm rằng dòng bản ghi ròng chỉ khác 0 khi cả hai kênh bật. Có làm tiếp không?
- Có nên đo "độ xung đột" bằng một đại lượng thay cho P_QK·P_HT, ví dụ tích của hai tốc độ hoặc độ không giao hoán giữa hai bộ sinh $[\mathcal L_B,\mathcal L_C]$, để tìm đại lượng tỉ lệ đúng với J không?
