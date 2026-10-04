# Thời gian tự sinh từ góc rô-to (Dirac–Kuchař/WKB): kiểm số

> **Kiểm cái gì:** đề xuất "bỏ đồng hồ ngoài $\hat K$, dùng góc rô-to mở cuộn $\hat\Phi\in\mathbb R$ làm thời gian". Đề xuất này nhằm giải hai giới hạn: sàn đo lường của $\hat K$ và chuẩn hóa phân kỳ của WdW liên tục.
> **Số liệu:** [giai_so_dong_ho_tu_sinh.py](giai_so_dong_ho_tu_sinh.py), kết quả thô [dong_ho_tu_sinh_ket_qua.txt](dong_ho_tu_sinh_ket_qua.txt), chạy khoảng 1.3 phút.
> **Nền:** mở rộng phép đo N3 trong [part3-latex-so-lieu-kiem-chung.md](part3-latex-so-lieu-kiem-chung.md) (đồng hồ rô-to có tác dụng J, sàn 1 − F ≈ 0.05). Nối tiếp [kim-ngoai-dong-ho-vat-ly-kiem-chung.md](kim-ngoai-dong-ho-vat-ly-kiem-chung.md).

## Tóm tắt

| Tuyên bố | Kiểm số | Kết quả |
|---|---|---|
| **I.1** Ràng buộc WdW thỏa tới 4·10⁻¹⁵ | Đúng (4.0–6.1·10⁻¹⁵), nhưng **lịch sử sai cũng thỏa** (3.0·10⁻¹⁵ khi τ′ = 2τ; 3.5·10⁻¹⁵ khi bỏ H) — [kim-chi-gio](kim-chi-gio-page-wootters-kiem-chung.md) | ⚠ đồng nhất thức của cách dựng, không phải bằng chứng riêng cho E_dual |
| **I.2** Lindblad nổi lên với sai số O(τ) | Tỉ số khi chia đôi τ là 2.00–2.04 ở 3 mô hình bể — [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md) | ✓ |
| **I.3** s_k không giảm, h = 5.23 nat/T0 | 0/6·10⁶ bước giảm — [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md) | ✓ |
| **I.4** ⟨ΔE_Carnot⟩ = 0.412 ħω0 | Là **mỗi va chạm** (mô hình chặn gãy p\* = 0.3); cả đợt ≈ 1.3–2.5 ħω0 — [part2 A1–A9](part2-cau-hoi-mo-A1-A9-kiem-chung.md) | ✓ (ghi rõ đơn vị) |
| **I.5** Qua KS ở 5/7 bộ EXFOR | Đúng (D = 0.042–0.066). Nhưng cùng tài liệu đó ghi **exciton chuẩn khớp tốt hơn động cơ ở 7/7 bộ** khi tính đủ cạnh tranh p/α/γ, hiệu ứng cặp và χ² công bằng; KS chưa gộp sai số thực nghiệm | ⚠ phải ghi kèm |
| **III.2** $[\hat T_{\text{geom}},\hat p_\theta]=i\hbar$ | Trên S¹ không có toán tử góc nào: $\min_\Phi\lVert[\Phi,L]-i\rVert_F=\sqrt d$, tức sai 100% (A). Trên ℝ: $[T_{\text{geom}},p]=i\hbar\,I/\langle p\rangle=i\hbar/3$ | ✗ sai hệ số; đúng ra là $[T,H_{\text{rot}}]\approx i\hbar$ gần $\langle p\rangle$ |
| **III.3** Gần đúng WKB/đoạn nhiệt cho Schrödinger của hạt | Ở J ≈ 3ħ của học thuyết: 1 − F = 0.054 (1 vòng), **0.31** (WdW, đồng hồ S¹) | ⚠ thô ở J của học thuyết |
| **IV.2** Mở cuộn Φ ∈ ℝ "khắc phục hoàn toàn" sàn 1 − F ≈ 0.05 | Xóa được phần **tái hồi** (0.31 → 0.043 ở J = 3). Còn lại sàn **phản lực** 0.043–0.060 ở J = 3. Ở giới hạn WdW (tw = 50 T_rot), sàn 0.048 của J = 30 **không đổi** khi chu kỳ đồng hồ còn ≤ 8 T_rot (B) | ✗ không "hoàn toàn"; còn đổi lại chuẩn hóa phân kỳ |
| **IV.3** "x vượt bán kính giam r_br ≈ 6.11 fm" khi gãy | Trên đa tạp $\lvert x\rvert\le\kappa=0.679$ fm = **0.11 r_br**. Gãy là ngưỡng **lực** $\lvert\lambda\rvert\ge R=k\,r_{\text{br}}$, không phải ngưỡng vị trí (C) | ✗ đọc sai ý nghĩa của r_br |

## A. Toán tử thời gian góc

**Trên S¹** (mômen $L=\hbar\ell$, ℓ nguyên), trong cơ sở ℓ ta có $[\Phi,L]_{ab}=\Phi_{ab}(\ell_b-\ell_a)$. Đường chéo luôn bằng 0 với **mọi** Φ, nên

$$\min_\Phi\lVert[\Phi,L]-i\mathbb 1\rVert_F=\sqrt d$$

| ℓ ∈ | d | min‖[Φ, L] − i·1‖_F | Tương đối |
|---|---|---|---|
| [−3, 3] | 7 | 2.6458 = √7 | 100% |
| [−10, 10] | 21 | 4.5826 = √21 | 100% |
| [−30, 30] | 61 | 7.8102 = √61 | 100% |

**Trên ℝ** (góc mở cuộn): $[\hat\Phi,\hat p_\theta]=i\hbar$, nên
$$[\hat T_{\text{geom}},\hat p_\theta]=\frac{I}{\langle p_\theta\rangle}\,i\hbar=\frac{i\hbar}{3}\quad(\langle p_\theta\rangle\approx3\hbar,\ I=1),$$
không phải iħ. Cặp liên hợp có ý nghĩa là $[\hat T_{\text{geom}},\hat H_{\text{rot}}]=i\hbar\,\hat p/\langle p\rangle$, và nó chỉ bằng iħ gần $\langle p\rangle$.

Định lý Pauli cũng áp dụng: $\hat H$ bị chặn dưới nên không có $\hat T$ tự liên hợp nào thỏa $[\hat T,\hat H]=i\hbar$ chính xác.

## B. "Mở cuộn" có xóa được sàn 1 − F không?

**Cách đo.** Rô-to được đặt trên một vòng có chu vi $2\pi/\delta$, nên mômen là bội của δ:
- δ = 1 là S¹ của học thuyết;
- δ → 0 là đường thẳng, tức "Φ ∈ ℝ";
- chu kỳ của đồng hồ là $T_{\text{rot}}/\delta$.

Trạng thái WdW được lọc qua cửa sổ năng lượng có độ rộng nghịch đảo tw. Đo độ trùng giữa trạng thái hạt điều kiện theo góc và nghiệm Schrödinger với θ = Ωt (lấy min trên t ∈ {0.1, …, 0.4} T_rot). Thông số giống N3: g = 1, Ω = 0.3.

**Kiểm lại N3:** với δ = 1, J = 3 cho 0.054 / 0.306 và J = 30 cho 0.015 / 0.048, khớp N3 gốc (0.054 / 0.307 và 0.015 / 0.048).

| J (ħ) | δ | Chu kỳ đồng hồ | 1 − F (tw = 0.35 T_rot) | (tw = 5) | (tw = 50, WdW) |
|---|---|---|---|---|---|
| 3 | 1 (S¹) | 1 T_rot | 0.054 | 0.255 | **0.307** |
| 3 | 1/2 | 2 T_rot | 0.060 | 0.043 | **0.043** |
| 3 | 1/4 | 4 T_rot | 0.060 | 0.043 | 0.043 |
| 3 | 1/8 | 8 T_rot | 0.055 | 0.059 | 0.059 |
| 30 | 1 (S¹) | 1 T_rot | 0.015 | 0.048 | **0.048** |
| 30 | 1/2 | 2 T_rot | 0.003 | 0.022 | 0.048 |
| 30 | 1/4 | 4 T_rot | 0.004 | **0.006** | **0.048** |
| 30 | 1/8 | 8 T_rot | 0.004 | 0.006 | 0.048 |

Sàn 1 − F gồm **hai phần** khác nhau:

1. **Tái hồi (do tô-pô):** trạng thái tại t bị trộn kết hợp với trạng thái tại $t+nT_{\text{clock}}$. Mở cuộn xóa được phần này, **nhưng chỉ khi chu kỳ đồng hồ dài hơn độ dài kết hợp tw**:
   - J = 30, tw = 5: sàn giảm từ 0.048 xuống 0.006 khi chu kỳ đạt 4 T_rot.
   - J = 30, tw = 50 (giới hạn WdW): sàn vẫn **0.048** dù chu kỳ đã là 8 T_rot.
   - Trong WdW thật (tw → ∞), mọi đồng hồ compact đều giữ sàn này. Chỉ đồng hồ **phi compact** thật sự (δ → 0) mới xóa được.
2. **Phản lực / đoạn nhiệt:** phần còn lại khi đã mở cuộn.
   - J = 3, tức J của học thuyết (tác dụng rô-to ≈ 3.2ħ): **0.043–0.060** ở mọi δ < 1.
   - Sàn này bằng đúng cỡ "0.05" mà đề xuất muốn xóa. Nó đến từ việc bỏ $\partial_\phi^2\chi$ (gần đúng WKB) khi J nhỏ, nên mở cuộn không chạm tới.
   - J = 30: chỉ còn khoảng 0.004.

## C. Biên độ x và bán kính r_br

| Đại lượng | Giá trị |
|---|---|
| $r_{\text{br}}=3\ell_0$ | 6.11 fm |
| $\lvert x\rvert$ trên đa tạp, $\le\kappa=\ell_0/3$ | ≤ **0.679 fm** = 0.11 r_br |
| Dao động điểm không $x_{\text{rms}}=\ell_0/\sqrt2$ | 1.44 fm |

$r_{\text{br}}$ là bán kính tại đó lực AdS $m\omega_0^2r$ bằng R ([ads-cft-mo-rong](ads-cft-mo-rong-kiem-chung.md)). Gãy xảy ra khi **phản lực** $\lvert\lambda\rvert$ đạt R, trong khi hạt vẫn ở trong vùng κ = 0.68 fm. Câu "biên độ x vượt ra ngoài r_br" là sai.

## Kết luận

1. **Phần I (năm số liệu cũ):** cả năm đều có nguồn và đúng, nhưng phải kèm hai lưu ý.
   - Con số 4·10⁻¹⁵ là đồng nhất thức của cách dựng: lịch sử sai cũng thỏa ràng buộc.
   - Câu "KS 5/7" phải đi cùng kết quả: exciton chuẩn khớp tốt hơn ở 7/7 bộ.
2. **Toán tử thời gian góc:**
   - Trên S¹ không tồn tại.
   - Trên ℝ thì hệ thức đúng là $[T_{\text{geom}},p_\theta]=i\hbar I/\langle p\rangle$, không phải iħ.
   - Định lý Pauli cấm một T tự liên hợp chính xác.
3. **Mở cuộn Φ ∈ ℝ không "khắc phục hoàn toàn" sàn 0.05.**
   - Nó xóa phần tái hồi, và chỉ khi đồng hồ thật sự phi compact.
   - Đổi lại, đồng hồ phi compact có phổ liên tục, nên $\langle\langle\Psi|\Psi\rangle\rangle$ phân kỳ. Đó chính là "Giới hạn 2" mà đề xuất muốn né.
   - Ở J ≈ 3ħ của học thuyết, sàn phản lực 0.04–0.06 vẫn còn.
4. **Có một sự đánh đổi tô-pô:**
   - Đồng hồ compact (S¹): chuẩn hóa được, nhưng tuần hoàn nên có sàn tái hồi.
   - Đồng hồ phi compact (ℝ): không có tái hồi, nhưng không chuẩn hóa được.
   - Bộ đếm vòng ghi bản ghi (V.C.4, bể Carnot) đứng giữa hai thái cực. Vì vậy học thuyết vẫn cần bản ghi Carnot như một bộ đếm, không thể dựa riêng vào góc rô-to.
5. **Phần đứng được của đề xuất:** ý "động cơ vừa sinh công vừa là đồng hồ hình học" đúng **ở J lớn** (J = 30: 1 − F ≈ 0.004 khi mở cuộn và tw hữu hạn). Nhưng hạt nhân của học thuyết có J ≈ 3ħ, nên đó là một đồng hồ thô.

## Câu hỏi còn mở

- Có muốn thử phương án "rô-to + bộ đếm vòng thuận nghịch" (Φ = θ + 2πn với n là một thanh ghi) không? Phương án này có thể né được sàn tái hồi mà vẫn giữ chuẩn hóa nếu thanh ghi hữu hạn. Cái giá là thanh ghi tự nó lại là một đồng hồ.
- Muốn giảm sàn phản lực ở J ≈ 3 thì cần giữ $\partial_\phi^2\chi$ (sửa lỗi sau WKB, kiểu Born–Oppenheimer bậc hai). Có muốn giải bậc tiếp theo không?
