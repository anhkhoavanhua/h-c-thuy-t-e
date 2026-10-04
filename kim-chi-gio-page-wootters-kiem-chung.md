# Kiểm chứng số chứng minh "kim chỉ giờ" Page–Wootters của $\mathcal E_{\text{dual}}$

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Part III. Vũ trụ C ⊗ S ⊗ E lấy từ [giai_so_wdw_pw_e7.py](giai_so_wdw_pw_e7.py):
> - S = rô-to ℓ ≤ 4 ⊗ hạt n ≤ 3, tức **36 chiều**;
> - mỗi tích tắc có một ancilla 3 mức {không, Carnot, bơm};
> - τ = T0/25.
>
> **Số liệu:** [giai_so_kim_chi_gio.py](giai_so_kim_chi_gio.py). Chạy khoảng 5 phút trên máy này, vì CPU đang bị giới hạn ở 26% tốc độ.
> **Đầu vào:** chứng minh 4 bước do người dùng dán vào (ràng buộc WdW → chiếu lên kim → Lindblad → mũi tên entropy).

## Tóm tắt

| Bước trong chứng minh | Kiểm số | Kết luận |
|---|---|---|
| 1. $\hat{\mathbb J}\vert\Psi\rangle\rangle=\vert\Psi\rangle\rangle$ tới 4·10⁻¹⁵ | ✓ ra 4.0–6.1·10⁻¹⁵. Nhưng **lịch sử sai cũng thỏa**: 3.0·10⁻¹⁵ khi τ' = 2τ, 3.5·10⁻¹⁵ khi bỏ hẳn H | ⚠ Là **đồng nhất thức của cách dựng**, không phải bằng chứng riêng cho E_dual |
| 2. $\rho_S(k)=\Lambda^k(\rho_0)$ | ✓ đúng nếu mỗi tích tắc dùng **ancilla mới** (mô hình va chạm). Cần chuẩn hóa lại vì ⟨k\|Ψ⟩⟩ có chuẩn 1/√N | ✓ (sửa nhỏ) |
| 3a. Khai triển Kraus | Bậc sai số đo được: τ^1.91 (K0), τ^1.39 (K_j). Lý thuyết là 2 và 1.5 | ✓ bậc đúng |
| | Nhưng ở τ = T0/25: sai số K0 là **3.2** theo chuẩn toán tử và **1.05** ngay trên ψ0, vì τ‖H‖ = 2.9 rad, τΔE = 1.35 rad | ✗ câu "τ = T0/25 rất nhỏ" sai |
| 3b. Hội tụ về Lindblad O(τ) | Sai số giảm **2.04 / 2.01 / 2.00 / 2.00** lần mỗi khi chia đôi τ (t = 1 T0), và 1.91–1.99 lần (t = 5 T0) | ✓ khớp "1.95–2.04". Phải lấy giới hạn **t = kτ cố định**, k → ∞ |
| 4. Mũi tên: $\Delta S(E)/\Delta k\ge0$ là định luật hai | Trong mô hình này S(E) đơn điệu (0/1500 bước giảm). Nhưng **phản ví dụ** (chỉ Carnot, T = 0) cho tới 2489 bước giảm | ✗ không phải định luật hai. Đại lượng đơn điệu đúng là entropy bản ghi thô Σᵢ S(Eᵢ) |
| Kim = bộ đếm Carnot $N_{\text{Carnot}}$ | Mỗi vạch n nhòe **σ_t = 2–7.6 T0**, trong khi cần 0.071 T0 | ✗ $\hat K$ lý tưởng và $N_{\text{Carnot}}$ vật lý **không phải cùng một kim** |

## Bước 1. Ràng buộc Wheeler–DeWitt rời rạc

Để Ĵ đóng vòng, vạch cuối nối về vạch đầu bằng $W=(U_N\cdots U_1)^{-1}$. Do đó **mọi** dãy unita đều thỏa ràng buộc:

| Lịch sử | ‖Ĵ\|Ψ⟩⟩ − \|Ψ⟩⟩‖ |
|---|---|
| E_dual đúng (7 tích tắc, E có 3⁷ chiều) | 4.0·10⁻¹⁵ – 6.1·10⁻¹⁵ |
| Sai: τ' = 2τ | 3.0·10⁻¹⁵ |
| Sai: H' = 0 (không có động cơ) | 3.5·10⁻¹⁵ |

- Ràng buộc thỏa tới sai số máy là **tự động**. Nó cho thấy vũ trụ tĩnh và người quan sát bên ngoài thấy phân bố đều trên C.
- Nó **không** chọn ra động lực E_dual. Động lực nằm ở việc chọn U (H và các L_j), không nằm ở phương trình Ĵ = 1.
- Bản thảo nên viết: "ràng buộc thỏa theo cách dựng". Không nên viết "được kiểm chứng tuyệt đối".

## Bước 3. Từ kim k đến phương trình Lindblad

**3a. Khai triển Kraus**, với $K_j=\langle j\vert U\vert 0\rangle$:

| τ | ‖K0 − (1 − iτH − τ/2 ΣL†L)‖ | ‖K_j + i√τ L_j‖ |
|---|---|---|
| T0/25 | **3.23** | 0.94 |
| T0/50 | 0.98 | 0.43 |
| T0/100 | 0.26 | 0.16 |
| T0/200 | 0.066 | 0.058 |
| T0/400 | 0.017 | 0.021 |

- Bậc hội tụ đúng. Σ K†K − 1 = 1.4·10⁻¹⁴, tức vết được bảo toàn chính xác.
- **Sửa một câu trong chứng minh:**
  - ‖H‖ = 11.6 ħω0, nên τ‖H‖ = 2.9 rad mỗi tích tắc. Phổ cắt ở ℓ ≤ 4, n ≤ 3.
  - Ngay trên trạng thái đầu, τΔE = 1.35 rad.
  - Vì vậy ở τ = T0/25, khai triển bậc một **không** hội tụ từng bước.
- Kết quả cuối vẫn đúng vì phần unita được tính chính xác (U = exp). Sai số O(τ) đến từ việc **tách** H khỏi phần tiêu tán (kiểu Trotter), không đến từ τH nhỏ.

**3b. Hội tụ về e^{tℒ}** (khoảng cách vết; nghiệm tham chiếu là Liouvillian chéo hóa chính xác):

| t | τ = T0/25 | T0/50 | T0/100 | T0/200 | T0/400 | Tỉ số khi chia đôi |
|---|---|---|---|---|---|---|
| 1 T0 | 3.96e-2 | 1.94e-2 | 9.64e-3 | 4.82e-3 | 2.41e-3 | 2.04, 2.01, 2.00, 2.00 |
| 5 T0 | 1.90e-2 | 9.92e-3 | 5.07e-3 | 2.56e-3 | 1.29e-3 | 1.91, 1.96, 1.98, 1.99 |

- Giữ **k = 1** cố định rồi cho τ → 0 thì chỉ được đạo hàm tại t = 0: ‖(Λ−1)ρ/τ − ℒρ‖ = 2.96 → 0.79 → 0.20.
- Phương trình chủ chỉ nổi lên khi **k → ∞ cùng lúc τ → 0, với t = kτ cố định**. Chứng minh nên viết đúng giới hạn này.

## Bước 4. Mũi tên thời gian: S(E) có thật sự đơn điệu?

Với ancilla mới ở mỗi tích tắc, trạng thái SE là thuần, nên **S(E) = S(ρ_S)** theo phân tích Schmidt. Hệ quả: S(E) bị chặn bởi ln d = 3.58.

| k | t/T0 | S(E) = S(ρ_S) | P(≥ 1 Carnot) | Σᵢ S(Eᵢ), các bản ghi riêng lẻ |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | 0 |
| 1 | 0.04 | 0.063 | 0.011 | 0.06 |
| 3 | 0.12 | 0.155 | 0.064 | 0.20 |
| 7 | 0.28 | 0.440 | 0.251 | 0.71 |
| 25 | 1 | 1.269 | 0.613 | 4.37 |
| 125 | 5 | 1.623 | 0.972 | 25.6 |
| 500 | 20 | 1.667 | 1.000 | 106.1 |
| 1500 | 60 | **1.667** (bão hòa) | 1.000 | **320.6** (tuyến tính, 5.36 nat/T0) |

- Các số ở k = 1, 3, 7 **khớp đúng** đoạn chứng minh bạn dán.
- Giá trị bão hòa trong mô hình này là **1.667**, không phải 1.615 như các báo cáo trước. Chưa rõ 1.615 đến từ lần chạy nào (xem câu hỏi còn mở).
- Trong mô hình này S(E) đơn điệu: 0/1500 bước giảm, bước tăng nhỏ nhất 2.6·10⁻¹¹.
- **Phản ví dụ** cùng hệ, chỉ giữ kênh Carnot (làm lạnh Davies ở T = 0), xuất phát từ một trạng thái riêng thuần:

| Trạng thái đầu | S(E) cực đại | Lúc k = 2500 | Số bước ΔS < 0 |
|---|---|---|---|
| Trạng thái riêng thứ 6 | 1.029 (0.4 T0) | 0.014 | **2489** |
| Trạng thái riêng thứ 12 | 1.060 (0.6 T0) | 0.886 | 110 |

**Kết luận bước 4:**
- "ΔS(E)/Δk ≥ 0" **không phải** định luật hai. S(E) của toàn môi trường chỉ là độ vướng víu giữa S và E. Nó tăng rồi **giảm** khi hệ được làm lạnh về trạng thái thuần.
- Tính đơn điệu trong mô hình E_dual có được là nhờ trạng thái dừng là hỗn hợp (do có bơm), chứ không nhờ một nguyên lý chung.
- **Đại lượng đơn điệu thật sự** là entropy bản ghi thô Σᵢ S(Eᵢ):
  - mỗi tích tắc cộng thêm S(E_k) ≥ 0;
  - tăng tuyến tính không giới hạn;
  - Σᵢ S(Eᵢ) ≥ S(E), dấu bằng không xảy ra.
- Đây chính là "đồng hồ entropy của bể" trong [nam-dong-ho-luong-tu-kiem-chung.md](nam-dong-ho-luong-tu-kiem-chung.md).

## Bước 5. Dùng bộ đếm Carnot làm kim

Đoạn dán đồng nhất hai vai trò: $\hat K$ là chuỗi xung lý tưởng τ ≤ T0/25, và $k_{\text{eff}}=N_{\text{Carnot}}$ là kim vật lý. Tôi tách hai vai trò này bằng cách điều kiện hóa ⟨N_C = n\| thay cho ⟨k\|. Phép tính dùng thống kê đếm chính xác ρ_n(k), không lấy mẫu.

Tốc độ bản ghi Carnot là **r = 0.468 mỗi T0**, so với 25 mỗi T0 của kim lý tưởng. Hệ số Fano của N_C ở 60 T0 là 0.81.

| n (số đọc) | ⟨t\|n⟩ ($T_0$) | **σ_t\|n** ($T_0$) | ‖ρ_S\|n − ρ_Lindblad(⟨t⟩)‖ | Kim lý tưởng k ở cùng ⟨t⟩ |
|---|---|---|---|---|
| 0 | 1.56 | 2.02 | **0.317** | 0.036 |
| 1 | 3.15 | 2.76 | **0.100** | 0.022 |
| 2 | 4.83 | 3.38 | 0.031 | 0.019 |
| 5 | 10.8 | 4.87 | 0.025 | 0.021 |
| 10 | 21.5 | 6.54 | 0.022 | 0.022 |
| 20 | 41.8 | 7.63 | 0.021 | 0.022 |

(Hàng n = 30 bị cắt bởi cửa sổ 60 T0 nên không đưa vào bảng.)

- Ở n = 0, 1, trạng thái điều kiện **lệch** khỏi Lindblad 0.1–0.3.
- Từ n ≥ 2, độ lệch bằng mức của kim lý tưởng. Nhưng đó chỉ vì ρ_S đã gần trạng thái dừng: khi hệ không còn đổi, đồng hồ nào cũng "khớp".
- **Phép thử quyết định là độ phân giải.** Mỗi vạch n nhòe **σ_t = 2–7.6 T0**, lớn hơn 0.071 T0 (yêu cầu 99%) tới 30–100 lần và dài hơn cả chu trình 5 T0. Độ nhòe còn tăng theo √n.

## Chứng minh viết lại cho đúng (đề xuất đưa vào bản thảo)

1. **Ràng buộc.** $\hat{\mathbb J}\vert\Psi\rangle\rangle=\vert\Psi\rangle\rangle$ được thỏa **theo cách dựng**: đóng vòng bằng $W=(U_N\cdots U_1)^{-1}$, kiểm số 4·10⁻¹⁵. Ràng buộc này mã hóa tính tĩnh của vũ trụ, không mã hóa động lực. Động lực nằm ở việc chọn $\hat U$.
2. **Kim lý tưởng $\hat K=\sum_k k\vert k\rangle\langle k\vert$** là chuỗi xung đều (đồng hồ Brown–Kuchař hình thức). Chiếu lên ⟨k\| với **ancilla mới mỗi tích tắc** cho $\rho_S(k)=\Lambda^k(\rho_0)$. Trong giới hạn τ → 0, k → ∞ với t = kτ cố định, ta có $\dot\rho=\mathcal L\rho$ với sai số O(τ) (tỉ số 2.00 khi chia đôi τ).
3. **Mũi tên thời gian** không đến từ S(E) toàn cục: S(E) = S(ρ_S), bị chặn bởi ln d và có thể giảm. Nó đến từ **entropy bản ghi thô** $\sum_i S(E_i)$. Đại lượng này tăng đơn điệu với mọi U và mọi ψ0, và tăng tuyến tính (5.36 nat/T0 trong E_dual).
4. **Kim vật lý $N_{\text{Carnot}}$** đơn điệu và không tuần hoàn, nên mang đúng mũi tên. Nhưng nó đọc thời gian với σ_t = 2–7.6 T0. Nó **không** tái lập động lực ở phân giải từng đợt. $\hat K$ lý tưởng và $N_{\text{Carnot}}$ là hai kim khác nhau, đúng như kết quả TUR/Erker trong [nam-dong-ho-luong-tu-kiem-chung.md](nam-dong-ho-luong-tu-kiem-chung.md).

**Câu tiếng Anh cho bản thảo:**

> The discrete Wheeler–DeWitt constraint is satisfied by construction (residual 4×10⁻¹⁵ for any unitary sequence closed by its inverse); dynamics enters only through the choice of the collision unitary. Conditioning on an ideal pulse-train clock with fresh environment ancillas reproduces the Lindblad equation with O(τ) error (ratio 2.00 under τ-halving at fixed t = kτ). The arrow of time is carried not by the global environment entropy, which equals S(ρ_S) and can decrease, but by the coarse-grained record entropy Σᵢ S(Eᵢ), which grows linearly (5.4 nat/T₀). The physical Carnot-record counter is monotonic but resolves time only to σ_t ≈ 2–8 T₀, so it cannot replace the ideal clock at the resolution of individual bursts.

## Câu hỏi còn mở
- Entropy bão hòa: lần chạy này ra 1.667 (ℓ ≤ 4, n ≤ 3, τ = T0/25, 60 T0). Các báo cáo trước ghi 1.615 (Part III N1, rồi [nam-dong-ho-luong-tu-kiem-chung.md](nam-dong-ho-luong-tu-kiem-chung.md)). Cần đối chiếu cấu hình của lần chạy N1 trước khi đưa số vào LaTeX.
