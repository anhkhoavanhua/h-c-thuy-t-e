# Wheeler–DeWitt "= 0": Học thuyết E tính ra được gì, và giải tiếp tới đâu

> **Câu hỏi:** phương trình WdW $\hat{\mathcal H}|\Psi\rangle\rangle=0$ làm vũ trụ "đóng băng". Giải qua toán thì học thuyết $\mathcal E_{\text{dual}}$ thu được gì? Nếu chưa đủ, giải tiếp để tìm điểm đột phá.
> **Số liệu mới:** [giai_so_thoi_gian_tu_sinh_bac2.py](giai_so_thoi_gian_tu_sinh_bac2.py), kết quả thô [thoi_gian_tu_sinh_bac2_ket_qua.txt](thoi_gian_tu_sinh_bac2_ket_qua.txt). Chạy theo từng phần: `python3 giai_so_thoi_gian_tu_sinh_bac2.py [edual|clock|width|comm|mech]`, mỗi phần 20–40 giây; chạy không tham số là bảng quét g theo J.
> **Nối tiếp:** [dong-ho-tu-sinh-goc-roto-kiem-chung.md](dong-ho-tu-sinh-goc-roto-kiem-chung.md), [kim-ngoai-dong-ho-vat-ly-kiem-chung.md](kim-ngoai-dong-ho-vat-ly-kiem-chung.md), [hai-vai-tro-kim-va-be-kiem-chung.md](hai-vai-tro-kim-va-be-kiem-chung.md).

## 1. "= 0" không có nghĩa là "không có gì"

$\hat{\mathcal H}|\Psi\rangle\rangle=0$ chỉ nói rằng **tổng** năng lượng (gồm cả đồng hồ) bằng 0. Toàn bộ nội dung vật lý nằm ở **tương quan giữa các phần**: điều kiện theo đồng hồ, bản ghi, nhánh. Đây là cách đọc Page–Wootters, và Học thuyết E đã khai thác nó.

## 2. Những gì học thuyết đã tính ra từ phương trình ràng buộc (có nguồn)

| # | Đại lượng | Giá trị | Lưu ý | Nguồn |
|---|---|---|---|---|
| 1 | Ràng buộc $\hat{\mathbb J}\lvert\Psi\rangle\rangle=\lvert\Psi\rangle\rangle$ | 4.0–6.1·10⁻¹⁵ | ⚠ Đồng nhất thức của cách dựng; lịch sử sai cũng thỏa | [kim-chi-gio](kim-chi-gio-page-wootters-kiem-chung.md) |
| 2 | Lindblad nổi lên khi điều kiện theo kim | sai số O(τ), tỉ số 2.00–2.04 | Đúng với mọi mô hình bể | [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md) |
| 3 | Mũi tên: surprisal bản ghi $s_k$ | 0/6·10⁶ lần giảm; h = 5.23 nat/T0 | S(E) bão hòa 1.667; mũi tên tắt khi động cơ tắt | [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md) |
| 4 | Năng lượng mỗi va chạm Carnot | 0.412 ħω0 | Mỗi va chạm, không phải mỗi đợt | [part2 A1–A9](part2-cau-hoi-mo-A1-A9-kiem-chung.md) |
| 5 | Phổ neutron EXFOR | KS qua ở 5/7 bộ | ⚠ Exciton chuẩn khớp tốt hơn ở 7/7 bộ | [part2 A1–A9](part2-cau-hoi-mo-A1-A9-kiem-chung.md) |
| 6 | Hai vai trò | Bể cho **chiều**; $\hat K$ cho **số đo** | Bể làm đồng hồ thô: σ_t ≈ 0.8–5.6 T0 | [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md) |
| 7 | $\hat K$ vật lý bên ngoài | cần ≤ 29 ys (xếp thứ tự), ≤ 1.2 ys (trạng thái) | Xung attosecond chậm 10⁶ lần; rung tích lũy thêm $\frac D2\mathcal L^2$ | [kim-ngoai](kim-ngoai-dong-ho-vat-ly-kiem-chung.md) |
| 8 | Góc rô-to mở cuộn làm đồng hồ | Xóa phần tái hồi, còn sàn phản lực | $[T_{\text{geom}},p]=i\hbar I/\langle p\rangle$; đồng hồ phi compact thì không chuẩn hóa được | [dong-ho-tu-sinh](dong-ho-tu-sinh-goc-roto-kiem-chung.md) |

**Điểm nghẽn còn lại:** muốn dựng được thời gian ($\hat K$) thì vẫn cần một đồng hồ **lý tưởng ở ngoài**. Câu hỏi đột phá: liệu chính động cơ (góc rô-to θ) có thay được $\hat K$ ở **đúng thông số của học thuyết** không?

## 3. Giải tiếp: đồng hồ rô-to ở đúng thông số E_dual

Các phép đo trước (N3, file dong-ho-tu-sinh) dùng mô hình thử với **g = 1, Ω = 0.3**. Học thuyết (Bảng I) có **g = 0.1, I = 1, θ̇ ≈ 3ω0, ΔE = 1, J ≈ 3ħ**. Phép đo ở đây: lấy trạng thái WdW (cửa sổ năng lượng tw = 50 T_rot, gần như một mặt năng lượng), điều kiện theo góc rô-to φ, rồi so với nghiệm Schrödinger của hạt khi rô-to cổ điển quay θ_c(t). Lấy max 1 − F trên 8 góc trong vòng đầu.

### 3.1 Sàn "0.05" là do g = 1 của mô hình thử

Hạt xuất phát ở trạng thái cơ bản, tức chỉ có chuyển động **cưỡng bức** bởi rô-to:

| g | 1 − F, S¹, WdW | 1 − F, mở cuộn ×4 | Nền: coi hạt đứng yên |
|---|---|---|---|
| **0.1 (E_dual)** | **0.0013** | 0.0014 | 0.0032 |
| 0.3 | 0.0113 | 0.0122 | 0.0282 |
| 1.0 (mô hình N3) | 0.1215 | 0.1169 | 0.2720 |

- Sai số theo luật **1 − F ≈ 0.12 g²**.
- Ở g = 0.1, mở cuộn **không cần thiết**: S¹ và mở cuộn ×4 cho cùng kết quả.
- Nhưng mức nền cũng nhỏ (0.0032). Đồng hồ rô-to chỉ thu được khoảng **55–60%** độ lệch cưỡng bức (1 − F nhỏ hơn nền khoảng 2.3 lần ở cả ba giá trị g).

### 3.2 Hạt chuyển động tự do: đồng hồ rô-to thất bại

Hạt xuất phát ở trạng thái kết hợp α = 0.707i (ψ0 của Part III), nên nó dao động với tần số ω0 riêng của mình:

| g | 1 − F, S¹, WdW | 1 − F, mở cuộn ×4 | Nền |
|---|---|---|---|
| 0.1 | **0.424** | 0.398 | 0.820 |
| 0.3 | 0.485 | 0.445 | 0.867 |
| 1.0 | 0.707 | 0.576 | 0.962 |

Ở thông số của học thuyết, góc rô-to **không đo được** dao động riêng của hạt: sai 42%.

### 3.3 Tìm nguyên nhân: quét tốc độ rô-to θ̇ và độ rộng bó sóng σ_ℓ

Các điều kiện chung: J = 30, ΔE = 0, α = 0.707i, g = 0.1, rô-to S¹.

| θ̇/ω0 | 3 | 2 | 1.5 | **1** | 0.7 | 0.5 | 1/3 | 0.3 | 0.25 |
|---|---|---|---|---|---|---|---|---|---|
| max 1 − F (tw = 50) | 0.416 | 0.435 | 0.469 | **0.056** | 0.312 | 0.354 | 0.371 | 0.374 | 0.378 |

- Ở J = 3, 10, 30 (σ_ℓ = 0.7): chỉ θ̇ = ω0 thành công (0.08–0.09 khi J ≥ 10, 0.27 khi J = 3). Các giá trị khác đều ≈ 0.37–0.42.
- **σ_ℓ không quyết định:** tăng từ 0.7 lên 10 không đổi kết quả khi θ̇ = 3 hoặc 0.3. Với θ̇ = 1, kết quả giảm từ 0.09 xuống 0.056.
- **Giả thuyết "khớp số học" k·θ̇ = ω0 bị bác bỏ:** θ̇ = 1/2 và 1/3 vẫn thất bại.
- Kết luận thực nghiệm: trong WdW nghiêm ngặt, đồng hồ rô-to S¹ chỉ tái lập dao động riêng của hạt khi **cộng hưởng θ̇ ≈ ω0**. Học thuyết có θ̇ ≈ 3ω0, nằm ngoài cộng hưởng.

## 4. Điểm đột phá rút ra (kết quả mô hình, có giới hạn)

**Thời gian tự sinh từ động cơ chỉ tồn tại cho phần chuyển động bị động cơ khóa.**

| Pha | Chuyển động của hạt | Đồng hồ rô-to (thông số E_dual) | Thời gian đến từ đâu |
|---|---|---|---|
| **Quá khứ** $\mathcal E_{\text{past}}$ (khóa, $x=\kappa\sin\theta$) | Bị rô-to kéo theo | Đo được (1 − F ≈ 10⁻³, sai số ∝ g²) | **Góc rô-to**: thời gian tự sinh |
| **Hiện tại** $\mathcal E_{\text{present}}$ (gãy, x tự do) | Dao động riêng ở ω0 ≠ θ̇ | **Thất bại** (1 − F ≈ 0.42) | Không thể từ rô-to. Phải từ **bản ghi Carnot** (bể) hoặc $\hat K$ ngoài |
| **Tương lai** $\mathcal E_{\text{future}}$ (tái khóa) | Lại bị khóa | Đo lại được | Góc rô-to, cộng thêm số bản ghi để đếm vòng |

Ý nghĩa:
- WdW "= 0" không làm vũ trụ của học thuyết đóng băng hoàn toàn. Trong pha khóa, động cơ **tự là đồng hồ**.
- Nhưng đúng lúc hạt thoát khỏi liên kết (gãy), đồng hồ nội tại **mất khả năng đo**. Đó cũng là lúc bể Carnot bắt đầu ghi bản ghi, tức lúc mũi tên thời gian sinh ra.
- Như vậy, ba pha của học thuyết trùng với **ba chế độ của thời gian**:
  - đo bằng hình học (rô-to);
  - đo bằng nhiệt động (bản ghi);
  - rồi quay về đo bằng hình học.
- Đây là một liên kết **mới** giữa cơ học ba pha của Part I và vấn đề thời gian của Part III, có số đo hỗ trợ.

**Câu tiếng Anh đề xuất:**

> Solving the timeless constraint at the model's own parameters (g = 0.1, θ̇ ≈ 3ω₀, J ≈ 3ħ), the rotor angle serves as an internal clock only for the motion it locks: the driven response is reproduced to 1 − F ≈ 10⁻³ (scaling as 0.12 g²), whereas the particle's free oscillation is not (1 − F ≈ 0.42), except at the resonance θ̇ ≈ ω₀. Self-emergent time therefore exists in the synchronized Past phase but is lost at breakdown, precisely when the Carnot bath begins to write records. The three phases correspond to geometric, thermodynamic and again geometric time-keeping.

## 5. Giới hạn (cần ghi rõ nếu đưa vào bài)

- Mô hình WdW ở đây là rô-to + dao động tử **không có bể**. Câu "bể tiếp quản thời gian ở pha Hiện tại" là suy luận nối với [hai-vai-tro](hai-vai-tro-kim-va-be-kiem-chung.md), chưa giải trong cùng một trạng thái WdW có cả rô-to, hạt và bể.
- Trường hợp "khóa" được mô phỏng bằng ghép yếu với hạt xuất phát ở trạng thái cơ bản, chứ không phải ràng buộc cứng $x=\kappa\sin\theta$. Trong trường hợp này mức nền nhỏ (0.0032), nên kết quả 10⁻³ chỉ tốt hơn nền khoảng 2.3 lần.
- Cơ chế đầy đủ của cộng hưởng θ̇ ≈ ω0 chưa được giải thích: giả thuyết khớp số học đã bị bác bỏ.
- Mọi kết quả dùng $J\le30$ và cửa sổ năng lượng Gauss hữu hạn (tw = 50 T_rot) thay cho mặt năng lượng chính xác.


## 6. Giải tiếp: cơ chế của cộng hưởng θ̇ ≈ ω0

Chạy `python3 giai_so_thoi_gian_tu_sinh_bac2.py mech`, cùng một đoạn quét tâm cộng hưởng (kết quả trong file txt).

**Giả thuyết.** Trên một mặt năng lượng WdW với phổ rời rạc, hạt chỉ "chuyển động" theo góc rô-to khi các cặp $(\ell,n)$ và $(\ell-h,n+1)$ **gần suy biến** và **được tương tác nối trực tiếp**. Ghép $\sin h\theta$ dời ℓ đi ±h, nên điều kiện cộng hưởng là $h\dot\theta\approx\omega_0$. Ngoài cộng hưởng, mỗi trạng thái riêng gần như tích $|\ell\rangle|n\rangle$, nên trạng thái hạt điều kiện theo góc **đứng yên**. Đây chính là "đóng băng thời gian" dưới dạng rời rạc.

**(a) Đổi bậc điều hòa của ghép, h** (J = 30, σ_ℓ = 5, g = 0.1, α = 0.707i, S¹, tw = 50 T_rot):

| Ghép | θ̇ = ω0 | θ̇ = ω0/2 | θ̇ = ω0/3 |
|---|---|---|---|
| sin θ (h = 1) | **0.056** | 0.354 | 0.371 |
| sin 2θ (h = 2) | 0.436 | **0.140** | — |
| sin 3θ (h = 3) | — | — | **0.264** |

Điểm thành công dịch theo h đúng như dự đoán. Tuy vậy, h càng cao thì độ trùng càng kém.

**(b) Tâm và độ rộng của cộng hưởng (ghép sinθ):**

| g | Tâm θ̇*/ω0 | 1 − F tại tâm | Vùng 1 − F < 0.05 | Độ rộng |
|---|---|---|---|---|
| 0.05 | 0.97 | 0.018 | [0.93, 1.00] | 0.07 |
| 0.1 | 0.95 | **0.008** | [0.91, 0.99] | 0.08 |
| 0.3 | 0.85 | **0.004** | [0.825, 0.90] | 0.075 |

- **Tâm dịch xuống khoảng 0.5g.** Độ dịch này cùng cỡ yếu tố ma trận nối hai trạng thái gần suy biến, $g\langle n{+}1|x|n\rangle\langle\ell{-}1|\sin\theta|\ell\rangle\approx0.35\text{–}0.5\,g$. Cách đọc khả dĩ: rô-to và hạt lai thành các mode chuẩn, và đồng hồ bám tần số của mode lai chứ không phải ω0 trần. Đây là diễn giải, chưa chứng minh.
- **Độ rộng gần như không đổi** (0.07–0.08), **không** tỉ lệ g. Nửa này của dự đoán bị bác bỏ. Độ rộng có thể do cửa sổ năng lượng hoặc do J quyết định; chưa kiểm.
- Ở đúng tâm, ghép mạnh hơn cho độ trùng tốt hơn (0.018 → 0.008 → 0.004).

**Hệ quả cho E_dual.** Rô-to của học thuyết quay ở θ̇ ≈ 3ω0 với ghép bậc h = 1. Như vậy nó lệch cộng hưởng gấp 3 lần, và cộng hưởng chỉ rộng khoảng 0.08ω0. Dao động riêng của hạt **không thể** được đồng hồ rô-to đo. Muốn đo được thì cần θ̇ ≈ 0.95ω0, thấp hơn nhiều so với điểm hút ở pha khóa của Part I (θ̇_d ≈ 0.82–0.89 θ̇_br ≈ 2.6–2.9ω0). Kết luận ở mục 4 vẫn đứng, và nay đã có cơ chế: **thời gian tự sinh cần cộng hưởng phổ giữa đồng hồ và hệ, và động cơ của học thuyết không cộng hưởng với dao động riêng của hạt.**

## 7. Đối chiếu một bài phân tích bên ngoài về file dong-ho-tu-sinh

| Ý trong bài phân tích | Đối chiếu với số liệu mới | Kết luận |
|---|---|---|
| A, C, I: toán tử góc sai; r_br bị đọc sai; 4·10⁻¹⁵ là đồng nhất thức; KS 5/7 kèm exciton chuẩn tốt hơn | Khớp file [dong-ho-tu-sinh](dong-ho-tu-sinh-goc-roto-kiem-chung.md) | ✓ |
| Đánh đổi tô-pô compact ↔ phi compact | Khớp | ✓ |
| "Sàn 0.04–0.06 là nhiễu lượng tử tự nhiên của hạt nhân ở J ≈ 3ħ" | Sàn đó đo ở **g = 1** của mô hình thử. Ở g = 0.1 của học thuyết chỉ còn **0.0013**, và sai số theo luật ≈ 0.12g² (mục 3.1) | ✗ **lỗi thời**: không phải đặc tính của hạt nhân |
| "Sàn chỉ biến mất khi J ≫ 30ħ" | Sàn do g quyết định, không do J: J = 3 với g = 0.1 cho 0.0013 | ✗ |
| "Hướng B (Born–Oppenheimer bậc hai) kéo sàn xuống" | Ở g của học thuyết, sàn phản lực đã nhỏ. Nút thắt thật là **cộng hưởng phổ** (dao động riêng sai 42%), thứ mà BO bậc hai không chạm tới. Ehrenfest tự hợp còn làm **tệ hơn** (0.064 so với 0.043 ở g = 1, J = 3) | Hướng B **không phải** đường đột phá |
| "Mỗi chu trình gãy–tái khóa ghi nhận một lượng tử 0.412 ħω0" | 0.412 là **mỗi va chạm**; cả đợt ≈ 1.3–2.5 ħω0 (khoảng 1.5–5 va chạm) | ✗ nhầm đơn vị |
| "Rô-to là hàng đơn vị, bản ghi Carnot là hàng chục" | Hợp lý cho pha khóa. Nhưng mục 3.2 và 6 cho thấy rô-to còn **không đo được** dao động riêng của hạt; khi đó bản ghi phải làm **toàn bộ** việc đo, không chỉ đếm vòng | ⚠ cần sửa: phân vai theo **pha**, không chỉ theo chữ số |
| Hướng A: rô-to S¹ + thanh ghi Carnot | Đúng hướng, và là bước kế tiếp: giải WdW có cả rô-to, hạt và ancilla bản ghi | ✓ chưa giải |

## Câu hỏi còn mở

- Giải một trạng thái WdW **có cả bể** (rô-to + hạt + vài ancilla bản ghi), để kiểm trực tiếp rằng sau khi gãy, việc điều kiện theo số bản ghi khôi phục được thời gian mà rô-to đã mất. Đây là bước kế tiếp tự nhiên.
- Đã giải một phần ở mục 6: cộng hưởng hθ̇ ≈ ω0 dịch khoảng 0.5g (cách đọc khả dĩ: mode lai). Còn mở: độ rộng khoảng 0.08 do đâu (cửa sổ năng lượng hay J)?
