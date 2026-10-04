# Năm loại đồng hồ lượng tử áp vào động cơ $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf). $\hbar\omega_0$ = 10 MeV, $T_0$ = 4.14·10⁻²² s, A ~ 70.
> **Số liệu:** [giai_so_nam_dong_ho.py](giai_so_nam_dong_ho.py), chạy vài giây. Dùng lại số liệu của Part III, mô hình chặn gãy và thế quang học KD.
> **Đòi hỏi chung** ([dong-ho-hat-nhan-page-wootters-kiem-chung.md](dong-ho-hat-nhan-page-wootters-kiem-chung.md)): xếp đúng 99% ba pha khi $\sigma_t\le$ 0.071 $T_0$; xếp đúng 90% khi $\sigma_t\le$ 0.27 $T_0$. Một chu trình khoảng 5 $T_0$.

## Tóm tắt

| Đồng hồ | Phép kiểm số | Dùng được cho học thuyết? |
|---|---|---|
| 1. Salecker–Wigner | Một **nucleon** giữ được độ chính xác 99% chỉ trong **3 $T_0$**; cả hạt nhân A = 70 giữ được 208 $T_0$ | Chỉ như **giới hạn**: đồng hồ phải nặng cỡ cả hạt nhân |
| 2. WKB bán cổ điển | Pha WKB của rô-to đạt chuẩn: \|d(1/p)/dθ\| ≤ 0.075. Nhưng J ≈ 2.4ħ nên phản lực lớn (1−F = 5–31%) và đồng hồ tuần hoàn | Gần đúng trong **một vòng**, không đếm vòng |
| 3. Brown–Kuchař (bụi) | Quãng đường tự do của nucleon 1.5–24 fm ≈ R (4.9 fm): hạt nhân "bán trong suốt" | Chỉ là **đồng hồ lý tưởng hình thức** |
| 4. Entropy | Entropy **của hệ bão hòa** sau ~20 $T_0$ (1.615 nat). Entropy **của bể** tăng tuyến tính mãi (2.5–10 $k_B$ mỗi $T_0$) | ✓ Dùng **entropy của bể**, không dùng entropy của hệ |
| 5. Tự hành (Erker) + TUR | Đồng hồ Carnot đạt 12–21% giới hạn TUR. Để σ_t ≤ 0.071 $T_0$ phải tỏa ≥ **1984 MeV mỗi $T_0$**, thiếu 66–99 lần | ✓ Đúng loại đồng hồ của học thuyết, nhưng **thô** do ngân sách năng lượng |

## 1. Đồng hồ Salecker–Wigner: $Mc^2\gtrsim\hbar T/\delta t^2$

Với $\hbar/T_0$ = 1.59 MeV:

| δt ($T_0$) | T ($T_0$) | $Mc^2$ tối thiểu | Số nucleon | Phần của hạt nhân A = 70 |
|---|---|---|---|---|
| 0.27 | 5 | 109 MeV | 0.1 | 0.002 |
| 0.27 | 300 | 6.6 GeV | 7 | 0.10 |
| 0.071 | 5 | 1.6 GeV | 1.7 | 0.02 |
| 0.071 | 50 | 15.8 GeV | 17 | 0.24 |
| 0.071 | 300 | **94.7 GeV** | **101** | **1.44** |

- Một nucleon (939 MeV) giữ độ chính xác 99% trong **3 $T_0$**, và độ chính xác 90% trong 43 $T_0$.
- Cả hạt nhân A = 70 giữ độ chính xác 99% trong 208 $T_0$.
- **Ý nghĩa:** muốn theo dõi ba pha trong ~300 chu kỳ, đồng hồ phải **nặng hơn chính hạt nhân**. Đồng hồ nằm trong hạt nhân không thể vừa nhẹ vừa chính xác lâu dài.

## 2. Đồng hồ WKB bán cổ điển

Tại $\dot\theta_d$ = 2.63, xung lượng rô-to dọc quỹ đạo là $p_\theta$ = 2.0–2.9 ħ. Bước sóng de Broglie theo góc 2.2–3.1 rad.

- **Tiêu chí pha WKB đạt:** $|d(1/p)/d\theta|$ max 0.075, trung bình 0.048, vì thế thay đổi chậm dọc quỹ đạo.
- **Nhưng hai điều kiện khác của "đồng hồ nặng" không đạt** (số liệu Part III, N3):
  - **phản lực:** với J = 3ħ, 1−F = 0.054 trong một vòng, và **0.31** ở giới hạn WdW;
  - **tuần hoàn:** sàn sai số do trộn vòng (0.014–0.05) không giảm khi J tăng.
- **Ý nghĩa:** thời gian WKB có nổi lên từ rô-to, nhưng chỉ chính xác vài % trong một vòng. Nó không đếm vòng và không đủ nặng để bỏ qua phản lực.

## 3. Đồng hồ Brown–Kuchař: nucleon có là "bụi" không?

Quãng đường tự do trung bình tính từ thế quang học KD 2003: λ = ħv/(2W).

| E (MeV) | $W_V$ | $W_V+W_D$ | λ (thể tích) | λ (cả bề mặt) | λ/R |
|---|---|---|---|---|---|
| 5 | 0.43 | 6.76 | 23.7 fm | 1.5 fm | 0.3–4.8 |
| 10 | 0.76 | 7.63 | 18.9 fm | 1.9 fm | 0.4–3.8 |
| 20 | 1.62 | 8.10 | 12.5 fm | 2.5 fm | 0.5–2.5 |
| 40 | 3.75 | 8.33 | 7.7 fm | 3.5 fm | 0.7–1.6 |

- Phía trong, nucleon đi được gần như xuyên hạt nhân. Ở bề mặt, nó bị hấp thụ sau 1.5–3.5 fm.
- Một nucleon 20 MeV đi hết bán kính trong khoảng 0.2 $T_0$. Trong vài $T_0$ nó chạm bề mặt nhiều lần và va chạm, nên **không giữ được "thời gian riêng" của mình** như bụi.
- **Ý nghĩa:** đồng hồ bụi Brown–Kuchař trong học thuyết chỉ dùng được như **đồng hồ lý tưởng hình thức** (Ĥ_C = p̂). Đó chính là đồng hồ đã dùng ở Part III, chứ không phải các nucleon thật.

## 4. Đồng hồ entropy

| Đại lượng | Theo thời gian | Làm đồng hồ? |
|---|---|---|
| Entropy hệ S(ρ_S) (Part III) | 0 → 1.25 (1 $T_0$) → 1.57 (5 $T_0$) → **1.615 (20 $T_0$) → 1.615 (60 $T_0$)** nat | ✗ bão hòa, đồng hồ dừng |
| Entropy bể, mô hình V.C.1 | 0.50 bản ghi/$T_0$ × 20.5 $k_B$ = **10.2 $k_B$/$T_0$** | ✓ tăng tuyến tính mãi |
| Entropy bể, mô hình chặn gãy (khớp Part I) | 0.61 × 4.1 $k_B$ = **2.5 $k_B$/$T_0$** | ✓ |

**Ý nghĩa:** "đồng hồ entropy" của học thuyết phải đo **entropy đã đổ vào bể Carnot**, tức số bản ghi nhân ΔE/k_BT, chứ không phải entropy von Neumann của động cơ. Entropy của động cơ chạm trần sau khoảng 20 $T_0$, khi hệ vào chu trình dừng.

## 5. Đồng hồ tự hành (Erker 2017) và bất đẳng thức bất định nhiệt động

TUR giới hạn độ chính xác: $N=1/\mathrm{CV}^2\le\Sigma_{\text{tích tắc}}/2k_B$. Mỗi tích tắc Carnot tỏa $\Sigma\approx$ 4.1 $k_B$, nên giới hạn là 2.05.

| Đồng hồ | N đo được | % giới hạn TUR |
|---|---|---|
| Chỉ bản ghi Carnot | 0.24 | 12% |
| Mọi bản ghi | 0.43 | 21% |

**Năng lượng cần để đạt yêu cầu** (r tích tắc mỗi $T_0$, sai số $\sigma_t=\sqrt{T/(rN)}$, nên cần $rN\ge T/\sigma_t^2$ và entropy ≥ $2rN$ mỗi $T_0$):

| σ_t | rN cần | Entropy tối thiểu mỗi $T_0$ | Năng lượng tỏa mỗi $T_0$ (k_BT = 1 MeV) | So với 20–30 MeV |
|---|---|---|---|---|
| 0.071 $T_0$ (99%) | 992 | 1984 $k_B$ | ≥ **1984 MeV** | thiếu **66–99 lần** |
| 0.27 $T_0$ (90%) | 69 | 137 $k_B$ | ≥ 137 MeV | thiếu 5–7 lần |

**Ý nghĩa:** đồng hồ của học thuyết (người quan sát / bể Carnot) **đúng là đồng hồ tự hành kiểu Erker**. Định luật hai nhiệt động chặn độ chính xác của nó: muốn phân giải từng đợt thì phải tỏa năng lượng gấp ~70–100 lần năng lượng hạt nhân có.

## Kết luận chung

1. **Cả năm loại đồng hồ cho cùng một kết quả theo năm con đường khác nhau:** khối lượng (Salecker–Wigner), phản lực và tính tuần hoàn (WKB), va chạm (Brown–Kuchař), bão hòa (entropy hệ), và định luật hai (TUR). Không có đồng hồ **vật lý** nào nằm trong hạt nhân phân giải được từng đợt ở mức 99%.
2. **Đồng hồ tự nhiên của học thuyết là loại 4 + 5:** đồng hồ entropy của bể Carnot, hoạt động như đồng hồ tự hành. Nó **đúng** (đơn điệu, không tuần hoàn, chính là mũi tên thời gian) nhưng **thô**: khoảng 12–21% giới hạn TUR, và cần thêm khoảng 70–100 lần năng lượng mới phân giải được từng đợt.
3. **Đồng hồ lý tưởng τ = $T_0$/25 của Part III** tương đương đồng hồ Brown–Kuchař hình thức. Nó hợp lệ toán học (thỏa ràng buộc WdW) nhưng không có hạt nhân thật nào hiện thực hóa được.

**Câu đề xuất cho bản thảo (tiếng Anh):**

> The natural Page–Wootters clock of the engine is its Carnot bath: an autonomous thermodynamic clock whose ticks are irreversible Carnot records. It is monotonic and aperiodic but operates at 12–21% of the thermodynamic-uncertainty bound; resolving individual bursts (σ_t ≤ 0.07 T₀) would require dissipating ~2 GeV per T₀, two orders of magnitude above the nuclear excitation energy. Salecker–Wigner, WKB and Brown–Kuchař clocks lead to the same conclusion through mass, back-action and collision constraints respectively.
