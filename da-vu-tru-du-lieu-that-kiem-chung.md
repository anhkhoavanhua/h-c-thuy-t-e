# Độ lệch thời gian chuẩn với dữ liệu thật: "đa hạt nhân" thay cho đa vũ trụ

> **Điều phải nói trước:** **không tồn tại** dữ liệu thực nghiệm về E năng lượng giữa các vũ trụ. Chưa có quan sát nào chạm tới một vũ trụ khác. Vì vậy ở đây tôi dùng dữ liệu thật **trong miền của chính Học thuyết E** (hạt nhân, A ≈ 70, ħω0 ≈ 10 MeV): **mỗi hạt nhân, hoặc mỗi thí nghiệm, là một động cơ E_dual**. Con số thu được là **σ_t giữa các hạt nhân thật**, không phải giữa các vũ trụ.
> **Nối tiếp:** [da-vu-tru-do-lech-thoi-gian-chuan-kiem-chung.md](da-vu-tru-do-lech-thoi-gian-chuan-kiem-chung.md), nơi con số mô hình là σ_t = 1.22 ± 0.06 (cửa sổ vận hành của bài báo).
> **Số liệu:**
> - [giai_so_da_vu_tru_du_lieu_that.py](giai_so_da_vu_tru_du_lieu_that.py), kết quả thô [da_vu_tru_du_lieu_that_ket_qua.txt](da_vu_tru_du_lieu_that_ket_qua.txt), chạy khoảng 30 giây.
> - Dữ liệu: `rl_data/iaea_ground_states.csv`, tải ngày 2026-10-02 từ IAEA LiveChart, `nds.iaea.org/relnsd/v1/data?fields=ground_states&nuclides=all`, 3388 hạt nhân.
> - Dữ liệu: `rl_data/exfor_data/`, 7 bộ EXFOR đã dùng ở [part2 A9](part2-cau-hoi-mo-A1-A9-kiem-chung.md).

## Ánh xạ dữ liệu thật vào động cơ

| Bước | Công thức | Nguồn |
|---|---|---|
| Thang thời gian riêng của mỗi hạt nhân | ħω0(A) = 41·A^(−1/3) MeV, T0(A) = 2πħ/ħω0(A) | hệ thống học vỏ hạt nhân |
| E năng lượng của động cơ (đơn vị ħω0) | **IAEA:** ε = S_n/ħω0(A), năng lượng hạt nhân giữ được trước khi neutron thoát. **EXFOR:** ε = E*/ħω0(A), với E* = E_tới + S_tách | S_n: IAEA; E_tới: EXFOR |
| ε → trạng thái động cơ | Tìm cường độ lái a sao cho ΔE(a) = ε, theo đường cong **chính xác** của động cơ Part III | mô hình |
| Nhịp thời gian | r(a) bản ghi/T0, nhịp vật lý = r(a)/T0(A) | mô hình |
| Độ lệch chuẩn | σ_t = std ln(r/⟨r⟩), chỉ là tỉ số giữa các đồng hồ, trạng thái chung vẫn tĩnh (WdW) | |

Đường cong của động cơ: khi a đi từ 10⁻³ tới 10³, ΔE đi từ 0.20 tới 8.53 ħω0 và r đi từ 0.0007 tới 1.079 /T0.

## I. IAEA: 286 hạt nhân có mặt trong tự nhiên

ħω0(A) = 6.6–32.5 MeV; S_n = 1.66–20.58 MeV; ε = 0.07–1.40. Có 2 hạt nhân nằm dưới dải của mô hình và được kẹp ở mép dưới.

| Tập hợp | Số hạt nhân | **σ_t** | Hệ số lệch điển hình | r_max / r_min | Phần do thang T0(A) | Phần do trạng thái động cơ r(ε) |
|---|---|---|---|---|---|---|
| Mọi hạt nhân bền | 286 | 0.703 | ×2.02 | 90.5 | 0.258 | 0.834 |
| **A = 40–100 (miền của học thuyết)** | **85** | **0.540** | **×1.72** | **10.7** | 0.091 | 0.558 |
| A ≥ 20 | 271 | 0.511 | ×1.67 | 30.8 | 0.189 | 0.543 |

- **Đồng hồ nhanh nhất:** Si-28 (S_n = 17.2 MeV, ε = 1.27). **Chậm nhất:** Be-9 (S_n = 1.66 MeV, ε = 0.08, nằm dưới dải mô hình nên đã bị kẹp).
- **Phần lớn độ lệch đến từ E năng lượng (S_n), không đến từ kích thước hạt nhân.** Riêng thang T0(A) chỉ cho σ = 0.09 trong miền A = 40–100; riêng trạng thái động cơ r(ε) cho 0.56. Đúng như học thuyết nói: tốc độ thời gian phụ thuộc chủ yếu vào E năng lượng ở đó.
- Hai phần không cộng thẳng: phần thang và phần trạng thái tương quan một phần với nhau.

## II. EXFOR: 7 thí nghiệm phản ứng hạt nhân thật

| Thí nghiệm | A | E_tới (MeV) | E* (MeV) | ε = E*/ħω0 | r(ε) (/T0) | Nhịp vật lý (/s) |
|---|---|---|---|---|---|---|
| Cu(n,xn) 14.1 — Baba 1988 | 64.5 | 14.1 | 21.7 | 2.12 | 0.204 | 5.04·10²⁰ |
| Cu(n,xn) 18.0 — Baba 1988 | 64.5 | 18.0 | 25.6 | 2.50 | 0.256 | 6.32·10²⁰ |
| Cu(n,xn) 14.1 — Takahashi 1988 | 64.5 | 14.1 | 21.7 | 2.12 | 0.204 | 5.04·10²⁰ |
| Ge(n,xn) 14.1 — Gotoh 1992 | 73.6 | 14.1 | 21.5 | 2.20 | 0.214 | 5.06·10²⁰ |
| Ni(n,xn) 14.1 — Yabuta 1986 | 59.7 | 14.1 | 23.1 | 2.20 | 0.214 | 5.44·10²⁰ |
| Ni(n,xn) 18.0 — Yabuta 1986 | 59.7 | 18.0 | 27.0 | 2.57 | 0.265 | 6.72·10²⁰ |
| ⁶⁰Ni(p,xn) 22 — Biryukov 1980 | 61.0 | 22.0 | 26.8 | 2.57 | 0.265 | 6.67·10²⁰ |

**σ_t (7 thí nghiệm) = 0.125** (lệch điển hình ×1.13), r_max/r_min = 1.33. Các thí nghiệm kích hạt nhân ở năng lượng gần nhau (E* ≈ 22–27 MeV), nên đồng hồ của chúng gần như cùng nhịp. Năng lượng tới cao hơn (18–22 MeV so với 14 MeV) thì nhịp nhanh hơn khoảng 25–30%.

## So sánh và con số chuẩn

| Nguồn | σ_t | Lệch điển hình | Lệch cực đại |
|---|---|---|---|
| Mô hình, cửa sổ vận hành f ∈ (0.90, 1.139) của bài báo | 1.22 ± 0.06 | ×3.4 | ×57 |
| **Dữ liệu thật IAEA, A = 40–100** | **0.54** | **×1.72** | **×10.7** |
| Dữ liệu thật IAEA, mọi hạt nhân bền | 0.70 | ×2.02 | ×90 |
| Dữ liệu thật EXFOR, 7 thí nghiệm 14–22 MeV | 0.125 | ×1.13 | ×1.33 |

**Con số chuẩn dựa trên dữ liệu thật: σ_t = 0.54**, tức các đồng hồ lệch nhau điển hình **×1.7**. Con số này tính cho 85 hạt nhân bền trong miền A = 40–100 của học thuyết, với E năng lượng của mỗi hạt nhân là S_n/ħω0 đo thật của IAEA. Trên toàn bảng tuần hoàn, σ_t = 0.70.

Con số mô hình trước đây (1.22) **lớn hơn** con số từ dữ liệu thật hơn hai lần. Cửa sổ lái của bài báo trải động cơ qua vùng dốc nhất của đường r(f), trong khi E năng lượng thật của các hạt nhân bền phân bố hẹp hơn trong không gian nhịp.

**Câu tiếng Anh cho bản thảo:**

> No experimental data exist on energies between universes; within the theory's own domain we instead treat each nucleus as an E_dual engine. Using measured neutron separation energies of the 85 naturally occurring nuclei with A = 40–100 (IAEA), mapped to the engine's drive through its exact energy curve and to physical time through ħω₀(A) = 41 A^{-1/3} MeV, the relational clock dispersion is σ_t = 0.54 (typical ratio 1.7, extreme 10.7); over all 286 naturally occurring nuclei it is 0.70. Most of the dispersion comes from the stored energy S_n/ħω₀, not from the nuclear size scale (0.56 versus 0.09). Seven EXFOR neutron-emission experiments at 14–22 MeV give σ_t = 0.125.

## Giới hạn

- **Đây không phải dữ liệu đa vũ trụ.** Đó là dữ liệu hạt nhân thật, được đọc qua lăng kính của học thuyết.
- **Ánh xạ là lựa chọn của học thuyết, không phải định luật:**
  - "E năng lượng = S_n/ħω0" (IAEA) và "= E*/ħω0" (EXFOR) là hai lựa chọn khác nhau, nên hai con số không so sánh trực tiếp được với nhau.
  - ħω0 = 41·A^(−1/3) MeV là công thức hệ thống học, không phải số đo riêng cho từng hạt nhân.
- **Giá trị S_n/ħω0 (0.07–1.40) nằm ở đoạn dốc đầu của đường ΔE(a)** (0.20–8.53), nên σ_t nhạy với cách ánh xạ. Có 2 hạt nhân rơi ra ngoài dải và bị kẹp.
- Động cơ dùng là mô hình Part III 36 chiều. Cùng ánh xạ đó áp vào một động cơ lớn hơn có thể cho đường r(ε) khác.

## Câu hỏi còn mở

- Có muốn thêm dữ liệu mật độ mức (level density) thật, ví dụ RIPL-3, để thay ħω0(A) bằng một thang năng lượng đo riêng cho từng hạt nhân không?
- Thử ánh xạ khác cho E năng lượng, ví dụ năng lượng liên kết trên nucleon hoặc năng lượng kích thích của mức 2⁺ đầu tiên, để xem σ_t nhạy với lựa chọn này đến đâu.
