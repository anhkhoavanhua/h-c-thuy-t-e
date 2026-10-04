# Xung đột Hiện tại ↔ Tương lai trong hạt nhân tán xạ tương lai, và nguồn năng lượng cần thiết

> **Ý của học thuyết:** xung đột giữa Hiện tại và Tương lai trong hạt nhân có thể **tán xạ cả tương lai** đi, nếu có đủ nguồn năng lượng.
> **Nối tiếp:** [xung-dot-hien-tai-tuong-lai-kiem-chung.md](xung-dot-hien-tai-tuong-lai-kiem-chung.md) (mỗi va chạm Carnot là một lần xung đột Hiện tại ↔ Tương lai; N ~ Geom(p)), [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md), [e-nang-luong-thoat-ly](e-nang-luong-thoat-ly-chan-khong-kiem-chung.md), [da-vu-tru-du-lieu-that](da-vu-tru-du-lieu-that-kiem-chung.md) (S_n thật của IAEA).
> **Số liệu:** [giai_so_tan_xa_tuong_lai.py](giai_so_tan_xa_tuong_lai.py), kết quả thô [tan_xa_tuong_lai_ket_qua.txt](tan_xa_tuong_lai_ket_qua.txt), chạy khoảng 10 giây.
> - Dữ liệu: mô phỏng Langevin đầy đủ của Part I (`rl_data/chu_trinh_f.json`, k_BT = 1 MeV, 4 × 120 T0 mỗi f).
> - S_n thật lấy từ IAEA (`rl_data/iaea_ground_states.csv`).

## "Tán xạ tương lai" nghĩa là gì trong mô hình

1. **Tương lai tách nhánh.** Ở lần xung đột thứ k của một đợt, có ba kết cục:
   - Tương lai thắng (tái khóa giữ được);
   - Hiện tại thắng nhưng năng lượng Carnot tích lũy E_exc vẫn dưới S_n;
   - Hiện tại thắng **và** E_exc ≥ S_n.

   Entropy H_k của ba kết cục này là lượng tương lai bị "tán" ra tại lần xung đột đó.
2. **Tán xạ thật trong hạt nhân.** Khi $E_{\text{exc}}=\sum_{i\le k}\Delta E_i\ge S_n$, một neutron **có thể** thoát ra. Phổ năng lượng của nó theo luật Ericson ở Part I: $\varepsilon\sim(E_{\text{exc}}-S_n)\,\mathrm{Beta}(2,n-1)$, với $n=2k\pm1$.
   - **Lưu ý:** mô hình Part I **không** lấy neutron ra khỏi hạt nhân. "Tán xạ" ở đây là **có thể xảy ra về năng lượng**, đúng như cách bài báo phát biểu.
3. **Nguồn năng lượng** là tỉ số lái f của động cơ, tức E năng lượng đổ vào.

S_n thật (IAEA) của các hạt nhân hợp phần: ⁶⁰Ni 11.39, ⁶²Ni 10.60, ⁶⁴Cu 7.92, ⁶⁶Cu 7.07, ⁷⁴Ge 10.20 MeV. Bài báo dùng S_n = 8.0 MeV.

## A. Nguồn năng lượng quyết định tương lai có bị tán xạ không

| f | Số đợt | ⟨N⟩ xung đột/đợt | ⟨E_exc⟩/đợt (MeV) | P(≥S_n) ⁶⁰Ni | ⁶²Ni | ⁶⁴Cu | ⁶⁶Cu | ⁷⁴Ge | bài báo (8) | Tán xạ / T0 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.95 | 16 | 2.12 | 9.0 | 0.19 | 0.25 | 0.50 | 0.50 | 0.31 | 0.50 | 0.017 |
| 1.05 | 69 | 4.51 | 16.8 | 0.45 | 0.46 | 0.62 | 0.67 | 0.51 | 0.62 | 0.093 |
| 1.13 | 75 | 8.61 | 36.2 | 0.69 | 0.72 | 0.80 | 0.81 | 0.75 | 0.80 | 0.130 |

- **Nguồn năng lượng mạnh hơn thì tương lai bị tán xạ nhiều hơn:** số lần tán xạ mỗi T0 tăng **7.6 lần** khi f đi từ 0.95 lên 1.13.
- **Ngưỡng nguồn năng lượng tùy hạt nhân.** Mức f tối thiểu để một nửa số đợt vượt S_n (nội suy tuyến tính):
  - Cu (S_n 7–8 MeV): f ≈ 0.95 là đủ;
  - Ge-74 (10.2 MeV): f ≈ 1.05;
  - Ni-62 (10.6 MeV): f ≈ 1.06;
  - Ni-60 (11.4 MeV): f ≈ 1.07.

  Hạt nhân giữ năng lượng chặt hơn (S_n lớn) thì cần nguồn mạnh hơn mới tán xạ được tương lai của nó.

**Lần xung đột thứ mấy thì vượt S_n = 8 MeV** (trong các đợt có vượt):

| f | Trung vị k | k = 1 | 2 | 3 | 4 | 5 | ≥ 6 |
|---|---|---|---|---|---|---|---|
| 0.95 | 2 | 0.12 | 0.88 | 0 | 0 | 0 | 0 |
| 1.05 | **3** | 0.21 | 0.16 | 0.28 | 0.23 | 0.07 | 0.05 |
| 1.13 | 2 | 0.28 | 0.25 | 0.17 | 0.10 | 0.05 | 0.15 |

Ở f = 1.05, trung vị k = 3, khớp đúng "median step k = 3" của bài báo (V.A).

## B. Tương lai tách nhánh theo từng lần xung đột (f = 1.05, S_n = 8 MeV)

| k | Đợt tới k | P(Tương lai thắng) | P(Hiện tại thắng, < S_n) | **P(Hiện tại thắng, ≥ S_n: tán xạ)** | **H_k (nat)** |
|---|---|---|---|---|---|
| 1 | 69 | 0.35 | 0.61 | 0.04 | 0.81 |
| 2 | 45 | 0.18 | 0.67 | 0.16 | 0.87 |
| **3** | 37 | 0.16 | 0.43 | **0.41** | **1.02** (cực đại) |
| 4 | 31 | 0.26 | 0.16 | 0.58 | 0.96 |
| 5 | 23 | 0.17 | 0.09 | 0.74 | 0.74 |
| 6 | 19 | 0.16 | 0 | 0.84 | 0.44 |
| 7 | 16 | 0.12 | 0 | 0.88 | 0.38 |
| 8 | 14 | 0 | 0 | 1.00 | 0 |

- **Tương lai bị tán xạ mạnh nhất ở lần xung đột thứ 3** (H = 1.02 nat). Đó đúng là lúc năng lượng tích lũy chạm S_n: tương lai chia đều giữa "tái khóa", "tiếp tục dưới ngưỡng" và "vượt ngưỡng phát".
- Sau đó tương lai **hội tụ lại**: khi đợt kéo dài quá 5–6 lần xung đột, gần như chắc chắn đã vượt S_n (H giảm về 0.4 rồi 0). Tán xạ là một **cửa sổ**, không tăng mãi.
- Lưu ý: số đợt giảm dần theo k (69 → 14), nên các giá trị ở k lớn có sai số lớn.

## C. Năng lượng bị tán xạ ra, và giá của nó

**Phổ neutron có thể phát ra** (Ericson, với S_n = 8 MeV):

| f | Số lần tán xạ | ⟨ε⟩ (MeV) | Trung vị | 95% dưới | P(ε > 3 MeV) |
|---|---|---|---|---|---|
| 0.95 | 8 | 1.32 | 1.18 | 2.64 | 0.12 |
| 1.05 | 43 | **2.21** | 1.62 | 5.52 | 0.30 |
| 1.13 | 60 | 2.88 | 2.38 | 7.18 | 0.33 |

Nguồn mạnh hơn thì neutron tán xạ ra mang năng lượng lớn hơn. Ở f = 1.05, ⟨ε⟩ = 2.2 MeV, nằm trên khoảng 0.9–1.6 MeV của bài báo; khác biệt chủ yếu do cách lấy mẫu n = 2k ± 1 ở đây. Phổ đo thật từ EXFOR đã so ở [part2 A9](part2-cau-hoi-mo-A1-A9-kiem-chung.md).

**Giá năng lượng để tán xạ tương lai** (lượng tử, Part III, dùng số đã đo ở [e-khuc](e-khuc-ba-pha-xung-dot-tuong-lai-mo-kiem-chung.md) mục D và [e-nang-luong-thoat-ly](e-nang-luong-thoat-ly-chan-khong-kiem-chung.md)):
- mỗi lần xung đột Hiện tại ↔ Tương lai mở **3.79 nat** tương lai;
- mỗi lần tốn khoảng **27.7 MeV** dòng năng lượng vào;
- tức khoảng **7.3 MeV cho mỗi nat** tương lai bị tán xạ.

Con số này là ước lượng: 27.7 MeV là toàn bộ năng lượng vào chia cho số bản ghi Carnot, nên đã gồm cả phần năng lượng phục vụ xung đột Quá khứ ↔ Hiện tại.

## Kết luận

| Câu hỏi | Trả lời bằng số |
|---|---|
| Xung đột Hiện tại ↔ Tương lai có tán xạ tương lai không? | Có. Mỗi lần xung đột chia tương lai thành ba nhánh, cực đại **1.02 nat ở lần xung đột thứ 3**, đúng lúc năng lượng chạm S_n |
| Có tán xạ thật trong hạt nhân không? | Có khả năng về năng lượng: 50–80% đợt vượt S_n = 8 MeV; phổ neutron ⟨ε⟩ = 1.3–2.9 MeV tùy nguồn |
| Cần bao nhiêu nguồn năng lượng? | Đủ để một nửa số đợt vượt S_n: f ≈ 0.95 (Cu), 1.05 (Ge), 1.06–1.07 (Ni). Giá khoảng 7.3 MeV mỗi nat tương lai |
| Nguồn càng mạnh thì sao? | Tán xạ càng thường xuyên (×7.6 mỗi T0 từ f = 0.95 lên 1.13), neutron ra càng mạnh, nhưng sau khoảng 6 lần xung đột tương lai lại hội tụ (chắc chắn vượt ngưỡng) |

**Câu tiếng Anh cho bản thảo:**

> Each Present–Future contest (a Carnot impact) splits the future into three outcomes: re-locking, continuation below the separation energy, or continuation above it, after which nucleon emission becomes energetically accessible. At f = 1.05 the branching entropy peaks at 1.02 nat at the third contest, exactly where the accumulated Carnot energy reaches S_n (median k = 3, as in the paper), and collapses after about six contests. The energy source decides whether the future scatters: with measured separation energies (IAEA), half of the bursts cross S_n at f ≈ 0.95 for Cu, 1.05 for Ge and 1.06–1.07 for Ni; raising the drive from f = 0.95 to 1.13 multiplies the scattering rate by 7.6 and hardens the emitted spectrum (⟨ε⟩ = 1.3 → 2.9 MeV). In the quantum model each contest opens 3.8 nat of future at ≈ 28 MeV of throughput, about 7 MeV per nat.

## Giới hạn

- **Mô hình cổ điển không lấy neutron ra khỏi hạt nhân.** "Tán xạ" là vượt S_n, tức có thể phát về năng lượng. Phát thật sẽ làm hạt nhân mất năng lượng và đổi động lực sau đó; chưa mô hình hóa.
- E_exc tích lũy toàn bộ năng lượng Carnot của đợt, theo cách bài báo dùng ở V.A. Một phần năng lượng có thể đã rời hạt nhân qua ma sát γx trước khi tích lũy.
- Ở f = 0.95 chỉ có 16 đợt và 8 lần tán xạ, nên sai số lớn.
- Giá 7.3 MeV/nat là phép ghép giữa hai tài liệu (năng lượng mỗi bản ghi và nat mỗi bản ghi). Đây chưa phải phép đo năng lượng chỉ dành riêng cho xung đột Hiện tại ↔ Tương lai.
- S_n thật áp vào động cơ có thang ħω0 = 10 MeV (A ≈ 70). Với Cu, Ni, Ge (A ≈ 60–74), ħω0 thật chỉ lệch 2–5%.

## Câu hỏi còn mở

- Thêm việc **phát thật**: khi E_exc ≥ S_n, cho neutron thoát (trừ năng lượng, kết thúc đợt), rồi xem tương lai sau tán xạ có "tái khóa" về Quá khứ nhanh hơn không.
- Đo giá năng lượng **riêng** cho kênh Carnot (Hiện tại ↔ Tương lai) trong mô hình lượng tử, thay cho phép chia toàn bộ dòng năng lượng vào.
