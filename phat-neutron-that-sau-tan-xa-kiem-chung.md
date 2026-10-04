# Phát neutron thật: tán xạ tương lai là cái van giữ cho động cơ không "cháy"

> **Nối tiếp:** [tan-xa-tuong-lai-hien-tai-tuong-lai-kiem-chung.md](tan-xa-tuong-lai-hien-tai-tuong-lai-kiem-chung.md). Ở bản đó "tán xạ" mới chỉ là vượt S_n (có thể phát về năng lượng). Ở đây neutron **thoát ra thật**.
> **Câu hỏi:** sau khi tán xạ, tương lai có quay về Quá khứ nhanh hơn không? Và việc phát làm gì cho động cơ?
> **Số liệu:** [giai_so_phat_that.py](giai_so_phat_that.py), kết quả thô [phat_that_ket_qua.txt](phat_that_ket_qua.txt). Mỗi nhiệt độ nền chạy khoảng 3 phút trên 8 lõi; mỗi biến thể 8 quỹ đạo × 100 T0, f = 1.05.

## Mô hình: thêm bể kích thích nội và phát thật vào Langevin của Part I

- **Bể kích thích nội** E_res nhận toàn bộ năng lượng Carnot của mọi va chạm.
- **Bể nung nóng hạt nhân** theo khí Fermi: $T_{\text{nuc}}=\sqrt{E_{\text{res}}/a}$, với $a=A/8=8.75$ MeV⁻¹ (A = 70). Nhiễu Langevin dùng $k_BT=\max(k_BT_{\text{nền}},T_{\text{nuc}})$: hạt nhân càng được kích nhiều thì càng "run" mạnh.
- **Phát thật:** khi $E_{\text{res}}\ge S_n=8$ MeV, một neutron thoát với $\varepsilon=(E_{\text{res}}-S_n)\,\mathrm{Beta}(2,n-1)$, $n=2k+1$ (Ericson). Sau đó $E_{\text{res}}\leftarrow E_{\text{res}}-S_n-\varepsilon$.

Ba biến thể:
- **(a)** k_BT cố định, không phản hồi (như Part I);
- **(b)** có nung nóng nhưng **không** phát;
- **(c)** có nung nóng **và** phát thật.

## Kết quả

**Nhiệt độ nền 1 MeV (như bài báo):**

| Biến thể | Gãy / T0 | Quá khứ (% thời gian) | T_nuc TB (MeV) | Phát / T0 | ⟨ε⟩ (MeV) | E_res cuối |
|---|---|---|---|---|---|---|
| (a) chuẩn, k_BT = 1 MeV | 0.940 | 52.4% | 1.00 | — | — | — |
| (b) nung nóng, không phát | **1.619** | 46.9% | **6.18** | — | — | **709 MeV** |
| (c) nung nóng + phát thật | 0.916 | 53.7% | 1.00 | **0.371** | 1.29 | < 8 MeV |

**Nhiệt độ nền 0.5 MeV:**

| Biến thể | Gãy / T0 | Quá khứ | T_nuc TB (MeV) | Phát / T0 | ⟨ε⟩ (MeV) | E_res cuối |
|---|---|---|---|---|---|---|
| (a) chuẩn, k_BT = 0.5 MeV | 0.464 | 63.1% | 0.50 | — | — | — |
| (b) nung nóng, không phát | **1.936** | 38.4% | **6.32** | — | — | **816 MeV** |
| (c) nung nóng + phát thật | 0.529 | 60.3% | 0.64 | 0.224 | 1.29 | < 8 MeV |

**Thời gian từ lúc tán xạ tới khi quay về Quá khứ** (khóa giữ được ≥ T0):

| Nhiệt độ nền | (a) mốc: E_exc vượt S_n, không phát | (c) mốc: neutron thoát thật | Tỉ số (c)/(a) |
|---|---|---|---|
| 1 MeV | TB 1.19 T0, trung vị 0.77 T0 (n = 237) | TB 1.14 T0, trung vị 0.72 T0 (n = 294) | **0.96** |
| 0.5 MeV | TB 0.76 T0, trung vị 0.34 T0 (n = 129) | TB 0.75 T0, trung vị 0.36 T0 (n = 179) | **0.99** |

## Đọc kết quả

1. **Không phát thì hạt nhân "cháy".**
   - E_res dồn tới 709–816 MeV sau 100 T0, và nhiệt độ lên khoảng 6 MeV.
   - Số lần gãy tăng 1.7 lần (nền 1 MeV) và 4.2 lần (nền 0.5 MeV); Quá khứ co lại còn 38–47% thời gian.
   - Hiện tại lấn át, và đi theo hướng **mất kiểm soát**: năng lượng càng dồn thì gãy càng nhiều, gãy càng nhiều thì năng lượng càng dồn.
2. **Phát thật là cái van.** Tán xạ tương lai ra ngoài (neutron mang S_n + ε) giữ E_res dưới S_n. Hạt nhân đứng yên trong chế độ ba pha: Quá khứ 54–60% thời gian, phát 0.22–0.37 lần/T0, ⟨ε⟩ = 1.29 MeV.
3. **Tham số của bài báo tự nhất quán.** Nhiệt độ của hạt nhân ngay tại ngưỡng phát là $\sqrt{S_n/a}=\sqrt{8/8.75}=$ **0.96 MeV**, gần đúng bằng k_BT = 1 MeV mà bài báo chọn. Vì vậy ở nền 1 MeV, (c) gần như trùng (a): bài báo đã ngầm đặt động cơ ở đúng nhiệt độ mà van phát giữ. Ở nền 0.5 MeV, phản hồi nhiệt mới lộ ra: hạt nhân tự nóng lên khoảng 0.64 MeV và gãy nhiều hơn 14%.
4. **Phát không làm tương lai quay về Quá khứ nhanh hơn.** Tỉ số thời gian là 0.96–0.99, tức bằng nhau trong sai số. Đợt đang diễn ra kết thúc theo cơ học va chạm–tái khóa của rô-to, không theo bể nội. Việc phát chỉ **ngăn đợt sau** bị nung quá mức, không rút ngắn đợt hiện tại.

## Kết luận

- Xung đột Hiện tại ↔ Tương lai tán xạ tương lai ra ngoài dưới dạng neutron, và đây **không** phải phần thừa: nó là **van điều nhiệt** của hạt nhân.
- Không có van, năng lượng dồn hàng trăm MeV, Hiện tại thắng mất kiểm soát, và chế độ ba pha (cùng với lời giải frozen time) bị phá.
- Có van, nhiệt độ bị giữ dưới $\sqrt{S_n/a}\approx$ 1 MeV, đúng nhiệt độ của Part I.
- Van **không** làm tương lai quay về Quá khứ nhanh hơn trong một đợt (tỉ số khoảng 1). Thứ nó giữ được là **chu trình dài hạn**.

**Câu tiếng Anh cho bản thảo:**

> Letting the nucleon actually leave closes the energy loop. With Carnot energy feeding an internal Fermi-gas reservoir that heats the nucleus (T = √(E_res/a), a = A/8), but no emission, the reservoir grows to 0.7–0.8 GeV within 100 T₀, the breakdown rate rises 1.7–4.2 times and the Present runs away. With emission at S_n, the reservoir is capped and the temperature stays at √(S_n/a) = 0.96 MeV, which is the k_BT = 1 MeV assumed in Part I, so the paper's bath temperature is self-consistent. Emission does not shorten the ongoing burst (return-to-Past time ratio 0.96–0.99); it acts as a thermostat that preserves the three-phase cycle.

## Giới hạn

- Mô hình nhiệt khí Fermi (T = √(E/a), a = A/8) là gần đúng chuẩn của vật lý hạt nhân, nhưng không phải số đo riêng của từng hạt nhân.
- Bể nội nhận **toàn bộ** năng lượng Carnot. Bể không mất năng lượng qua kênh nào khác ngoài phát neutron: không có phát γ, không có ma sát ra ngoài hạt nhân. Thêm phát γ sẽ hạ E_res, làm (b) bớt cực đoan.
- Phát không làm thay đổi trạng thái cơ học (θ, θ̇, x) ngay lúc phát; nó chỉ tác động qua nhiệt độ.
- Mỗi biến thể chỉ có 8 quỹ đạo × 100 T0.

## Câu hỏi còn mở

- Thêm kênh phát γ (cạnh tranh với neutron dưới S_n) để có một cân bằng năng lượng đầy đủ hơn.
- Cho việc phát **giật lùi** rô-to (bảo toàn động lượng góc khi neutron mang đi một phần ℓ), để xem khi đó tương lai có quay về Quá khứ nhanh hơn không.
