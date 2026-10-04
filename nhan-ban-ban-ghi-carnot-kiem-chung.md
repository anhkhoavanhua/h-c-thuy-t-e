# Tăng độ nhân bản bản ghi Carnot để hai người quan sát luôn đồng ý: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), mục V.A (ánh xạ exciton Ericson). Nối tiếp [ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md) (B2: bất đồng tới 13%) và [part2-cau-hoi-mo-A1-A9-kiem-chung.md](part2-cau-hoi-mo-A1-A9-kiem-chung.md) (A8: R_δ ≈ 1.4–1.7).
> **Số liệu:** [giai_so_nhan_ban_ban_ghi.py](giai_so_nhan_ban_ban_ghi.py), chạy ≈ 3 s.

## Cơ chế trong học thuyết: dòng thác exciton

- Mỗi va chạm Carnot tạo một cặp hạt–lỗ (Δn = +2), đúng ánh xạ Part I.
- Va chạm nucleon–nucleon bên trong tiếp tục tạo cặp mới với tốc độ Williams λ₊ (K = 400 MeV³), cho tới khi đạt $n_{\text{eq}}=\sqrt{1.39gE}$.
- Mỗi nucleon bị kích thích là **một bản sao** của sự kiện "đã có va chạm".
- **Nhiễu nền:** nucleon gần mặt Fermi có thể bị kích thích nhiệt dù không có va chạm, với xác suất q (0.02–0.2 ở k_BT ≈ 1 MeV).

**Số nucleon mang bản sao theo thời gian:**

| Nguồn | $n_{\text{eq}}$ | ⟨n⟩ lúc 0.5 / 1 / 2 / 5 $T_0$ |
|---|---|---|
| Một va chạm (E = 4.1 MeV) | 6 | 2.5 / 2.9 / 3.6 / 4.9 |
| Cả đợt (E = 20 MeV) | 12 | 4.0 / 5.4 / 7.4 / 10.8 |

Năng lượng giới hạn số bản sao: $n\le n_{\text{eq}}\propto\sqrt E$. Một va chạm chỉ đủ năng lượng tạo khoảng 6 bản sao.

## Độ nhân bản đạt được

$R_\delta=N_{\text{act}}/m_\delta$, với $m_\delta$ là số nucleon tối thiểu cần đọc để biết ≥ 0.9 bit (N_act là số nucleon gần mặt Fermi mà môi trường đọc được).

| Trường hợp | q | t ($T_0$) | Thông tin trên 1 nucleon (bit) | $R_\delta$ |
|---|---|---|---|---|
| Một va chạm, N_act = 30 | 0.02–0.2 | 0.5–5 | 0.004–0.06 | **không đạt 0.9 bit** |
| Cả đợt, N_act = 30 | 0.02 | 2 / 5 | 0.10 / 0.17 | 1.4 / **2.5** |
| Cả đợt, N_act = 30 | 0.1 | 5 | 0.10 | 1.4 |
| Cả đợt, N_act = 30 | 0.2 | 5 | 0.07 | không đạt |
| Cả đợt, **N_act = 12** (dải hẹp) | 0.02 / 0.1 | 5 | 0.69 / 0.55 | **4.0** / 2.4 |
| Cả đợt, N_act = 60 | 0.02 / 0.1 | 5 | 0.07 / 0.03 | 1.9 / không đạt |
| *Cơ chế hiện tại (A8)* | | | | *1.4–1.7* |

**Hai người quan sát, mỗi người đọc m nucleon riêng** (cả đợt, t = 2 $T_0$; ô ghi P(cùng kết luận) / P(cả hai đúng)):

| q | m = 1 | m = 3 | m = 5 | m = 10 |
|---|---|---|---|---|
| 0.02 | 78% / 50% | 70% / 62% | 74% / 72% | **80% / 78%** |
| 0.2 | 60% / 41% | 58% / 43% | 57% / 49% | 65% / 58% |

## Kết luận

1. **Dòng thác exciton có tăng độ nhân bản, nhưng ít:** cao nhất R_δ ≈ 2.5–4, chỉ khi cả đợt va chạm (≈ 20 MeV), nhiễu nhiệt thấp (q = 0.02), và dải nucleon hoạt động hẹp. Một va chạm đơn lẻ (≈ 4 MeV) thì **không** đạt được tính khách quan.
2. **Hai giới hạn vật lý ngăn việc tăng thêm:**
   - **Năng lượng:** số bản sao ≤ $n_{\text{eq}}\approx\sqrt{1.39gE}$, tức 6 cho một va chạm và 12 cho cả đợt.
   - **Pha loãng và nhiễu:** bản sao rải ngẫu nhiên trong vài chục nucleon, trộn với kích thích nhiệt. Mỗi nucleon chỉ mang 0.004–0.7 bit.
3. **Bên trong hạt nhân, hai người quan sát không thể luôn đồng ý.** Ngay cả trường hợp tốt nhất cũng chỉ đạt khoảng 80% cùng kết luận.

## Cơ chế đạt được đồng ý: bản ghi ra ngoài hạt nhân

Một **nucleon được phát ra** ngoài hạt nhân là bản ghi bền: nó không bị hạt nhân hấp thụ lại. Nếu nó tương tác với môi trường vĩ mô (như trong máy dò), nó được sao chép vào rất nhiều nguyên tử. Phần sao chép vĩ mô là **giả thiết** về môi trường bên ngoài, không tính ở đây. Lấy xác suất một đợt có phát nucleon từ mục A9 (85–90%):

| Xác suất đợt phát nucleon | Tỉ lệ thời gian có đợt | Hai người cùng kết luận | Cả hai đúng |
|---|---|---|---|
| 0.852 | 0.1 | **100%** | 98.5% |
| 0.852 | 0.5 | **100%** | 92.6% |
| 0.896 | 0.1 | **100%** | 99.0% |
| 0.896 | 0.5 | **100%** | 94.8% |

- Đọc cùng một bản ghi vĩ mô thì hai người **luôn đồng ý**.
- Nhưng họ **cùng sai** ở 10–15% số đợt không phát nucleon: đợt đó tiêu tán bên trong bằng gamma hoặc exciton, không để lại dấu vết bên ngoài.
- **Đề xuất cho học thuyết:** tính khách quan của "Hiện tại" không đến từ bên trong hạt nhân. Nó đến từ **kênh phát xạ**: nucleon (85–90% số đợt) hoặc gamma đi ra ngoài và được khuếch đại. Người quan sát mà mọi người cùng đồng ý là **môi trường bên ngoài hạt nhân**, không phải các nucleon bên trong.

## Câu hỏi còn mở
- Muốn 100% số đợt được ghi bền, cần tính cả **gamma** phát ra ở các đợt không phát nucleon (10–15%). Mô hình phát xạ ở mục R2 đã có kênh γ (GDR Brink–Axel). Tính xác suất mỗi đợt để lại ít nhất một lượng tử bay ra ngoài (nucleon hoặc gamma) là bước tiếp theo.
