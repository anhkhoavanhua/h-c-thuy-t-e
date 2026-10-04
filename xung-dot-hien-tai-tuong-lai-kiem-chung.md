# Xung đột Hiện tại ↔ Tương lai: phần còn lại của lời giải "frozen time"

> **Nối tiếp:** [xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md). Bài đó cho thấy xung đột Quá khứ ↔ Hiện tại tạo ra dòng bản ghi J > 0.
> **Cơ sở:** HocThuyetE (Part I: gãy, đợt Carnot $N\sim\mathrm{Geom}(p)$, tái khóa $\mathcal E_{\text{future}}$), cùng các tài liệu trong `hoc-thuyet-e/`.
> **Số liệu:** [giai_so_xung_dot_hien_tai_tuong_lai.py](giai_so_xung_dot_hien_tai_tuong_lai.py), kết quả thô [xung_dot_hien_tai_tuong_lai_ket_qua.txt](xung_dot_hien_tai_tuong_lai_ket_qua.txt), chạy khoảng 1.3 phút trên 8 lõi. Mô phỏng Langevin đầy đủ của Part I (γx = 0.70); mỗi trường hợp 4 quỹ đạo × 120 T0.

## Xung đột Hiện tại ↔ Tương lai là gì trong mô hình

- Sau khi gãy, hạt bay tự do (**Hiện tại**). Mỗi lần nó cắt lại đa tạp, một va chạm dẻo Carnot buộc nó khóa lại. Đó là **Tương lai thử thắng**.
  - Nếu khóa giữ được (θ̇ đã dưới ngưỡng, ≥ T0 không gãy lại), **Tương lai thắng** và đợt kết thúc.
  - Nếu phản lực lập tức vượt R lần nữa, **Hiện tại thắng lại** và hạt lại bay tự do.
- **Mỗi bản ghi Carnot là đúng một lần xung đột Hiện tại ↔ Tương lai.** Đây là cấu trúc của mô hình: năng lượng Carnot chỉ được ghi tại va chạm tái khóa. Số bản ghi trong một đợt, N, là số lần hai bên giằng co. Luật $N\sim\mathrm{Geom}(p)$ nghĩa là Tương lai thắng với xác suất p ở mỗi lần.
- Hệ quả trực tiếp: **nếu cấm tái khóa (chỉ có Hiện tại), sẽ không có va chạm nào, nên 0 bản ghi Carnot và mũi tên dừng.** Thiếu Tương lai cũng đóng băng thời gian, giống như thiếu Quá khứ ở bài trước.

## 1. Cán cân xung đột

Mặc định k_BT = 1 MeV, lái bật cả khi hạt tự do.

| Trường hợp | Số đợt | ⟨N⟩ (xung đột/đợt) | **p = P(Tương lai thắng)** | χ²(Geom) / bậc tự do | Hiện tại chiếm | Bản ghi Carnot / T0 |
|---|---|---|---|---|---|---|
| f = 0.95 (lái yếu) | 16 | 2.12 | **0.47** | 1.1 / 0 | 1.5% | 0.09 |
| f = 1.05 (bài báo) | 69 | 4.51 | **0.22** | 9.0 / 6 | 17.4% | 0.69 |
| f = 1.13 (lái mạnh) | 75 | 8.61 | **0.12** | 10.6 / 8 | 35.7% | 1.55 |
| f = 1.05, tắt lái khi tự do | 65 | 3.80 | 0.26 | **11.9 / 4** | 14.7% | 0.59 |
| f = 1.05, k_BT = 0.5 MeV | 54 | 2.85 | 0.35 | 2.5 / 3 | 10.1% | 0.38 |

- **Cán cân dịch đúng chiều.** Lái mạnh hơn (Hiện tại mạnh hơn) thì Tương lai thắng ít hơn (p giảm từ 0.47 xuống 0.12). Hiện tại chiếm nhiều thời gian hơn (1.5% → 36%), và bản ghi sinh ra nhanh hơn 17 lần (0.09 → 1.55/T0).
- Tắt lái khi hạt tự do, hoặc giảm nhiệt độ, nghiêng cán cân về Tương lai (p tăng lên 0.26 và 0.35).
- **Nhịp của mũi tên thời gian do xung đột Hiện tại ↔ Tương lai quyết định.** Mỗi lần giằng co cho ra một bản ghi, và Hiện tại càng mạnh thì càng nhiều lần giằng co trước khi Tương lai thắng.
- Luật Geom khớp tốt (χ² ≲ bậc tự do) ở 4/5 trường hợp. Riêng "tắt lái khi tự do" thì lệch (11.9/4).

## 2. Xung đột có "không nhớ" không? (hàm nguy cơ)

$h_k$ là xác suất Tương lai thắng ở va chạm thứ k, với điều kiện đợt đã tới va chạm k.

| Trường hợp | k = 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| f = 1.05 | **0.35** | 0.18 | 0.16 | 0.26 | 0.17 | 0.16 | 0.12 |
| f = 1.13 | 0.16 | 0.19 | 0.10 | 0.17 | 0.03 | 0.14 | 0.16 |
| k_BT = 0.5 | 0.39 | 0.42 | 0.26 | 0.21 | 0.27 | 0.38 | — |

- Từ va chạm thứ 2 trở đi, h_k gần như hằng số: xung đột **không nhớ** và đợt là Geom.
- Ở f = 1.05, **lần xung đột đầu tiên khác**: Tương lai thắng 35%, cao gấp đôi các lần sau (khoảng 17%). Phân bố N vì thế là **hai thành phần**: một phần đợt chỉ có một bản ghi (Hiện tại co thành một khoảnh khắc), phần còn lại là chuỗi Geom.
  - Bản lượng tử cho cùng hiện tượng nhưng mạnh hơn: **40%** đợt chỉ có một bản ghi, ⟨N⟩ = 1.51, tương ứng p ≈ 0.66 ([ba-pha-xung-dot](ba-pha-xung-dot-dong-ho-nguoi-quan-sat-kiem-chung.md)).
- Ở f = 1.13 không thấy hiệu ứng lần đầu này (0.16 so với 0.10–0.19).

## 3. Tương lai thắng thì trí nhớ bị xóa (chu trình tái lập)

Tương quan giữa hai chu trình liên tiếp:

| Trường hợp | corr(N_i, N_{i+1}) | corr(lưu_i, lưu_{i+1}) | corr(N_i, lưu_{i+1}) | Số cặp |
|---|---|---|---|---|
| f = 1.05 | +0.09 | −0.06 | −0.07 | 65 |
| f = 1.13 | −0.08 | +0.01 | +0.16 | 71 |
| tắt lái khi tự do | −0.15 | +0.08 | −0.09 | 61 |
| k_BT = 0.5 | +0.02 | +0.03 | +0.27 | 50 |

Sai số chuẩn khoảng ±0.13. Mọi tương quan đều nằm trong khoảng 2σ, nên **các chu trình gần như độc lập**. Khi Tương lai thắng, hệ quay về θ̇ ≈ 2.15 bất kể đợt trước dài hay ngắn, và quá trình trở thành **quá trình tái lập**. Mỗi chu trình là một "bây giờ" mới, không mang theo ký ức của chu trình cũ, chỉ để lại bản ghi trong bể. Riêng f = 0.95 chỉ có 12 cặp nên không dùng được.

## Bức tranh đầy đủ: hai xung đột cùng phá đóng băng

| Xung đột | Đại lượng sinh ra | Nếu một bên thắng tuyệt đối | Số đo |
|---|---|---|---|
| **Quá khứ ↔ Hiện tại** (lái lên ngưỡng ↔ giữ khóa) | **Khởi đầu** mỗi "bây giờ": tần suất gãy | Quá khứ thắng: không gãy (k_BT = 0: 0 lần; lượng tử: J = 4·10⁻¹⁵). Hiện tại thắng: kẹt trên ngưỡng, J = 0 | gãy 0.09–1.55 lần/T0 |
| **Hiện tại ↔ Tương lai** (gãy lại ↔ tái khóa Carnot) | **Viết bản ghi**: mỗi lần giằng co là 1 bản ghi; **kết thúc** "bây giờ" và xóa trí nhớ | Tương lai bị cấm: 0 va chạm, 0 bản ghi, mũi tên dừng. Tương lai luôn thắng ngay (p = 1): mỗi đợt 1 bản ghi, Hiện tại thành một khoảnh khắc | p = 0.12–0.47; ⟨N⟩ = 2.1–8.6 |

Đọc chung:
- Trong vũ trụ WdW tĩnh, **xung đột thứ nhất mở** mỗi khoảnh khắc Hiện tại.
- **Xung đột thứ hai lấp đầy** khoảnh khắc đó bằng bản ghi (mũi tên), rồi **đóng** nó lại, xóa trí nhớ để chu trình sau bắt đầu mới.
- Thời gian quan hệ có hướng cần **cả hai**.
- **Tốc độ** của mũi tên (bản ghi/T0) tăng khi Hiện tại mạnh lên trong xung đột thứ hai. **Độ hạt** của thời gian (số bản ghi mỗi "bây giờ", ⟨N⟩ = 1/p) do p quyết định.

**Câu tiếng Anh cho bản thảo:**

> Each Carnot record is one contest between Present and Future: the plastic re-locking impact either holds (Future wins, probability p) or the constraint breaks again (Present wins). The burst size N ~ Geom(p) is the statistics of these contests; the hazard is constant after the first impact (memoryless), with an enhanced first-impact win (0.35 vs 0.17 at f = 1.05). Strengthening the drive shifts the balance toward the Present (p = 0.47 → 0.12, Present time fraction 1.5% → 36%) and accelerates the record current seventeenfold (0.09 → 1.55 per T₀). When the Future wins, successive cycles are uncorrelated (|r| ≲ 0.16), so each cycle is a fresh renewal. The Past–Present conflict opens each moment; the Present–Future conflict writes its records and closes it. Forbidding re-locking yields no Carnot records at all.

## Giới hạn

- "Mỗi bản ghi là một lần xung đột" và "cấm tái khóa thì 0 bản ghi" là **cấu trúc của mô hình Part I** (năng lượng Carnot chỉ ghi tại va chạm tái khóa), không phải kết quả đo độc lập.
- Mỗi trường hợp chỉ có 50–75 đợt (f = 0.95 chỉ 16), nên h_k ở k lớn và các tương quan có sai số khoảng ±0.13.
- Ở f = 1.05 lần này đo được p = 0.22, ⟨N⟩ = 4.51, so với 0.195 và 5.12 của bài báo. Chênh lệch nằm trong mức thống kê của 69 đợt (≈ ±12%).
- Phân loại pha dùng quy tắc của `giai_so_burst_g.analyse`: đợt kết thúc ở lần khóa đầu tiên giữ được ≥ T0.

## Câu hỏi còn mở

- Vì sao lần xung đột đầu tiên dễ thắng hơn ở f = 1.05 mà không ở f = 1.13? Một giả thuyết: lần gãy đầu thường xảy ra sát ngưỡng, nên va chạm đầu tiên đủ kéo θ̇ xuống dưới. Có thể kiểm bằng phân bố θ̇ ngay trước va chạm đầu tiên.
- Nối hai bài thành một trạng thái WdW hữu hạn (rô-to + hạt + thanh ghi), để kiểm cả hai xung đột trong cùng một |Ψ⟩⟩.
