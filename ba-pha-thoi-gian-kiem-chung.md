# Ba pha $\mathcal E_{\text{quá khứ}}\to\mathcal E_{\text{hiện tại}}\to\mathcal E_{\text{tương lai}}$ và các đoạn văn đề xuất cho bài báo: kiểm chứng số liệu

> Tiếp nối [ngan-sach-ericson-kiem-chung.md](ngan-sach-ericson-kiem-chung.md).
> Số liệu tính bằng [giai_so_ba_pha_thoi_gian.py](giai_so_ba_pha_thoi_gian.py) (chạy `python3 giai_so_ba_pha_thoi_gian.py`, ≈ 8 phút trên 8 lõi). Script dùng vòng lặp riêng có **chế độ ngắt theo đợt**, và **không sửa** các file mô phỏng cũ.
> **Bộ thông số:** $f=1.05$, $\tau=10$, $\gamma_x=0.70$, $g=0.1$, $\hbar\omega_0=10$ MeV, $S_n=8$ MeV, $\dot\theta_{\text{gãy}}=3.2048$. Mỗi cấu hình 12 quỹ đạo × $60T_0$.
> **Định nghĩa đợt:** từ lần gãy đến khi một lần khóa kéo dài ≥ $T_0$. Lần khóa đó là điểm bắt đầu của $\mathcal E_{\text{tương lai}}$.
> Mã LaTeX trong tin nhắn đầu **bị cắt**, chỉ còn bảng ánh xạ ba pha, nên phần LaTeX chưa kiểm được.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| **Bảng ba pha.** $\mathcal E_{\text{quá khứ}}$: $n=1$, $\dot\theta\approx0.82\,\dot\theta_{\text{gãy}}$, $\tau_{\text{lưu}}\approx4.5T_0$ | ✓ (0.8208; 4.57 $T_0$). Lưu ý: ở pha này, nguồn lái **giữ** $\dot\theta$ quanh điểm dừng chứ không "tích lũy" dần lên ngưỡng. Gãy xảy ra do nhiễu (xem [recharge-kiem-chung.md](recharge-kiem-chung.md)) |
| **Bảng ba pha.** $\mathcal E_{\text{hiện tại}}$: $n=2$, bội số hình học $p\approx0.195$, ~5 bước | ✓ ($\chi^2=10.9$; TB 5.12) |
| **Bảng ba pha.** $\mathcal E_{\text{tương lai}}$: sau tái khóa $\dot\theta<\dot\theta_{\text{gãy}}$ | **✓ 100%** ở cả 6 cấu hình |
| **§IV.A.** Ngắt $Q_\theta$ trong toàn đợt ⟹ ngân sách đóng, $E_{\text{budget}}\le25$ MeV, **không còn đuôi năng lượng cao** | **Chỉ đúng ở nhiễu 1 MeV.** Ngắt lái, 1 MeV: max **19.9** MeV, **100%** đợt ≤ 25 MeV. Ngắt lái, **5 MeV**: max **271.5** MeV, chỉ **59%** đợt ≤ 25 MeV, vì nhiễu vẫn bơm năng lượng. **Ngắt cả lái lẫn nhiễu** mới đạt 100% ≤ 25 MeV ở cả hai mức (max 11.6 và 20.4 MeV). **Đoạn văn cần nói rõ là nhiễu cũng được tách khỏi hệ trong đợt** |
| $E_{\text{budget}}=\tfrac12I(\dot\theta_{\text{break}}^2-\dot\theta_{\text{quenched}}^2)\le25$ MeV | **Gần đúng nhưng không phải cận chặt:** TB 19.3–24.4 MeV, max 24.3–51.8 MeV; tỉ lệ ≤ 25 MeV là 53–100% tùy cấu hình. Tổng Carnot thực tế **nhỏ hơn nhiều** so với ngân sách (5.9 so với 19.3 MeV khi ngắt cả hai), nên ngân sách chỉ là cận trên lỏng |
| $E^*=S_n+E_{\text{beam}}$ | **Chưa có trong mô hình.** Năng lượng vào hệ đến từ nguồn lái ở $\mathcal E_{\text{quá khứ}}$, không đến từ một đạn bắn. Nếu viết $E^*=S_n+E_{\text{beam}}$ thì cần một cách gán $E_{\text{beam}}$ vào điều kiện đầu của đợt |
| **§IV.B.** Tích lũy tới $E_{\text{exc}}\ge S_n$ thì mở kênh thoát | **Mâu thuẫn với §IV.A.** Khi ngân sách được đóng như §IV.A, **phần lớn đợt không tới $S_n$**: tổng TB 5.9–8.0 MeV < 8 MeV. Tỉ lệ đợt mở kênh: 71% (không ngắt) → **36%** (ngắt lái) → **6.7%** (ngắt lái + nhiễu, 1 MeV). Phát nucleon trở thành biến cố hiếm |
| **§IV.B.** $n=3$, $E_{\max}\approx4.0$ MeV ⟹ $\langle\varepsilon\rangle=2.00$ MeV | **Số học đúng**, nhưng **giả thiết không nhất quán với mô hình**: (1) bước mở kênh $k$ trung vị 2–3, tức $n=2k\pm1\approx3$–7 chứ không cố định 3; (2) $E_{\max}\approx4$ MeV chỉ đúng khi **không** ngắt. Khi ngắt đúng như §IV.A: $E_{\max}$ TB **2.2–3.0 MeV**, $\langle\varepsilon\rangle=$ **0.7–1.5 MeV** (1 MeV nhiễu) |
| "…aligning with experimentally observed pre-equilibrium **evaporative** plateaus **without parameter fitting**" | **Hai vấn đề:** (1) "evaporation" là phát xạ **cân bằng** (từ hạt nhân hợp phần), ngược nghĩa với "pre-equilibrium"; (2) mô hình có ít nhất 8 tham số chọn tay ($f,\tau,\gamma_x,g,E_{\text{vt}}/\hbar\omega_0,\chi,e,k_BT$), và chưa có phổ thực nghiệm nào được so. Nên bỏ cụm "without parameter fitting" và "aligning with experimentally observed" cho tới khi có phép so thật |
| **§III.C.** $0.4\,\hbar\omega_0$ là kỳ vọng thống kê của các chuyển dời rời rạc | **Chấp nhận được về cách diễn giải.** Đo được: xung TB **0.37–0.56 $\hbar\omega_0$**, $P(\text{xung}\ge1\hbar\omega_0)=0$–21%, $P(\ge2\hbar\omega_0)=0$–6.5% |
| **§III.C.** "chuyển dời rời rạc **Δn = +2**" | *(Làm rõ sau khi đọc bản LaTeX đầy đủ: nếu $n$ là **số exciton** $p+h$ thì Δn = +2, tức tạo một cặp hạt–lỗ mỗi bước, là quy ước **đúng** của mô hình exciton. Nhận xét dưới đây chỉ áp dụng khi $n$ là số lượng tử của dao động tử — xem [bai-bao-latex-kiem-chung.md](bai-bao-latex-kiem-chung.md).)* **Nếu $n$ là số lượng tử dao động tử: sai quy tắc chọn và sai dấu.** Liên kết $x=\kappa\sin\theta$ **tuyến tính** theo $x$, nên yếu tố ma trận $\langle n'\lvert x\rvert n\rangle\ne0$ chỉ khi $n'=n\pm1$. Va chạm Carnot là **mất** năng lượng, nên chuyển dời là $\Delta n=-1$ (phát một lượng tử). $\Delta n=\pm2$ chỉ xuất hiện với ghép bậc hai $x^2$. Đoạn văn tiếng Anh đề xuất ở §III.C **không** chứa "Δn = +2" và có thể dùng được |

## 1. Số liệu theo chế độ ngắt

| Chế độ trong đợt | Nhiễu | Số đợt | Tổng Carnot/đợt: TB / trung vị / 90% / max | ≤ 25 MeV | Ngân sách: TB / max | Tổng ≤ ngân sách | $\dot\theta$ cuối < $\dot\theta_{\text{gãy}}$ | Mở kênh | $E_{\max}$ TB | $\langle\varepsilon\rangle$ ($n=2k\pm1$) |
|---|---|---|---|---|---|---|---|---|---|---|
| Không ngắt | 1 MeV | 100 | 19.4 / 15.1 / 41.9 / 73.1 | 75% | 18.2 / 26.3 | 58% | 100% | 71% | 4.23 | 1.22–1.85 |
| Không ngắt | 5 MeV | 108 | 44.8 / 34.7 / 94.0 / 285.9 | 36% | 20.8 / 48.0 | 30% | 100% | 87% | 8.28 | 2.30–3.67 |
| **Ngắt lái** | 1 MeV | 86 | **7.6 / 6.8 / 14.8 / 19.9** | **100%** | 19.7 / 29.2 | 95% | 100% | **36%** | 2.18 | **0.66–1.01** |
| **Ngắt lái** | 5 MeV | 100 | 30.0 / 19.9 / 61.1 / **271.5** | **59%** | 24.4 / 51.8 | 54% | 100% | 76% | 8.20 | 2.36–3.88 |
| **Ngắt lái + nhiễu** | 1 MeV | 104 | **5.9 / 5.6 / 6.6 / 11.6** | **100%** | 19.3 / 24.3 | 94% | 100% | **6.7%** | 3.01 | 1.01–1.53 |
| **Ngắt lái + nhiễu** | 5 MeV | 141 | **8.0 / 5.8 / 14.0 / 20.4** | **100%** | 19.7 / 36.1 | 89% | 100% | **34%** | 2.60 | 0.87–1.31 |

## 2. Hệ quả cho bài báo

1. **§IV.A cần sửa:** tách **cả** $Q_\theta$ lẫn nhiễu Langevin khỏi hệ trong đợt, hoặc giới hạn phát biểu cho $k_BT\lesssim1$ MeV. Khi đó "không còn đuôi năng lượng cao" đứng được: max ≤ 20.4 MeV.
2. **§IV.B phải dùng số của cấu hình §IV.A**, không dùng số của mô hình không ngắt. Ngân sách đóng làm chỉ 7–36% đợt đạt $S_n$, với $\langle\varepsilon\rangle\approx0.7$–1.5 MeV. Nếu giữ "$E_{\max}\approx4$ MeV, $n=3$, $\langle\varepsilon\rangle=2.00$ MeV", phần đó không nhất quán với §IV.A.
3. **Bảng ba pha** có thể giữ nguyên, kèm hai điều chỉnh:
   - $\mathcal E_{\text{quá khứ}}$: đổi "tích lũy năng lượng" thành "duy trì ở điểm dừng; gãy do thăng giáng".
   - $\mathcal E_{\text{tương lai}}$: ghi rõ "với ngân sách đóng, phần lớn đợt kết thúc **không** phát nucleon; năng lượng tiêu tán nội tại".

## Câu hỏi còn mở

- **Mã LaTeX trong tin nhắn bị cắt.** Nếu bạn gửi lại toàn bộ, mình kiểm từng công thức và từng con số trong đó.
- **Chọn một cấu hình chuẩn cho bài.** Với ngân sách đóng (ngắt lái + nhiễu), mô hình dự đoán **phát nucleon hiếm** (7–34%). Đây có phải là hành vi học thuyết mong muốn không? Nếu học thuyết cần phát xạ thường xuyên hơn, phải tăng ngân sách đầu vào (tức $E_{\text{beam}}$), và điều đó cần một điều kiện đầu mới.
