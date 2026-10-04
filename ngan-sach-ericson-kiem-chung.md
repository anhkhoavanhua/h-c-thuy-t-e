# Ngân sách năng lượng, phổ Ericson và khung "mô hình bán cổ điển hiệu dụng": giải và kiểm chứng số liệu

> Tiếp nối [exciton-kiem-chung.md](exciton-kiem-chung.md).
> Số liệu tính bằng [giai_so_ngan_sach_ericson.py](giai_so_ngan_sach_ericson.py) (chạy `python3 giai_so_ngan_sach_ericson.py`, ≈ 5.5 phút trên 8 lõi).
> [giai_so_dong_co_langevin.py](giai_so_dong_co_langevin.py) được thêm tùy chọn `drive_free=False`, tắt mô-men lái khi hạt ở pha tự do (nhiễu ξθ vẫn giữ). Mặc định giữ nguyên hành vi cũ.
> **Bộ thông số:** $f=1.05$, $\tau=10$, $\gamma_x=0.70$, $g=0.1$, $\hbar\omega_0=10$ MeV, $S_n=8$ MeV. Mỗi cấu hình 16 quỹ đạo × $60T_0$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| Ngân sách $E_{\text{budget}}\approx\tfrac12I(\dot\theta_{\text{gãy}}^2-\dot\theta_{\text{cuối}}^2)\approx15$–$25$ MeV | **✓ Đúng cỡ.** Đo được: TB **19.0–20.5 MeV**, 10–90% ≈ 12–26 MeV (nhiễu 1 MeV) |
| Cắt lái trong pha tự do ⟹ chuỗi dừng khi cạn ngân sách, **loại đuôi trăm MeV** | **Chỉ đúng một phần ở 1 MeV, sai ở 5 MeV.** Ở 1 MeV: tổng mỗi đợt TB 24.7 → **17.4** MeV, max 134.6 → **80.6** MeV. Đuôi giảm nhưng **vẫn vượt ngân sách** (80.6 so với ngân sách 90% ≈ 26 MeV). Ở 5 MeV: **không giảm** (TB 50.9 → 62.3, max 216.6 → **321.5** MeV); tổng/ngân sách TB **2.5**. Lý do: nhiễu $\xi_\theta,\xi_x$ và nguồn lái trong các pha khóa ngắn giữa đợt vẫn tiếp tục bơm năng lượng |
| Mật độ Ericson $\rho=g^nE^{n-1}/(p!h!(n-1)!)$; $P_n(\varepsilon)\propto\varepsilon(1-\varepsilon/E_{\max})^{n-2}$ | **✓ Đúng dạng** (với $\sigma_{\text{inv}}$ hằng số) |
| $\langle\varepsilon\rangle=E_{\max}/(n+1)$ | **Sai hệ số 2.** $P(\varepsilon)$ là phân bố **Beta(2, n−1)**, nên **$\langle\varepsilon\rangle=2E_{\max}/(n+1)$**. Tích phân số khớp: $E_{\max}=4$ cho **2.00 / 1.60 / 1.33 MeV** với $n=3/4/5$, trong khi đề xuất ghi 1.00 / 0.80 / 0.67 |
| $E_{\text{exc}}(3)\approx11$–$13$ MeV | **Hơi cao.** Đo được TB **10.4–11.4**, trung vị 8.5–10.5 MeV |
| $E_{\text{exc}}-S_n\approx4$ MeV tại bước mở kênh | **✓ ở 1 MeV:** TB **3.2–3.8**, trung vị 2.8 MeV. Ở 5 MeV: TB **7.9–8.9** MeV |
| $\langle\varepsilon\rangle\approx0.8$–$1.0$ (bảng cuối ghi 0.8–1.5) MeV | **Con số cuối gần đúng ở 1 MeV, nhưng do hai sai số triệt tiêu nhau** (thiếu hệ số 2, và $n$ quá nhỏ). Tính đúng trên từng đợt mô phỏng (với $n=2k\pm1$): **0.9–1.6 MeV** ở 1 MeV; **2.4–4.1 MeV** ở 5 MeV. Vẫn phụ thuộc nhiệt độ môi trường |
| "$\langle\varepsilon\rangle<E_{\text{vt}}$ nên triệt tiêu nghịch lý 100 MeV" | **Không triệt tiêu.** $\langle\varepsilon\rangle$ nhỏ là **trung bình**. Tổng năng lượng các đợt (tức năng lượng kích thích) vẫn có đuôi tới 80–320 MeV. Ngoài ra, động năng của nucleon **đã phát ra** được đo bên ngoài hạt nhân, nên so với độ sâu hố thế $E_{\text{vt}}$ không phải là tiêu chí đúng |
| Bảng cuối: phi Hamilton, $f\in(0.90,1.139)$ | ✓ Đã cập nhật theo kiểm chứng trước |
| Đóng khung là "Effective Classical Toy Model / Semiclassical Analog" | **Đồng ý, và nên ghi thêm các số liệu hợp lệ sau** (mục 3). Các số này cho thấy mô hình nằm ở giới hạn **lượng tử mạnh**, nên chỉ nên xem là mô hình tương tự định tính |

## 1. Phổ Ericson

$$P_n(\varepsilon)\propto\varepsilon\,(E_{\max}-\varepsilon)^{n-2}\ \Longrightarrow\ \frac{\varepsilon}{E_{\max}}\sim\text{Beta}(2,\,n-1),\qquad \langle\varepsilon\rangle=\frac{2E_{\max}}{n+1}$$

| $E_{\max}$ | $n$ | Tích phân số | $2E_{\max}/(n+1)$ | Đề xuất $E_{\max}/(n+1)$ |
|---|---|---|---|---|
| 4 | 3 | 2.0000 | 2.0000 | 1.0000 |
| 4 | 4 | 1.6000 | 1.6000 | 0.8000 |
| 4 | 5 | 1.3333 | 1.3333 | 0.6667 |

## 2. Ngân sách năng lượng: lái bật và tắt trong pha tự do

| Lái ở pha tự do | Nhiễu | Số đợt | Tổng Carnot/đợt: TB / trung vị / 10–90% / max | Ngân sách quay: TB / 10–90% | Tổng/ngân sách | $E_{\text{exc}}(3)$ TB | Mở kênh | $E_{\text{exc}}-S_n$ tại bước mở: TB | $\langle\varepsilon\rangle$ Ericson đúng |
|---|---|---|---|---|---|---|---|---|---|
| Bật | 1 MeV | 108 | 24.7 / 17.1 / 5.1–53.3 / **134.6** | 19.0 / 12.2–25.0 | 1.30 | 10.87 | 90/108 | 3.78 | 1.05–1.56 |
| Bật | 5 MeV | 119 | 50.9 / 43.2 / 8.8–116.4 / **216.6** | 20.6 / 3.6–35.3 | 2.47 | 11.25 | 110/119 | 8.90 | 2.37–3.82 |
| **Tắt** | 1 MeV | 127 | **17.4** / 11.3 / 4.2–40.4 / **80.6** | 20.5 / 14.4–26.5 | **0.85** | 10.38 | 82/127 | 3.23 | **0.89–1.33** |
| **Tắt** | 5 MeV | 114 | **62.3** / 38.1 / 6.9–173.4 / **321.5** | 24.5 / 8.3–36.7 | **2.54** | 11.37 | 96/114 | 7.92 | 2.49–4.14 |

**Đọc bảng:**
- Tắt lái ở pha tự do chỉ có tác dụng khi nhiễu nhỏ: trung bình giảm khoảng 30%, đuôi giảm khoảng 40%.
- Kể cả khi đó, tổng năng lượng mỗi đợt **vẫn không bị chặn** bởi ngân sách quay: 10% đợt lớn nhất vượt 40 MeV, trong khi ngân sách chỉ khoảng 26 MeV.
- Ở nhiễu 5 MeV, **nhiễu là nguồn năng lượng chính**, nên việc tắt lái không còn tác dụng.

**Muốn có một ngân sách năng lượng thật sự bị chặn**, cần đồng thời:
1. Tắt lái trong **toàn bộ** đợt, kể cả các pha khóa ngắn giữa đợt.
2. Tắt hoặc giảm mạnh nhiễu trong đợt, hoặc cân bằng năng lượng bằng một bể nhiệt có nhiệt độ hữu hạn.
3. Cố định năng lượng ban đầu $E^*$ theo phản ứng.

## 3. Mức hợp lệ bán cổ điển (để ghi vào Introduction và Conclusion)

| Đại lượng | Giá trị | Ý nghĩa |
|---|---|---|
| Mô-men động lượng pha tại ngưỡng $I\dot\theta_{\text{gãy}}/\hbar$ | **3.2** | Số lượng tử quay ~3. Hiệu chỉnh bán cổ điển cỡ $1/n\sim30\%$ |
| Năng lượng hạt đồng bộ cực đại $\tfrac12m(\kappa\dot\theta_{\text{gãy}})^2$ | **0.57 $\hbar\omega_0$** | Xấp xỉ năng lượng **điểm không** $\hbar\omega_0/2$. Theo cơ học lượng tử, hạt ở trạng thái cơ bản hoặc kích thích thứ nhất; quỹ đạo cổ điển không có nghĩa rõ ràng |
| Xung Carnot TB | **0.40 $\hbar\omega_0$** | **Nhỏ hơn một lượng tử** $\hbar\omega_0$. Một dao động tử lượng tử không thể hấp thụ hay mất 0.4 lượng tử. Trong mô tả lượng tử, các xung này phải được thay bằng chuyển dời rời rạc bội số $\hbar\omega_0$ |

**Đề xuất câu chữ cho bài báo:**

> The model is an effective classical analog: the phase coordinate carries only $\sim3\hbar$ of angular momentum and the bound particle's energy is comparable to its zero-point energy, so quantitative results (pulse energies $\approx0.4\,\hbar\omega_0$, dwell times, spectra) should be read as qualitative trends. A quantum treatment would replace the continuous Carnot losses by transitions in multiples of $\hbar\omega_0$.

## Câu hỏi còn mở

- **Ngân sách cố định.** Muốn so với phổ tiền cân bằng thực nghiệm, cần một lần chạy với năng lượng ban đầu $E^*$ cố định (ví dụ $E^*=S_n+E_{\text{đạn}}$), không có lái và nhiễu trong đợt. Khi đó phân bố $\varepsilon$ mới là một dự đoán có thể so sánh.
- **Lượng tử hóa xung.** Vì xung Carnot nhỏ hơn $\hbar\omega_0$, nên cân nhắc thay va chạm Carnot liên tục bằng quy tắc rời rạc: mất đúng $j\hbar\omega_0$ với xác suất tỉ lệ $\tfrac12m_{\text{eff}}u^2$. Cần kiểm lại bội số hình học và thời gian lưu dưới quy tắc này.
