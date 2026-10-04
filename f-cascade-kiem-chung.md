# Thời gian lưu theo khoảng cách $f$ và "đợt vi chấn = chuỗi gamma cascade": giải và kiểm chứng số liệu

> Tiếp nối [burst-g-kiem-chung.md](burst-g-kiem-chung.md).
> Số liệu tính bằng [giai_so_f_cascade.py](giai_so_f_cascade.py) (chạy `python3 giai_so_f_cascade.py`, ≈ 5.5 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \tau=10,\ \gamma_x=0.70$, nhiễu 1 MeV, $\hbar\omega_0=10$ MeV ($T_0=4.136\times10^{-22}$ s).

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\dot\theta_d=1.05\times0.782\,\dot\theta_{\text{gãy}}=0.821\,\dot\theta_{\text{gãy}}$, khoảng cách 18%, $\tau_{\text{rel}}=1.38\,T_0$ | **Đúng số** (0.8208; 0.1792; 1.375 $T_0$) |
| $\tau_{\text{lưu}}\approx(3$–$4)\,\tau_{\text{rel}}\approx4.8\,T_0$ | **Hệ số 3–4 là khớp, không phải suy ra.** Đo ở $f=1.05$: $\tau_{\text{lưu}}=4.57\,T_0=3.33\,\tau_{\text{rel}}$, nằm trong khoảng. Nhưng theo $f$, tỉ số $\tau_{\text{lưu}}/\tau_{\text{rel}}$ thay đổi **2.05–8.70**, nên "3–4" chỉ đúng quanh $f=1.0$–$1.05$ |
| Tăng $f$ thì $\tau_{\text{lưu}}$ giảm mạnh | **Đúng xu hướng, mức độ vừa phải:** $f=1.05\to1.2$ làm $\tau_{\text{lưu}}$ đổi $4.57\to2.82\,T_0$ (−38%), trong khi khoảng cách tới ngưỡng giảm 3 lần ($0.179\to0.062$) |
| $f<1$ thì $\tau_{\text{lưu}}\to\infty$ | **Gần đúng nhưng không tuyệt đối.** $f=0.95$: vẫn gãy 121 lần, $\tau_{\text{lưu}}=11.96\,T_0$ (hữu hạn). $f=0.90$: chỉ 7 lần gãy trong $480\,T_0$, gần như tắt. Chỗ "tăng vọt" nằm quanh $f\approx0.9$–$0.95$, không phải đúng $f=1$ |
| Cơ chế: thăng giáng vượt khoảng cách $\Omega_0-\dot\theta_d$ | **Mâu thuẫn với kết quả trước:** nếu thời gian lưu do thăng giáng nhiệt phải vượt một khoảng cách, thì nó phải phụ thuộc $k_BT$. Nhưng tốc độ đo được **gần như không đổi** từ 1 tới 8 MeV (xem [kramers-A-kiem-chung.md](kramers-A-kiem-chung.md)). Cơ chế đề xuất giải thích được sự phụ thuộc vào $f$ nhưng không giải thích được việc **không phụ thuộc $k_BT$** |
| Bội số xung mỗi đợt theo phân bố "Poisson/hình học" | **Hình học: ✓** ($\chi^2=10.9$ trên 12 ô, dưới ngưỡng 5% ≈ 18–20; TB 5.12, $p=0.195$). **Poisson: ✗** ($\chi^2=425$). Hai phân bố này khác nhau, nên cần chọn đúng tên: **hình học** |
| Đợt vi chấn = chuỗi gamma cascade của hạt nhân | **Không tương thích về thang thời gian.** Khoảng cách giữa hai xung trong đợt là **$1.47\times10^{-22}$ s**. Thời gian sống Weisskopf của chuyển dời gamma 3–5 MeV: E1 **$2\times10^{-18}$–$6\times10^{-17}$ s**, M1 $1$–$7\times10^{-16}$ s, E2 $4\times10^{-15}$–$1\times10^{-12}$ s. Mô hình nhanh hơn chuyển dời điện từ **4–10 bậc độ lớn** |
| Chuỗi gamma "3–6 photon, mỗi photon 3–5 MeV" | **Khái quát quá mức.** Trong chuỗi gamma thật, tổng năng lượng **bằng năng lượng kích thích** (sau bắt neutron ≈ $S_n\approx6$–$9$ MeV), nên năng lượng mỗi photon thường chỉ vài MeV. Mô hình cho tổng mỗi đợt TB 18 MeV, trải 5–36 MeV ở nhiễu 1 MeV, lớn hơn năng lượng kích thích điển hình |
| "Bức xạ hãm hạt nhân (bremsstrahlung)" | Bremsstrahlung là **phổ liên tục** do gia tốc hạt mang điện. Xung Carnot trong mô hình là tiêu tán rời rạc khi va chạm dẻo. Hai khái niệm không trùng nhau |
| Bảng chốt: $E_{\text{quay}}=51.35$ MeV, $\Delta t_{\text{burst}}=1.2$–$2.6\,T_0$, $\langle\Delta E_{\text{xung}}\rangle=4.0$–$4.6$ MeV | ✓ Khớp các số đo trước |

## 1. Thời gian lưu theo $f=\Omega_0/\dot\theta_{\text{gãy}}$

**Thiết lập:** 8 quỹ đạo × $60T_0$. $\tau_{\text{rel}}=1.375\,T_0$ **cho mọi $f$**, vì $\lvert F'\rvert=2\gamma_\theta(1-\gamma_x\kappa^2/2\gamma_\theta)$ không phụ thuộc $f$.

| $f$ | $\dot\theta_d/\dot\theta_{\text{gãy}}$ | Khoảng cách | Số gãy | Số pha dài | $\tau_{\text{lưu}}$ ($T_0$) | $\tau_{\text{lưu}}/\tau_{\text{rel}}$ | CV phần dư |
|---|---|---|---|---|---|---|---|
| 0.90 | 0.7036 | 0.2964 | **7** | 2 | 7.65 (chỉ 2 mẫu) | 5.57 | — |
| 0.95 | 0.7426 | 0.2574 | 121 | 27 | **11.96** | 8.70 | 0.85 |
| 1.00 | 0.7817 | 0.2183 | 190 | 65 | 5.51 | 4.01 | 1.07 |
| 1.05 | 0.8208 | 0.1792 | 472 | 77 | **4.57** | **3.33** | 1.12 |
| 1.10 | 0.8599 | 0.1401 | 746 | 73 | 3.97 | 2.88 | 1.09 |
| 1.20 | 0.9381 | 0.0619 | 1 279 | 49 | **2.82** | 2.05 | 0.86 |

**Đọc bảng:**
- Chiều phụ thuộc đúng như đề xuất dự đoán: $f$ tăng thì $\tau_{\text{lưu}}$ giảm và số lần gãy tăng mạnh.
- Mức giảm của $\tau_{\text{lưu}}$ yếu hơn nhiều so với mức giảm khoảng cách tới ngưỡng: khoảng cách giảm 3 lần mà $\tau_{\text{lưu}}$ chỉ giảm 1.6 lần.
- Dạng mũ vẫn giữ ở mọi $f$ (CV ≈ 1).

## 2. Bội số xung mỗi đợt ($f=1.05$, 16 quỹ đạo, 121 đợt)

| $N$ | Đếm | Hình học ($p=0.195$) | Poisson dịch (TB−1)+1 |
|---|---|---|---|
| 1 | 23 | 23.6 | 2.0 |
| 2 | 16 | 19.0 | 8.1 |
| 3 | 13 | 15.3 | 16.6 |
| 4 | 18 | 12.3 | 22.9 |
| 5 | 15 | 9.9 | 23.6 |
| 6 | 10 | 8.0 | 19.5 |
| 7 | 3 | 6.4 | 13.4 |
| 8 | 6 | 5.2 | 7.9 |
| 9 | 5 | 4.2 | 4.1 |
| 10 | 1 | 3.3 | 1.9 |
| 11 | 2 | 2.7 | 0.8 |
| ≥12 | 9 | 11.1 | 0.4 |
| **$\chi^2$** | | **10.9** ✓ | **425.1** ✗ |

**Cách hiểu bội số hình học:** sau mỗi lần khóa, có một xác suất gần như cố định $p\approx0.2$ để hệ ở lại pha đồng bộ dài. Với xác suất còn lại, hệ gãy lại ngay.

## 3. So với chuyển dời gamma thật

| | Mô hình | Chuyển dời gamma (Weisskopf, 3–5 MeV, A = 16–208) |
|---|---|---|
| Khoảng cách/thời gian sống giữa hai quanta | $1.47\times10^{-22}$ s | E1: $2.3\times10^{-18}$–$5.8\times10^{-17}$ s; M1: $1.4$–$6.6\times10^{-16}$ s; E2: $3.6\times10^{-15}$–$1.4\times10^{-12}$ s |
| Tổng năng lượng một chuỗi | TB 18 MeV, 5–36 MeV (nhiễu 1 MeV); TB 55 MeV (nhiễu 5 MeV) | Bằng năng lượng kích thích (≈ 6–9 MeV sau bắt neutron) |
| Phụ thuộc nhiệt độ môi trường | Mạnh (tăng ×3 từ 1 lên 5 MeV) | Không, vì tổng năng lượng do trạng thái đầu quyết định |

## 4. Bảng chốt đã sửa

| Thành phần | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Bậc tự do, liên kết, tiên đề tỉ lệ, E gãy, tái khóa, $\Delta E_{\text{xung}}$ | ✓ | ✓ (đã kiểm ở các tài liệu trước) |
| Pha đồng bộ dài | $\tau_{\text{lưu}}\approx(3$–$4)\tau_{\text{rel}}$ do thăng giáng vượt khoảng cách | Dạng mũ ✓; phụ thuộc $f$ ✓ (2.8–12 $T_0$ khi $f=1.2$–0.95). Hệ số "3–4" chỉ đúng ở $f\approx1$–1.05. **Không phụ thuộc $k_BT$, nên cơ chế thăng giáng nhiệt chưa giải thích được** |
| Đợt xả vi chấn | 3–8 xung, 1.2–2.6 $T_0$ | ✓; bội số **hình học** ($p\approx0.2$), không phải Poisson |
| Quanta phát xạ | Tương thích chuỗi gamma cascade | ✗. Nhanh hơn chuyển dời điện từ 4–10 bậc độ lớn; tổng năng lượng phụ thuộc nhiệt độ môi trường |

## Câu hỏi còn mở

- **Vì sao $\tau_{\text{lưu}}$ phụ thuộc $f$ mà không phụ thuộc $k_BT$?** Một khả năng: gãy xảy ra khi **gợn tất định** (do thế $-\Delta E\cos\theta$) cộng với một nhiễu rất nhỏ đưa $\lvert\lambda\rvert$ chạm $R$. Khi đó, xác suất gãy mỗi vòng do hình học của gợn quyết định, chứ không do biên độ nhiễu. Hướng kiểm: giảm nhiễu tới 0.1 MeV và 0.01 MeV. Nếu $\tau_{\text{lưu}}$ vẫn giữ nguyên, khả năng này đúng; nếu $\tau_{\text{lưu}}$ tiến ra vô hạn, thì nhiễu vẫn đóng vai trò quyết định.
- **Nếu không phải chuỗi gamma, thì xung Carnot $\sim10^{-22}$ s ứng với quá trình nào?** Thang $10^{-22}$ s là thang của tương tác mạnh (phát nucleon, va chạm nucleon–nucleon), không phải điện từ. Có thể so với phổ phát nucleon hoặc phổ tiền cân bằng, nhưng cần một kênh phát cụ thể trong mô hình.
