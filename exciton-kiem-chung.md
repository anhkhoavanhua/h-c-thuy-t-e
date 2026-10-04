# Tích lũy exciton, cửa sổ vận hành $1.0<f<1.14$ và bảng khung cuối: giải và kiểm chứng số liệu

> Tiếp nối [saturation-kiem-chung.md](saturation-kiem-chung.md).
> Số liệu tính bằng [giai_so_exciton.py](giai_so_exciton.py) (chạy `python3 giai_so_exciton.py`, ≈ 4 phút trên 8 lõi).
> **Bộ thông số:** $I=1,\ m=k=1,\ \Delta E=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \tau=10,\ \gamma_x=0.70$, $\hbar\omega_0=10$ MeV, $S_n=8$ MeV.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $M_{\det}(f=1.14)=1.0$ | **✓** Lưới mịn cho $M_{\det}=1$ tại **$f=1.1392$** |
| $f\ge1.14$: "liên kết không tồn tại ổn định dù chỉ một chu kỳ" | **Sai ngay trên biên.** Không nhiễu, cho phép gãy: tại **$f=1.15$**, hệ gãy tuần hoàn nhưng mỗi pha khóa kéo dài **1.93 $T_0$** (19/19 pha ≥ $T_0$). Chỉ từ $f\approx1.2$ trở lên pha khóa mới ngắn hơn $T_0$ (0.47 $T_0$ tại 1.20; 0.015 $T_0$ tại 1.30). $f\ge1.14$ là vùng **gãy tất định tuần hoàn**, không phải "không tồn tại" |
| Cận dưới $f>1.0$ | **Không có căn cứ từ mô phỏng.** Tại $f=0.95$ vẫn có 121 lần gãy và $\tau_{\text{lưu}}=12\,T_0$ ([f-cascade-kiem-chung.md](f-cascade-kiem-chung.md)). Vùng gãy do nhiễu kéo xuống tới $f\approx0.9$ |
| Mỗi xung $\approx4.0$–$4.6$ MeV; "$E_{\text{exc}}(1)\approx4.2$ MeV $<S_n$" | **Chỉ đúng về trung bình.** Phân bố xung rất rộng: ở 1 MeV, trung vị 3.91, khoảng 10–90% là **0.00–9.15 MeV**. **14.5%** số xung **đơn lẻ đã ≥ $S_n$**. Ở 5 MeV: trung vị **0.03**, khoảng 10–90% **0–14.8 MeV**, và 23.7% số xung ≥ $S_n$ |
| "$k=2$: $E_{\text{exc}}(2)\approx8.4$ MeV ≥ $S_n$, kênh thoát mở" | **Sai.** $E_{\text{exc}}(2)$ trung bình **6.43 MeV** (trung vị 6.44). Chỉ **32.7%** đợt vượt $S_n$ sau 2 xung. Bước mở kênh **điển hình là $k=3$** (trung vị 3, TB 2.69) |
| Đợt $N\ge3$ chiếm "≈ 50%" | **Sai số học.** Với phân bố hình học $p=0.195$: $(1-p)^2=$ **65%**. Đo được **61.7%** (1 MeV), **79.4%** (5 MeV) |
| $\varepsilon_{\text{kin}}=E_{\text{exc}}-S_n\approx4$–$15$ MeV | **Phân bố rộng hơn nhiều.** Ở 1 MeV: TB **18.3**, khoảng 10–90% **1.6–44.6 MeV**. Ở 5 MeV: TB **61.4**, khoảng 10–90% **4.4–142 MeV**. Ngoài ra, $E_{\text{exc}}-S_n$ chỉ là **cận trên**: trong mô hình exciton, năng lượng được chia cho nhiều exciton, nên nucleon phát ra chỉ nhận một phần |
| "Khớp tự nhiên với phổ neutron tiền cân bằng thực nghiệm" | **Chưa kiểm.** Không có phổ thực nghiệm nào được đưa ra để so. Thêm vào đó: (1) phổ mô hình **phụ thuộc nhiệt độ môi trường** (TB tăng ×3.4 từ 1 lên 5 MeV); (2) động năng tới 142 MeV **vượt xa** độ sâu hố thế ($E_{\text{vt}}=30$ MeV) và năng lượng kích thích điển hình. Trong phản ứng thật, năng lượng kích thích do **đạn bắn** quyết định. Ở đây nó do nguồn lái và nhiễu quyết định |
| "Khép kín, bảo toàn cấu trúc symplectic" | **Sai.** Mô hình có mô-men lái Van der Pol, cản Rayleigh, nhiễu Langevin và va chạm dẻo. Tất cả đều **phá** cấu trúc symplectic (hệ không còn Hamilton). Chỉ phần Lagrange gốc (1)–(2) là symplectic |
| Các mục còn lại của bảng khung (liên kết, tiên đề, E gãy, $\tau_{\text{lưu}}$ dạng mũ, bội số hình học, $\Delta E_{\text{xung}}$ TB) | ✓ như các tài liệu trước |

## 1. Lề ứng suất tất định quanh biên

| $f$ | 1.10 | 1.11 | 1.12 | 1.13 | **1.14** | 1.15 | 1.16 | 1.17 | 1.18 | 1.19 | 1.20 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $M_{\det}$ | 0.938 | 0.954 | 0.969 | 0.985 | **1.001** | 1.017 | 1.034 | 1.050 | 1.067 | 1.084 | 1.101 |

Nội suy tuyến tính: $M_{\det}=1$ tại $f=1.1392$.

## 2. Cho phép gãy, không nhiễu ($40T_0$ sau khởi động)

| $f$ | Số gãy | Pha khóa: dài nhất / TB | Số pha ≥ $T_0$ |
|---|---|---|---|
| 1.05 | 0 | khóa suốt | — |
| 1.10 | 0 | khóa suốt | — |
| **1.15** | 19 | 1.926 / 1.926 $T_0$ | **19/19** |
| 1.20 | 58 | 0.471 / 0.470 $T_0$ | 0 |
| 1.30 | 146 | 0.056 / 0.015 $T_0$ | 0 |

Ngay trên biên, hệ vào một **chu trình gãy tất định** có chu kỳ ≈ 1.93 $T_0$, gần như tuần hoàn chính xác. Pha đồng bộ vẫn tồn tại, chỉ ngắn dần khi $f$ tăng.

## 3. Tích lũy năng lượng exciton trong đợt ($f=1.05$, 16 quỹ đạo × $60T_0$)

| Chỉ số | Nhiễu 1 MeV | Nhiễu 5 MeV |
|---|---|---|
| Số đợt | 128 | 107 |
| Tỉ lệ $N\ge3$ (đo / hình học) | 61.7% / 64.0% | 79.4% / 86.7% |
| Xung đơn: TB / trung vị / 10–90% | 4.01 / 3.91 / 0.00–9.15 MeV | 4.35 / 0.03 / 0.00–14.78 MeV |
| $P(\text{xung đơn}\ge S_n)$ | **14.5%** | **23.7%** |
| $E_{\text{exc}}(2)$: TB / trung vị | 6.43 / 6.44 MeV | 7.11 / 4.24 MeV |
| $P(E_{\text{exc}}(2)\ge S_n)$ | 32.7% | 37.5% |
| Đợt có tổng ≥ $S_n$ | 70.3% | 90.7% |
| Bước mở kênh $k$: TB / trung vị | 2.69 / **3** | 3.68 / **3** |
| $\varepsilon_{\max}=E_{\text{exc}}-S_n$: TB / 10–90% | 18.3 / 1.6–44.6 MeV | 61.4 / 4.4–141.5 MeV |

**Đọc bảng:**
- Bức tranh "mỗi xung ~4.2 MeV, hai xung vừa đủ vượt $S_n$" không đúng.
- Phân bố xung **rất lệch**: có nhiều xung gần bằng 0, xen với các xung lớn ≥ 8 MeV.
- Một số đợt mở kênh thoát ngay ở xung đầu tiên (14.5%). Phần lớn đợt cần 3 xung.
- Năng lượng khả dụng cho phát xạ **không bị chặn** bởi một năng lượng kích thích xác định. Nó tăng theo nhiễu và theo độ dài đợt.

## 4. Bảng khung đã sửa

| Hạng mục | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Không gian pha | "Khép kín, bảo toàn symplectic" | Phần Lagrange gốc là symplectic. Mô hình đầy đủ (lái + cản + nhiễu + va chạm dẻo) **không** symplectic |
| Vùng vận hành | $1.0<f<1.14$, chọn 1.05 | Biên trên 1.139 ✓ (bắt đầu gãy tất định; pha khóa vẫn ~1.9 $T_0$ ngay trên biên). **Biên dưới 1.0 không có căn cứ**: gãy do nhiễu còn tới $f\approx0.9$ |
| Cơ chế phát xạ | Mỗi xung một bước exciton; ≥ 2 xung mở kênh | Mỗi xung một bước ✓ về cách hiểu. **Điển hình cần 3 xung**; 14.5% xung đơn đã đủ. Phổ $\varepsilon$ rộng 1.6–45 MeV và phụ thuộc nhiệt độ. Chưa có so sánh với phổ thực nghiệm |

## Câu hỏi còn mở

- **Ràng buộc năng lượng.** Trong phản ứng hạt nhân thật, tổng năng lượng kích thích bị cố định bởi đạn bắn (hoặc bởi $S_n$ khi bắt neutron). Mô hình hiện nhận năng lượng liên tục từ nguồn lái Van der Pol, nên $E_{\text{exc}}$ không bị chặn. Nếu muốn so với phổ tiền cân bằng, cần thêm một **ngân sách năng lượng cố định** cho mỗi phản ứng, rồi mới tính phổ.
- **Chia năng lượng giữa các exciton.** Để có phổ nucleon thật sự, cần mô hình hóa cách năng lượng $E_{\text{exc}}$ chia giữa $n$ exciton, ví dụ theo mật độ trạng thái Ericson. Chỉ dùng cận trên $E_{\text{exc}}-S_n$ là chưa đủ.
