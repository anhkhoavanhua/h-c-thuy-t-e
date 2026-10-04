# Kim $\hat K$ là đồng hồ ngoài hạt nhân: kiểm số đề xuất "chuỗi xung attosecond"

> **Nối tiếp:** [hai-vai-tro-kim-va-be-kiem-chung.md](hai-vai-tro-kim-va-be-kiem-chung.md) (câu hỏi còn mở) và [nam-dong-ho-luong-tu-kiem-chung.md](nam-dong-ho-luong-tu-kiem-chung.md).
> **Số liệu:** [giai_so_kim_ngoai.py](giai_so_kim_ngoai.py), kết quả thô [kim_ngoai_ket_qua.txt](kim_ngoai_ket_qua.txt), chạy khoảng 1.5 phút.
> **Mô hình:** cùng hệ S 36 chiều và Lindblad $\mathcal E_{\text{dual}}$ (Carnot + bơm) như hai file trên. Liouvillian 1296×1296 được chéo hóa chính xác. Khe phổ là 0.0503 ω0, ứng với thời gian hồi phục 3.17 $T_0$.

## Câu hỏi được kiểm

Bản đề xuất gồm hai phần:
1. Trong hình thức luận, $\hat K$ là đồng hồ chuẩn hình thức (fiducial), tương tự thời gian tọa độ trong thuyết tương đối rộng.
2. Trong thực nghiệm, $\hat K$ là một nguồn kết hợp **bên ngoài**, "ví dụ chuỗi xung laser quang học/tia X attosecond".

Phần 1 là một lựa chọn diễn giải nên không kiểm bằng số được. Phần 2 kiểm được. Câu hỏi cụ thể: một đồng hồ ngoài có độ rung $\sigma_K$ thì còn cho ra phương trình Lindblad không, và cần $\sigma_K$ nhỏ đến đâu?

## Tóm tắt

| Mệnh đề | Kiểm số | Kết quả |
|---|---|---|
| Chuỗi xung attosecond làm được $\hat K$ | Yêu cầu 0.071 $T_0$ = **29.4 ys** (yoctogiây). Xung ngắn nhất đã đo là 43 as, dài gấp **1.46×10⁶** lần yêu cầu. Chu kỳ chuỗi xung khoảng 1.33 fs, gấp **4.5×10⁷** lần | ✗ **sai 6–7 bậc** |
| Đồng hồ ngoài có rung **không tích lũy** vẫn cho Lindblad | Nửa nhóm $\tilde\rho(\hat t+T_0)=e^{T_0\mathcal L}\tilde\rho(\hat t)$ đúng tới **≤ 8×10⁻¹⁶** với σ_K ≤ 0.071 $T_0$. Rung chỉ làm nhòe điều kiện đầu, **không đổi ℒ** | ✓ luật chuyển động giữ nguyên |
| Đồng hồ ngoài có rung **tích lũy** (dao động tự chạy) | $\tilde\rho(t)=e^{t(\mathcal L+\frac D2\mathcal L^2)}\rho_0$ đúng tới **10⁻¹⁵**. Đồng hồ thật thêm một số hạng giải kết hợp $\frac D2\mathcal L^2$ (Gambini–Porto–Pullin) | ✓ ℒ **bị sửa** |
| σ_t ≤ 0.071 $T_0$ là đủ | Đủ để **xếp thứ tự** ba pha (99%). Để dựng lại **trạng thái** ρ_S(t) trong nhiễu thống kê 0.016 thì cần **σ_K ≤ 0.003 $T_0$ = 1.2 ys**, chặt hơn 24 lần | yêu cầu thật còn khắt khe hơn |
| Năng lượng của đồng hồ ngoài | Băng thông $\hbar/2\sigma$: **11.2 MeV** (mức xếp thứ tự), **265 MeV** (mức trạng thái). Photon có ít nhất một chu kỳ trong nhịp τ = T0/25: **250 MeV**, vượt ngưỡng sinh pion 140 MeV | ✗ không thể là ánh sáng quang học/tia X |

## A. Yêu cầu của $\hat K$ theo đơn vị vật lý

| Đại lượng | Giá trị |
|---|---|
| $T_0$ | 4.136×10⁻²² s = 413.6 ys |
| σ_t cho 99% xếp đúng ba pha (0.071 $T_0$) | 2.94×10⁻²³ s = **29.4 ys** |
| Nhịp τ = $T_0$/25 | 1.65×10⁻²³ s = 16.5 ys |
| Băng thông tối thiểu ħ/(2σ_t) | **11.2 MeV** |
| Photon có ≥ 1 chu kỳ trong τ: h/τ | **250 MeV** (ngân sách kích thích 20–30 MeV; ngưỡng pion ≈ 140 MeV) |

| Nguồn "attosecond" | Độ dài | Theo $T_0$ | Gấp mấy lần yêu cầu |
|---|---|---|---|
| Xung cô lập ngắn nhất đã đo (43 as) | 4.3×10⁻¹⁷ s | 1.04×10⁵ $T_0$ | 1.46×10⁶ |
| Chu kỳ chuỗi xung (nửa chu kỳ laser 800 nm) | 1.33×10⁻¹⁵ s | 3.2×10⁶ $T_0$ | 4.5×10⁷ |
| Photon XUV 100 eV | 10⁻⁴ MeV | — | băng thông chỉ bằng 9×10⁻⁶ lần mức cần |

Thang thời gian của học thuyết là **yoctogiây**, không phải attogiây. Chuỗi xung attosecond quét qua khoảng 10⁵ chu kỳ $T_0$ trong mỗi xung.

## B. Rung không tích lũy: luật Lindblad giữ nguyên, điều kiện đầu bị nhòe

**Mô hình.** Mỗi lần đọc kim lệch $\varepsilon\sim N(0,\sigma_K^2)$, ví dụ do độ dài hoặc độ trễ của xung. Nếu xung đến trước khi kích thì hệ vẫn ở ρ0. Trạng thái gán cho số đọc $\hat t$ là:

$$\tilde\rho(\hat t)=\mathbb E_\varepsilon\big[\rho(\max(0,\hat t+\varepsilon))\big].$$

Kỳ vọng được tính bằng cầu phương Gauss–Hermite 120 nút. Số đọc lấy mỗi 0.04 $T_0$ trên 5 $T_0$.

| σ_K / $T_0$ | σ_K (s) | max‖ρ̃ − ρ‖ | Biên độ động lực max‖ρ̃(t̂₁) − ρ̃(t̂₂)‖ (lý tưởng 0.981) | Nửa nhóm ‖ρ̃(t̂+T0) − e^{T0ℒ}ρ̃(t̂)‖ |
|---|---|---|---|---|
| 0.001 | 4.1×10⁻²⁵ | 0.0055 | 0.981 | 8×10⁻¹⁶ |
| 0.002 | 8.3×10⁻²⁵ | 0.0109 | 0.981 | 7×10⁻¹⁶ |
| **0.003** | **1.2×10⁻²⁴** | **0.0164** ≈ nhiễu thống kê | 0.980 | 6×10⁻¹⁶ |
| 0.010 | 4.1×10⁻²⁴ | 0.054 | 0.978 | 8×10⁻¹⁶ |
| 0.071 | 2.9×10⁻²³ | 0.418 | 0.890 | 6×10⁻¹⁶ |
| 0.27 | 1.1×10⁻²² | 0.745 | 0.618 | 1×10⁻⁷ |
| 1 | 4.1×10⁻²² | 0.895 | 0.507 | 1×10⁻⁴ |
| 5 | 2.1×10⁻²¹ | 0.875 | 0.288 | 6×10⁻⁴ |
| 104 (= 43 **zs**, vẫn ngắn hơn 43 as khoảng 1000 lần) | 4.3×10⁻²⁰ | 0.875 | **0.000** | 2×10⁻³ |

1. **ℒ không phụ thuộc chất lượng kim.** Với $\hat t>3\sigma_K$, ta có $\tilde\rho(\hat t)=e^{\hat t\mathcal L}\,\mathbb E[e^{\varepsilon\mathcal L}]\rho_0$, đúng tới máy tính khi σ_K ≤ 0.071 $T_0$. Phần dư nhỏ ở σ lớn đến từ việc kẹp t ≥ 0. Kim rung không làm sai **phương trình chuyển động**, chỉ làm sai **điều kiện đầu**. Điều này củng cố mệnh đề V1a của tài liệu gốc: ℒ đến từ bể, tham số t đến từ kim.
2. **Nhưng nhận ra trạng thái thì cần kim rất mịn.** Dao động kết hợp nhanh của rô-to (cỡ vài ω0) bị trung bình hóa, nên sai số trạng thái vượt nhiễu thống kê ngay khi σ_K > 0.003 $T_0$.
3. **Ở σ_K = 104 $T_0$, biên độ động lực bằng 0.000.** Mọi số đọc cho cùng một trạng thái, nên kim không còn thông tin thời gian nào. Một xung 43 as dài gấp 1000 lần mức này.

## C. Rung tích lũy: đồng hồ thật sửa phương trình Lindblad

**Mô hình.** Đồng hồ tự chạy có nhiễu pha, $\mathrm{Var}[\hat t-t]=Dt$. Khi đó

$$\tilde\rho(t)=\mathbb E\big[e^{(t+W_t)\mathcal L}\big]\rho_0=e^{t(\mathcal L+\frac D2\mathcal L^2)}\rho_0,$$

vì ℒ giao hoán với chính nó và hàm sinh mô-men của biến Gauss cho đúng dạng này.

| D / $T_0$ | σ_t tại 5 $T_0$ | ‖cầu phương − e^{t(ℒ+Dℒ²/2)}ρ0‖ | max‖ρ̃ − ρ‖ |
|---|---|---|---|
| 10⁻⁴ | 0.022 $T_0$ | 1.0×10⁻¹⁵ | 0.0074 |
| 10⁻³ | 0.071 $T_0$ | 1.0×10⁻¹⁵ | 0.062 |
| 3×10⁻³ | 0.122 $T_0$ | 1.8×10⁻¹⁰ | 0.141 |
| 10⁻² | 0.224 $T_0$ | 9.5×10⁻⁵ | 0.276 |
| 10⁻¹ | 0.707 $T_0$ | 7.8×10⁻² | 0.594 |

- Đẳng thức đúng tới máy tính. Lệch ở D lớn chỉ do kẹp t ≥ 0.
- Một đồng hồ vật lý tự chạy **luôn** thêm số hạng $\frac D2\mathcal L^2$. Đây là một kênh giải kết hợp mới, không đến từ bể Carnot.
- Muốn giữ sai số ≤ 0.016 trên 5 $T_0$ cần D ≲ 2.5×10⁻⁴ $T_0$.

## Kết luận: $\hat K$ thực chất là gì

1. **Ví dụ "chuỗi xung laser/tia X attosecond" phải bỏ.** Nó quá thô 6–7 bậc. Với độ rung như vậy, biên độ động lực đọc được bằng 0 (mục B), và photon XUV thiếu năng lượng 5 bậc. Thang của học thuyết là yoctogiây và MeV.
2. **Một $\hat K$ vật lý bên ngoài, nếu có, phải là xung γ cỡ zeptô/yoctogiây với băng thông ≥ 11 MeV.**
   - Mức 11 MeV đủ để xếp thứ tự ba pha.
   - Muốn dựng lại trạng thái thì cần ≥ 265 MeV, vượt ngưỡng sinh pion và phá hạt nhân.
   - Nếu dùng sơ đồ bơm–dò, yêu cầu này đặt lên cả độ dài xung lẫn độ chính xác của độ trễ.
   - Công nghệ hiện nay chưa có nguồn nào như vậy.
3. **Vì vậy, ở cấp phương trình Lindblad, $\hat K$ chỉ có thể là đồng hồ chuẩn hình thức (tham số toán học).** Diễn giải "thời gian tọa độ" của bản đề xuất đứng được, và mục B cho nó một cơ sở số: rung không tích lũy không đổi ℒ. Nhưng một đồng hồ vật lý thật thì đổi ℒ thành ℒ + (D/2)ℒ² (mục C). Do đó "σ_t = 0" chỉ đúng với kim lý tưởng, không đúng với bất kỳ thiết bị nào.
4. **Phân vai giữ nguyên.** Bể cho chiều; $\hat K$ cho số đo. Nhưng "số đo" ở độ phân giải của học thuyết hiện chỉ tồn tại như tham số toán học.

**Câu tiếng Anh đề xuất (thay cho ví dụ attosecond):**

> K̂ is a fiducial reference clock, the analogue of coordinate time: timing jitter that does not accumulate leaves the conditioned generator ℒ unchanged (semigroup property verified to 10⁻¹⁵) and only smears the initial state, while a free-running physical clock with phase diffusion D replaces ℒ by ℒ + (D/2)ℒ². Resolving the engine's states requires σ_t ≲ 0.003 T₀ ≈ 1 ys (bandwidth ≳ 0.27 GeV), and ordering its three phases requires σ_t ≲ 0.07 T₀ ≈ 30 ys (≳ 11 MeV); existing attosecond sources are six orders of magnitude too slow. Inside the nucleus only the coarse Carnot clock exists: the bath fixes the direction, K̂ the measure.

## Câu hỏi còn mở

- Có muốn thêm số hạng $\frac D2\mathcal L^2$ vào Part III như "giới hạn đồng hồ thực" không? Nếu có, cần chọn D từ một mô hình nguồn γ cụ thể.
- Ngưỡng 0.003 $T_0$ phụ thuộc trạng thái đầu ρ0, ở đây là bó sóng rô-to ℓ ≈ 3. Có cần kiểm thêm với ρ0 là trạng thái riêng hoặc trạng thái dừng không? Với trạng thái dừng thì ngưỡng sẽ lỏng hơn nhiều.
