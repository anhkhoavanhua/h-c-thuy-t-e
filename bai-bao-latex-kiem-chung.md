# Kiểm chứng bản thảo LaTeX "Classical Analytical Mechanics Formulation of a Self-Sustained Dual-Engine Coupling…"

> Đối chiếu từng công thức và con số với các tài liệu kiểm chứng trong thư mục này. Không có số liệu mô phỏng mới, trừ vài phép đổi đơn vị.
> Máy không cài LaTeX, nên chưa biên dịch thử. Các lỗi biên dịch dự kiến được nêu ở mục 5.

## Tóm tắt

- **Phần giải tích (Mục II–III, các phương trình (1)–(11)): đúng toàn bộ.** Đã kiểm bằng sympy và mô phỏng ở các tài liệu trước.
- **Bảng I có 3 lỗi đơn vị.**
- **Mục IV–V có 7 chỗ cần sửa:** số liệu lấy từ các lần chạy khác nhau, một chỗ lập luận vòng tròn, $n$ không khớp với $k$, và bỏ sót phần giới hạn đã biết.
- **Phần hình thức:** vài lỗi LaTeX và placeholder cần thay.

## 1. Mục II–III: phương trình

| Chỗ trong bài | Kiểm chứng | Nguồn |
|---|---|---|
| (1) Lagrangian | ✓ | [he-e-kep-hat-vat-chat-lagrange.md](he-e-kep-hat-vat-chat-lagrange.md) |
| (2) $f=x-\kappa\sin\theta$, $\kappa=\ell_0\Delta E/E_{\text{vt}}$ | ✓ ($=1/3$ khi $E_{\text{vt}}=3\Delta E$) | |
| (3)–(4) Lagrange loại một | ✓ (sympy) | [lien-ket-dong-bo-e-gay-kiem-chung.md](lien-ket-dong-bo-e-gay-kiem-chung.md) |
| (5) phương trình rút gọn | ✓ (sympy) | nt |
| (6) $\lambda(\theta,\dot\theta)$ dạng đóng | ✓ (sympy; khớp DAE $4\times10^{-12}$) | nt |
| (7) $\lvert\lambda\rvert\ge R$ | ✓ | |
| (8) $\dot\theta_{\text{gãy}}\approx\sqrt{R/m\kappa}=E_{\text{vt}}/(\ell_0\sqrt{m\Delta E})$ | ✓ đại số. **Nên ghi thêm giá trị chính xác:** tiệm cận cho **3.000**, chính xác là **3.205** (lệch 6.4%) | [chi-va-dai-hoi-phuc-kiem-chung.md](chi-va-dai-hoi-phuc-kiem-chung.md), [coherence-va-g-kiem-chung.md](coherence-va-g-kiem-chung.md) |
| (9) Liên tục vị trí và xung lượng; $\Delta\ddot x=-\lambda/m$, $\Delta\ddot\theta=\lambda\kappa\cos\theta/I$ | ✓ (khớp số tới 6 chữ số) | [lien-ket-dong-bo-e-gay-kiem-chung.md](lien-ket-dong-bo-e-gay-kiem-chung.md) |
| (10)–(11) $\Lambda=m_{\text{eff}}u$, $m_{\text{eff}}$, Carnot $\tfrac12m_{\text{eff}}u^2$ | ✓ (sympy; khớp mô phỏng tới $10^{-14}$) | [va-cham-tai-lien-ket-kiem-chung.md](va-cham-tai-lien-ket-kiem-chung.md) |

## 2. Bảng I

| Dòng | Trong bài | Kiểm chứng | Sửa thành |
|---|---|---|---|
| Caption "$A=56$, $\hbar\omega_0=10.0$ MeV" | | $41\cdot56^{-1/3}=$ **10.72 MeV**. Giá trị 10 MeV ứng với $A\approx69$ | Bỏ "$A=56$", hoặc ghi "$\hbar\omega_0=10$ MeV (≈ $A\sim70$)" |
| $\hbar\omega_0$, $T_0$, $\Delta E$, $I$, $E_{\text{vt}}$, $\ell_0$, $\kappa$, $R$, $g$ | | ✓ tất cả (4.136e-22 s; 6.941e-57 kg·m²; 2.037 fm; 0.679 fm; 2359.4 N; 78.65 N) | — |
| **Drive damping $\tau=I/\gamma_\theta$ = "10.0 $T_0$"** | | **Sai đơn vị.** Giá trị chuẩn hóa 10 là theo $1/\omega_0$, nên bằng $10/2\pi=$ **1.59 $T_0$** | `\(10/\omega_0 = 1.59\,T_0\)` |
| **Nucleon friction $\gamma_x$ = "$m/(1.27\,T_0)$", 0.700** | | **Sai hai chỗ.** (1) 1.27 là theo $1/\omega_0$ (= **0.20 $T_0$**), không phải $T_0$. Nếu thật là $m/(1.27T_0)$ thì giá trị chuẩn hóa là 0.125. (2) $m/\tau_{\text{relax}}$ với $\tau_{\text{relax}}=1.267/\omega_0$ cho **0.789**, không phải 0.700. Giá trị 0.70 đến từ giả thiết $\lambda_{\text{mfp}}\propto A^{1/3}$ (không có căn cứ, xem [kramers-A-kiem-chung.md](kramers-A-kiem-chung.md)) | Ghi `\(0.70\)` là **tham số chọn**, hoặc dùng `\(m/\tau_{\rm relax}=0.789\)` rồi chạy lại số liệu |
| Thiếu | | Bảng không có $f$ (=1.05), $k_BT$ (=1 MeV), $\Omega_0$, $S_n$. Mọi thống kê ở Mục IV phụ thuộc mạnh vào $f$ và $k_BT$ | Thêm 4 dòng này |

## 3. Mục IV–V: số liệu và lập luận

| # | Chỗ trong bài | Vấn đề | Đề xuất sửa |
|---|---|---|---|
| 1 | (12)–(13) phương trình Langevin | Đây là phương trình của **pha tự do**. Trong pha khóa còn có $\lambda$ và liên kết (2), nhưng bài không nói | Thêm câu: "in the locked phase Eqs. (12)–(13) are supplemented by $\lambda$ and constraint (2)" |
| 2 | "$f\in[1.00,1.14]$" | Biên trên đo được là **1.139** ✓. **Biên dưới 1.00 không có căn cứ**: gãy do nhiễu còn xảy ra tới $f\approx0.9$ ([f-cascade-kiem-chung.md](f-cascade-kiem-chung.md)). Giá trị dùng trong bài ($f=1.05$) không được nêu | `\(f\in(0.9,\,1.139)\)`, ghi rõ $f=1.05$ |
| 3 | "$\tau_{\text{lưu}}\approx4.57\,T_0$ (CV ≈ 1.006)" | **Ghép số từ hai lần chạy khác nhau.** 4.57 $T_0$ là từ [f-cascade-kiem-chung.md](f-cascade-kiem-chung.md) (CV phần dư 1.12). CV 1.006 là từ [rice-poisson-kiem-chung.md](rice-poisson-kiem-chung.md) (TB 5.54 $T_0$). Ngoài ra CV ≈ 1 là của **phần dư** $w-T_0$, không phải của $\tau_{\text{lưu}}$ (CV bản thân $\tau_{\text{lưu}}$ là 0.82). Các lần chạy cho $\tau_{\text{lưu}}=4.5$–$6.1\,T_0$ (sai số thống kê ~15%) | "dwell time $T_0+\mathrm{Exp}$, mean $4.5$–$6\,T_0$, residual $CV\approx1.0$–$1.1$" |
| 4 | "sub-impacts occur at short intervals ($<0.2\,T_0$)" | 69–82% số lần gãy lại xảy ra trong $<0.2T_0$, nhưng khoảng cách **trung bình** giữa hai lần khóa trong một đợt là **0.355 $T_0$** | "predominantly within $0.2\,T_0$ (mean spacing $0.36\,T_0$)" |
| 5 | "**strictly** follows geometric… $\langle N\rangle=1/p=5.12$ ($\chi^2=10.9$)" | **Vòng tròn:** $p$ được **tính từ** $\langle N\rangle$, nên $1/p=5.12$ là đồng nhất thức. "strictly" cũng quá mạnh: kiểm định chỉ cho thấy dữ liệu **không bác bỏ** phân bố hình học | "is consistent with a geometric law ($\chi^2=10.9$ over 12 bins); the maximum-likelihood $p=1/\langle N\rangle=0.195$" |
| 6 | "$\dot\theta_{\text{end}}\approx2.15$" | ✓ Đo được TB 2.153 ([recharge-kiem-chung.md](recharge-kiem-chung.md)). "dưới ngưỡng" đúng 100% | — |
| 7 | "Δn = +2" | ✓ **Đúng theo nghĩa số exciton** $n=p+h$ | Nên viết rõ "exciton number $n=p+h$ increases by 2" để không nhầm với số lượng tử dao động tử |
| 8 | "$k=3$, $n=4$" | **Không khớp với Δn = +2.** Xuất phát từ $n=1$ ($1p0h$) thì $k=3$ cho $n=7$; xuất phát từ $n=2$ ($1p1h$) thì $n=6$. $n=4$ chỉ ứng với $k=2$ từ $n=2$. Con số 0.9–1.6 MeV lấy từ tính trên từng đợt với $n=2k\pm1$ ([ngan-sach-ericson-kiem-chung.md](ngan-sach-ericson-kiem-chung.md)), không phải từ $n=4$ cố định (với $n=4$, $E_{\max}=3.78$ MeV: 1.51 MeV) | Ghi $n=2k+1$ (hoặc $2k$) và "$0.9$–$1.6$ MeV averaged over simulated bursts" |
| 9 | "typically achieved at step $k\ge3$" | Trung vị $k=3$ ✓. Nhưng 14.5% xung đơn đã ≥ $S_n$, và 32.7% đợt đạt $S_n$ sau 2 xung | "median $k=3$" |
| 10 | Thiếu | Bài không nói **tổng năng lượng mỗi đợt không bị chặn** (đuôi tới 73–286 MeV khi không ngắt), và các thống kê **phụ thuộc $k_BT$** (không nhiễu thì không có gãy khi $f<1.139$; tốc độ bão hòa từ ~1 MeV). Phản biện sẽ hỏi | Thêm một đoạn "Limitations", dựa trên [ba-pha-thoi-gian-kiem-chung.md](ba-pha-thoi-gian-kiem-chung.md) và [saturation-kiem-chung.md](saturation-kiem-chung.md) |
| 11 | Mục V "derived dynamical cycles… represent **exact** features", Kết luận "rigorous", "**exact** mathematical foundation" | Các phương trình là chính xác, nhưng mọi thống kê ($\tau_{\text{lưu}}$, $p$, $\langle\varepsilon\rangle$) là **kết quả số hữu hạn mẫu** (64–1400 sự kiện) phụ thuộc tham số chọn tay | "numerically established (finite-sample) features" |
| 12 | Kết luận "Poisson-governed relaxation cascades" | Đợt (cascade) có bội số **hình học**. Chỉ thời gian lưu mới có dạng Poisson (mũ) | "exponentially distributed dwell times and geometrically distributed cascades" |
| 13 | Tóm tắt "internal exciton dissipation **quantum**" | Mâu thuẫn với chính câu sau, vốn nói xung nhỏ hơn một lượng tử ($0.4\hbar\omega_0$) | "dissipation step" |

## 4. Tài liệu tham khảo

| Mục | Kiểm |
|---|---|
| Ericson 1960, Adv. Phys. **9**, 425 | ✓ có thật |
| Griffin 1966, PRL **17**, 478 | ✓ có thật |
| Goldstein, Poole, Safko 3rd ed. | ✓ (năm in thường ghi 2002 hoặc 2001) |
| Brogliato 1999, Springer | ✓ (bản 2) |
| **N. V. Dao, *Analytical Mechanics*, 1985** | ⚠ Cần đối chiếu với **bìa sách bạn dùng** ("Cơ học giải tích", Nguyễn Văn Đạo): tên tiếng Việt, nhà xuất bản, năm. Nên ghi tên gốc tiếng Việt kèm dịch |
| Thiếu | Van der Pol (dao động tự kích), Langevin/Kubo (thăng giáng–tiêu tán), Blann hoặc Gadioli–Hodgson (mô hình exciton tiền cân bằng). Các cơ chế này dùng trong bài nhưng chưa được trích dẫn |

## 5. Hình thức và biên dịch

| Chỗ | Vấn đề | Sửa |
|---|---|---|
| `\text{lưu}`, `\text{gãy}`, `\mathcal{E}_{\text{kép}}` | Ký tự tiếng Việt có dấu trong revtex/pdflatex **không có** `inputenc`/`fontenc` T5 sẽ **lỗi biên dịch** hoặc mất dấu. Tạp chí quốc tế cũng yêu cầu ký hiệu tiếng Anh | `\tau_{\rm dw}`, `\dot\theta_{\rm br}`, `\mathcal{E}_{\rm dual}` |
| `prd` | Physical Review **D** là tạp chí hạt/trường/hấp dẫn. Nội dung hạt nhân phù hợp **PRC** | `prc` |
| `showpacs` | PACS đã ngừng từ 2016 | Bỏ |
| `contact@example.com` | Placeholder | Email thật |
| Tiêu đề và tóm tắt | "Self-Sustained" và "autonomous" không đúng hoàn toàn: hệ có nhiễu ngoài, và không có nhiễu thì không có gãy khi $f<1.139$ | Cân nhắc "Noise-Activated" |

## Việc tiếp theo đề xuất

Mình có thể xuất một file `.tex` đã áp các sửa ở mục 2–5, và giữ nguyên phần giải tích đúng. Bạn chỉ cần báo nếu muốn.
