# Rà soát bản thảo Part I (`paper_prc.tex`) so với số liệu kiểm chứng

> **Nguồn:** [gemini-code-1790740027736.txt](gemini-code-1790740027736.txt) (giữ nguyên, không sửa).
> **Bản đã sửa:** [paper_prc.tex](paper_prc.tex), gồm 15 chỗ sửa (30 dòng khác bản gốc).
> Máy không cài LaTeX nên **chưa biên dịch thử**. Đã kiểm tĩnh: ngoặc `{}` cân, `\( \)` cân, các cặp `\begin`/`\end` khớp, không còn ký tự ngoài ASCII.

## 1. Các lỗi cũ đã được bạn sửa đúng ✓

- Tùy chọn `prc`; đã bỏ `showpacs`; email thật.
- Ký hiệu tiếng Anh `dw`, `br`, `dual`.
- $\tau=10/\omega_0=1.59\,T_0$; chú thích $A\sim70$.
- Đã thêm $f$, $k_BT$, $S_n$ vào Bảng I.
- $\dot\theta_{\text{asymp}}=3.000$ và giá trị chính xác 3.205.
- $f\in(0.90,1.139)$.
- Thời gian lưu dạng $T_0+\mathrm{Exp}$ với CV phần dư.
- Phân bố hình học "consistent with" và $p=1/\langle N\rangle$.
- $\dot\theta_{\text{end}}=2.15$; 91.3% thời gian quay thuận; 7/8 quỹ đạo.
- Δn = +2 được nói rõ là số exciton; $n=2k\pm1$.
- Outlook: 16 quỹ đạo, 1.57 bước, 2.05 $\hbar\omega_0$, tính phụ thuộc toán tử ghép.
- Đã thêm tài liệu KD03, Kubo, van der Pol, Gadioli, Breuer.

## 2. Các chỗ mình đã sửa trong `paper_prc.tex`

| # | Vị trí | Bản gốc | Đã sửa | Căn cứ |
|---|---|---|---|---|
| 1 | Tóm tắt | "autonomous … framework" | "driven, noise-activated … framework" | Không nhiễu thì không gãy ([saturation-kiem-chung.md](saturation-kiem-chung.md)) |
| 2 | Tóm tắt | "exact analytical baseline" | "analytical baseline" | Các thống kê là kết quả số với mẫu hữu hạn |
| 3 | §III.A, trước Eq. (8) | "maintains … provided … does not exceed" (trong khi phương trình là $\ge R$) | "…stays below the threshold; breakdown occurs when" | Câu văn và phương trình mâu thuẫn |
| 4 | §III.B | "Total energy balance … drift below $10^{-14}$" | Năng lượng cơ học mất tại mỗi lần khóa bằng Eq. (carnot) tới $10^{-14}$ | Đúng nội dung đã kiểm ([va-cham-tai-lien-ket-kiem-chung.md](va-cham-tai-lien-ket-kiem-chung.md)) |
| 5 | §IV.A | (thiếu) | Thêm một đoạn: hệ **kích hoạt bởi nhiễu** (lề 13.7%, $M_{\det}=0.863$, bão hòa từ ~1 MeV, rào ≈ 0.4 MeV); công lái cân bằng ma sát + Carnot, trong đó **Carnot ≈ 20%** | [saturation-kiem-chung.md](saturation-kiem-chung.md), [phan-bien-part-I-kiem-chung.md](phan-bien-part-I-kiem-chung.md) |
| 6 | §IV.B.1 | "exponentially distributed dwell time $T_0+\mathrm{Exp}$" | Nói rõ $T_0$ là **ngưỡng định nghĩa** của một pha dài | Cách định nghĩa đã dùng trong phân tích |
| 7 | §V.A | $(0.37$–$0.46)\hbar\omega_0$ | $(0.40$–$0.46)\hbar\omega_0$ | 4.0 MeV bằng 0.40 $\hbar\omega_0$ |
| 8 | §V.B.3 (ngân sách) | "decoupling … during the free burst phase enforces $E_{\text{budget}}\le20.4$" | Chỉ ngắt ở pha bay tự do là **không đủ** (80 MeV). Ngắt lái **suốt đợt** cho max 19.9 MeV (1 MeV); nhiễu vài MeV thì ngắt cả bể (max 20.4 MeV). Đổi "E_budget" thành "năng lượng tiêu tán mỗi đợt" | [ba-pha-thoi-gian-kiem-chung.md](ba-pha-thoi-gian-kiem-chung.md), [ngan-sach-ericson-kiem-chung.md](ngan-sach-ericson-kiem-chung.md) |
| 9 | Outlook 2 ($\gamma_x$) | "$\gamma_x=0.70$ corresponds physically to surface absorption $W_D\approx3.5$" | Thêm 3 lưu ý: trung bình theo quỹ đạo còn **0.22–0.45**; vùng liên kết nằm ngoài miền khớp của KD; thay đổi 2.5 lần trong 5–10 MeV. Ghi $\gamma_x=0.70$ là **tham số hiệu dụng**, ứng với $E-\varepsilon_F\approx11$–13 MeV | [wd-trung-binh-kiem-chung.md](wd-trung-binh-kiem-chung.md) (chọn hướng (a)) |
| 10 | Outlook 3 | $\langle x^2\rangle=0.473$ (không nêu tần số cắt) | $0.43$ với $\hbar\omega_c=50$–100 MeV (và $0.473$ với 4.5 MeV) | Để nhất quán với Outlook 2 ([bom-pha-omega-c-kiem-chung.md](bom-pha-omega-c-kiem-chung.md)) |
| 11 | Kết luận | "the framework **proves** …" | "numerical simulations show that, **with $k_BT\gtrsim1$ MeV**, …" | Kết quả số và phụ thuộc nhiễu |
| 12–14 | §I, §IV.A, §V.A | 5 tài liệu có trong danh mục nhưng **không được trích** | Thêm `\cite` cho van der Pol, Kubo, Griffin, Gadioli, Ericson | Tạp chí yêu cầu mọi tài liệu phải được trích |
| 15 | Tài liệu Dao1985 | `C\textohorn{} h\d{o}c gi\h{a}i t\'{\i}ch` | `Co hoc giai tich` | `\textohorn`, `\h` cần font T5, nên revtex/pdflatex sẽ **báo lỗi biên dịch** |

## 3. Việc bạn cần làm hoặc quyết định

- **$\gamma_x$:** mình đã sửa theo hướng (a): giữ 0.70 và ghi là tham số hiệu dụng kèm lưu ý. Nếu muốn hướng (b), tức hạ về khoảng 0.4, thì toàn bộ thống kê ở §IV phải chạy lại.
- **Tài liệu Dao1985:** kiểm lại **nhà xuất bản và năm** trên bìa sách bạn dùng. Nếu muốn giữ dấu tiếng Việt, thêm `\usepackage[T5,T1]{fontenc}` (cần gói vntex) rồi khôi phục dòng gốc.
- **Biên dịch thử** bằng pdflatex hoặc Overleaf để bắt các lỗi mà kiểm tĩnh không thấy.
- **Bảng độ nhạy theo $\gamma_x$** (0.5 / 0.7) nên đưa vào bài, nếu phản biện hỏi.
