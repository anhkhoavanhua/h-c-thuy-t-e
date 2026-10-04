# Thế quang học, "rời rạc hóa lượng tử của chuỗi" và đoạn Outlook: kiểm chứng số liệu

> Tiếp nối [bom-pha-omega-c-kiem-chung.md](bom-pha-omega-c-kiem-chung.md).
> Số liệu tính bằng [giai_so_quang_hoc_pho.py](giai_so_quang_hoc_pho.py) (cần numpy; chạy vài giây).
> Thế quang học dùng **dạng thể tích Koning–Delaroche (2003)**: $W_V(\Delta E)=w_1\Delta E^2/(\Delta E^2+w_2^2)$, với $w_1=12.195+0.0167A$ MeV và $w_2=73.55+0.0795A$ MeV. **Chưa tính** số hạng hấp thụ mặt, vốn lớn hơn ở năng lượng thấp; cần đối chiếu bảng KD03 đầy đủ.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $W(E)\propto(E-\varepsilon_F)^2$ gần mặt Fermi (chặn Pauli) | ✓ Đúng hiện tượng học |
| $W$ bão hòa với thang $E_c\approx80$–100 MeV | ✓ **Khớp số với tham số KD:** $w_2=$ **78.0 MeV** (A=56), **90.1 MeV** (A=208) |
| Vậy $\hbar\omega_c\approx80$ MeV "đến từ thế quang học" | **Chỉ là trùng thang, và có hai vấn đề.** (1) $W$ là hàm của **năng lượng nucleon** $E$, còn $\omega_c$ là **tần số** của thăng giáng trong bể. Đồng nhất hai đại lượng này là một giả thuyết, không phải suy dẫn. (2) **Quan trọng hơn:** ở năng lượng của mô hình (hạt ~6 MeV, xung ~4 MeV trên đáy), KD cho $W_V=0.05$–$0.21$ MeV, tức **$\gamma=2W/\hbar\omega_0=$ 0.011–0.042**. Con số này **nhỏ hơn 17–64 lần** so với $\gamma_x=0.70$ mà Part I dùng. Muốn đạt $\gamma=0.7$, nucleon phải ở **~47 MeV** trên mặt Fermi |
| "Cơ học cổ điển giả định chuỗi Zeno **vô hạn bước**, $\Delta E_{\text{Carnot}}\to0$" | **Sai với chính mô hình Part I.** Với $e=0$, mỗi đợt có **hữu hạn** bước (TB 5.12, phân bố hình học), mỗi bước TB **0.40 $\hbar\omega_0$**, không vô cùng bé. Chuỗi Zeno vô hạn chỉ xuất hiện khi $e>0$ |
| "**Không tồn tại** chuyển dời 0.4 $\hbar\omega_0$ giữa các trạng thái riêng của $H$" | **Sai.** Trong cửa sổ $E\in(2,6)$ (vùng rô-to quay ℓ ≈ 2–3) có 21 mức riêng. **Cả 20/20 khe liền kề đều < 0.5 $\hbar\omega_0$** (trung vị 0.21, lớn nhất 0.34) |
| Nên hệ "**bắt buộc**" nhảy các bước 1.5–2 $\hbar\omega_0$ | **Không bắt buộc. Đó là do chọn toán tử $\hat u$.** Với trọng số $\lvert\langle j\vert\hat u\vert i\rangle\rvert^2$, các chuyển dời đi xuống từ cửa sổ trên có TB **1.04 $\hbar\omega_0$**. **0%** trọng số rơi vào $\Delta E<0.5$. Tức là $\hat u$ **không nối** các mức gần nhau, theo quy tắc chọn của nó ($\hat p_x$ đổi $n$ ±1, còn $\cos\theta\,\hat p_\theta$ đổi $\ell$ ±1) |
| "Sai lệch lượng tử–cổ điển (1–2 bước lớn so với ~5 bước nhỏ) là **dự đoán trung tâm** của Part II" | **Chưa đủ vững để gọi là dự đoán.** Kết quả phụ thuộc vào việc chọn $\hat u$ làm toán tử ghép. Một toán tử ghép khác nối các mức gần nhau có thể cho các bước nhỏ. Trước khi gọi là dự đoán, cần chỉ ra **toán tử ghép do vật lý quyết định** và cho thấy kết quả ổn định khi đổi các chi tiết của nó |
| Đoạn Outlook: "Monte Carlo … (**300 trajectories**) … 1.2–1.6 jumps of ~2.05 ħω₀" | **Ghép số từ hai lần chạy khác nhau.** Số 1.2–1.6 bước và 2.05 $\hbar\omega_0$ đến từ lần chạy **có bơm**: 16 quỹ đạo × ~240 $T_0$. 300 quỹ đạo là lần chạy **chỉ làm nguội** (TB 2.08 bước, 1.70 $\hbar\omega_0$). Tại ⟨ℓ⟩ ≈ 2.4 thì số bước mỗi đợt là **1.57**, không phải cả khoảng 1.2–1.6 |
| Đoạn Outlook: $\hat u$ dấu "+", tính Hermite, nung nóng | ✓ |
| Đoạn Outlook: $\hbar\omega_c\gtrsim80$ MeV, $\langle p^2\rangle\approx0.79$–0.90 (tăng 60–80%) | ✓ (0.794 và 0.904, tức +59% và +81%) |
| "Part I đã chuẩn PRC và sẵn sàng xuất bản" | **Chưa kiểm được.** Không có file `paper_prc.tex` nào trong dự án. Bản LaTeX gần nhất đã kiểm ([bai-bao-latex-kiem-chung.md](bai-bao-latex-kiem-chung.md)) còn lỗi đơn vị ở Bảng I, ký tự tiếng Việt trong công thức, tùy chọn `prd`, email placeholder, v.v. |

## 1. Thế quang học (thể tích, KD03)

| A | ΔE (MeV) | 5 | 10 | 20 | 26.3 | 50 | 100 |
|---|---|---|---|---|---|---|---|
| 56 | $W_V$ (MeV) | 0.054 | 0.212 | 0.810 | 1.340 | 3.824 | 8.163 |
| 56 | $\gamma=2W/\hbar\omega_0$ | **0.011** | **0.042** | 0.162 | 0.268 | 0.765 | 1.633 |
| 208 | $W_V$ (MeV) | 0.048 | 0.191 | 0.736 | 1.231 | 3.690 | 8.649 |
| 208 | $\gamma$ | 0.010 | 0.038 | 0.147 | 0.246 | 0.738 | 1.730 |

$\gamma=0.70$ (tức $W=3.5$ MeV) cần $\Delta E=$ **47 MeV** (A=56) hoặc **48 MeV** (A=208).

## 2. Phổ và chuyển dời $\hat u_-$

| Chỉ số | Giá trị |
|---|---|
| Số mức riêng trong $E\in(2,6)$ | 21 |
| Khe liền kề: min / trung vị / max | 0.000 / 0.212 / 0.341 $\hbar\omega_0$ |
| Số khe < 0.5 $\hbar\omega_0$ | **20/20** |
| Năng lượng chuyển dời TB theo trọng số của $\hat u$ | **1.040 $\hbar\omega_0$** |
| Trọng số cho $\Delta E<0.5$ | **0.0%** |
| Trọng số cho $\Delta E>1.5$ | 3.1% |

## 3. Đoạn Outlook đã sửa (thay mục (ii) và (iii))

```latex
\item \textbf{Quantum discretization of the cascade (operator-dependent):} Combining a
Davies-type relaxation channel $\hat u_-$ with an incoherent saturating phase pump
$L_p=\sqrt{\Gamma_0}\,e^{i\theta}f(\hat\ell)$ yields a stationary quantum limit cycle with
geometrically distributed relaxation bursts. In Monte Carlo wave-function simulations
(16 trajectories of $\approx240\,T_0$), at $\langle\ell\rangle\approx2.4$ the bursts contain
$1.57$ jumps on average, each removing $\approx2.05\,\hbar\omega_0$, compared with $5.12$
classical impacts of $\approx0.40\,\hbar\omega_0$. Although the spectrum of $H$ contains
many gaps below $0.5\,\hbar\omega_0$, the selection rules of $\hat u$ give them zero
weight; the discrepancy is therefore a property of the chosen coupling operator and
becomes a prediction only once that operator is fixed by the underlying physics.
\item \textbf{Environmental cutoff:} Retaining $\ge90\%$ of the classical friction at
$\Omega=2.63\,\omega_0$ requires $\hbar\omega_c\gtrsim80$~MeV, numerically close to the
saturation scale $w_2\simeq78$--$90$~MeV of the Koning--Delaroche volume absorption.
However, at the nucleon energies relevant here (5--10~MeV above $\varepsilon_F$) that
absorption corresponds to $\gamma\approx0.01$--$0.04\,\omega_0$, one to two orders of
magnitude below the $\gamma_x=0.70\,\omega_0$ used in this work.
```

## Câu hỏi còn mở

- **Toán tử ghép vật lý.** Toán tử nào (thay cho $\hat u$) được suy ra từ vật lý của liên kết và va chạm? Cần có câu trả lời này thì "1–2 bước lớn" mới là dự đoán của Part II.
- **$\gamma_x=0.70$ ở Part I** lớn hơn hấp thụ thể tích KD ở năng lượng của mô hình 17–64 lần. Cần kiểm thêm hấp thụ mặt KD03 ở 5–10 MeV. Nếu vẫn nhỏ, $\gamma_x$ cần được hạ, hoặc được biện minh bằng một cơ chế khác.
- **File bản thảo.** Nếu `paper_prc.tex` đã được sửa, hãy đặt nó vào thư mục `hoc-thuyet-e/` để mình kiểm lại toàn bộ trước khi nộp.
