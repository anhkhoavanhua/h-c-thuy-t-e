# Page–Wootters lượng tử có môi trường và đồng hồ đếm vòng cho $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** bản thảo [HocThuyetE (5).pdf](../../HocThuyetE%20(5).pdf). Bài này trả lời hai câu hỏi mở cuối [dong-ho-photon-kiem-chung.md](dong-ho-photon-kiem-chung.md).
> Số liệu tính bằng [giai_so_pw_luong_tu_mo.py](giai_so_pw_luong_tu_mo.py). Lệnh `python3 giai_so_pw_luong_tu_mo.py 1` chạy ≈ 30 s, lệnh `… 2` chạy ≈ 30 s trên 8 lõi. Script chỉ import các module cũ, **không sửa** file nào.

## Mô hình

**Hệ S** (như phần Outlook của bản thảo, cơ sở thu nhỏ để chạy nhanh):
- Rô-to $\ell\in[-5,5]$ ⊗ dao động tử $n\le4$, tức 55 chiều.
- $\hat H$ với ghép $-g\,x\sin\theta$ (dấu đúng), $g=0.1$.

**Môi trường E** nhận hai loại lượng tử:
- **Làm nguội Carnot:** $L_c=\sqrt{0.5}\,\hat u_-$ (Davies, $T=0$).
- **Bơm pha:** $L_p=\sqrt{0.2}\,e^{i\theta}f(\hat\ell)$, với $\ell_0=4$.

**Đồng hồ C** đếm nhịp. **Trạng thái lịch sử** là $\lvert\Psi\rangle\rangle=\sum_k\lvert k\rangle_C\otimes\lvert\psi_{SE}(k)\rangle$. Mỗi nhịp, môi trường E ghi lại một trong ba kết quả: không nhảy, nhảy Carnot, hoặc nhảy bơm.

Đơn vị: $\hbar=m=\omega_0=1$, $T_0=2\pi$.

## Câu 1. Phiên bản lượng tử: độ tinh khiết của hạt khi tiêu tán

### 1a. Kiểm tường minh trạng thái lịch sử C⊗S⊗E

| Nhịp $k$ | Số nhánh của E | $\lVert\rho_S^{\text{lịch sử}}-\rho_S^{\text{Lindblad}}\rVert$ | $\operatorname{Tr}\rho_S^2$ |
|---|---|---|---|
| 1 | 3 | 7.5·10⁻¹⁹ | 0.9802 |
| 2 | 9 | 1.1·10⁻¹⁶ | 0.9609 |
| 3 | 27 | 4.6·10⁻¹⁷ | 0.9422 |

- Lấy vết trên E của trạng thái điều kiện $\langle k\vert\Psi\rangle\rangle$ cho **đúng nghiệm Lindblad**, sai khác chỉ ở mức làm tròn máy tính.
- Toàn bộ C⊗S⊗E vẫn thuần, nên "vũ trụ tĩnh" được giữ nguyên.
- S mất độ tinh khiết **chỉ vì** nó vướng víu với bản ghi của E.

### 1b. Đồng hồ đếm xung lý tưởng: $\rho_S(t)=\operatorname{Tr}_E\langle t\vert\Psi\rangle\rangle$

Xuất phát từ $\lvert\ell=0,n=0\rangle$:

| $t/T_0$ | $\langle\ell\rangle$ | $\operatorname{Tr}\rho^2$ (rô-to + hạt) | $S_{vN}$ | $\operatorname{Tr}\rho^2$ (riêng hạt) | $S_{vN}$ (hạt) |
|---|---|---|---|---|---|
| 0 | 0.00 | 1.000 | 0.000 | 1.000 | 0.000 |
| 0.1 | 0.12 | 0.790 | 0.392 | 0.998 | 0.007 |
| 0.5 | 0.30 | 0.486 | 1.002 | 0.990 | 0.032 |
| 1 | 0.63 | 0.292 | 1.440 | 0.981 | 0.055 |
| 4 | 1.60 | **0.189** | 1.769 | **0.963** | 0.096 |
| 16 (dừng) | 2.31 | 0.246 | 1.599 | 0.972 | 0.076 |

**Kết quả:**
1. **Độ tinh khiết giảm trơn, không có bước sụt.** Tốc độ giảm lớn nhất là 0.40 mỗi đơn vị thời gian, ngay tại $t=0$. Con số này đúng bằng $2\Gamma_0$, tức phần lớn độ tinh khiết mất đi do **kênh bơm**, không phải do kênh Carnot.
2. **Riêng hạt gần như vẫn thuần.** $\operatorname{Tr}\rho_x^2\ge0.96$ suốt quá trình. Lý do: ghép $g=0.1$ yếu, và làm nguội Carnot kéo hạt về gần trạng thái cơ bản.
3. **Không đơn điệu.** Độ tinh khiết chạm đáy ở khoảng $4T_0$ rồi hồi lên ở trạng thái dừng. "Mất dần mãi" không đúng. Nó hội tụ về một trạng thái dừng hỗn hợp.

### 1c. Điều kiện hóa cả bản ghi của môi trường (quỹ đạo lượng tử)

Mô phỏng một quỹ đạo dài 400 đơn vị thời gian (≈ 64 $T_0$), có 41 bước nhảy Carnot.

| Đại lượng | Giá trị |
|---|---|
| $\operatorname{Tr}\rho^2$ (rô-to + hạt) | **1 tại mọi thời điểm** |
| $\operatorname{Tr}\rho_x^2$ giữa các bước nhảy | trung bình 0.977 (min 0.890) |
| … ngay **trước** bước Carnot | trung bình 0.980 |
| … ngay **sau** bước Carnot | trung bình **0.994** |

**Điểm cần sửa trong hình ảnh "mỗi lần tiêu tán hạt mất tinh khiết":**
- Nếu đọc được bản ghi của E, một lần phát lượng tử Carnot **tăng** độ tinh khiết của hạt (0.980 → 0.994). Lượng tử Carnot mang năng lượng đi, hạt rơi về gần trạng thái cơ bản, và vướng víu rô-to–hạt bị cắt.
- Mất độ tinh khiết là tính chất của **tập hợp** (không biết E đã ghi gì). Nó không phải hệ quả của từng biến cố tiêu tán.
- Câu đúng cho bản thảo: *"trạng thái điều kiện của hạt, lấy trung bình trên các bản ghi của môi trường, mất độ tinh khiết một cách trơn tru; dọc một bản ghi xác định, hạt vẫn thuần và mỗi lần phát lượng tử Carnot làm nó thuần hơn."*

## Câu 2. Photon đồng hồ đếm vòng bằng cách nào?

### 2a. Mô-men quỹ đạo (OAM)

Trạng thái góc của photon OAM là $\lvert\varphi\rangle=\sum_{\ell=-L}^{L}e^{-i\ell\varphi}\lvert\ell\rangle/\sqrt{2L+1}$.

| Dải OAM | $\lvert\langle\varphi\vert\varphi+2\pi\rangle\rvert^2$ | Số thời điểm trực giao trong một vòng |
|---|---|---|
| $\ell\in[-1,1]$ | **1.000000** | 3 |
| $\ell\in[-5,5]$ | **1.000000** | 11 |
| $\ell\in[-20,20]$ | **1.000000** | 41 |

- OAM làm **mịn** số đọc trong một vòng: $2L+1$ thời điểm phân biệt được, so với 2 thời điểm của photon phân cực.
- Nhưng $\ell$ nguyên thì góc liên hợp **luôn tuần hoàn 2π**. Vì vậy OAM **không đếm vòng được**, bất kể $L$ lớn đến đâu.
- Về cấu trúc, đây chính là rô-to θ của $\mathcal E_{\text{dual}}$. Thay photon phân cực bằng photon OAM không giải quyết được vấn đề.

### 2b. Chuỗi photon liên tiếp (bộ đếm xung)

Đồng hồ phát một photon mỗi $T_{\text{rep}}$. Số đọc là **số photon đã phát**, một số nguyên tăng đơn điệu nên đếm vòng được. Câu hỏi còn lại: $T_{\text{rep}}$ phải nhỏ đến đâu để phân biệt được Hiện tại (lần gãy) với Tương lai (lần tái khóa)?

Mô phỏng Langevin với tham số Bảng I (1 MeV): 8 quỹ đạo × 40 $T_0$, gồm 31 đợt và 179 cặp gãy→khóa.
- Thời lượng đợt: trung bình 1.02 $T_0$, ngắn nhất 0.16 $T_0$.
- Khoảng từ một lần gãy đến lần khóa kế tiếp: trung bình 0.15 $T_0$, ngắn nhất 0.0003 $T_0$.

| $T_{\text{rep}}$ ($T_0$) | Xung / vòng rô-to | Đợt không xếp được thứ tự (gãy mở đợt và khóa kết thúc rơi cùng một xung) | Cặp gãy→khóa kề nhau không phân biệt |
|---|---|---|---|
| 0.01 | 244 | 0.0% | 36.2% |
| 0.05 | 49 | 0.0% | 38.5% |
| 0.10 | 24 | **0.0%** | 38.8% |
| 0.20 | 12 | 0.9% | 40.0% |
| 0.50 | 4.9 | 17.9% | 70.3% |
| 1.00 | 2.4 | 36.6% | 84.9% |

**Đọc bảng:**
1. **Ba pha của một đợt** (Quá khứ → Hiện tại → Tương lai) được xếp đúng 100% khi $T_{\text{rep}}\le0.1\,T_0$, tức khoảng **≥ 25 photon mỗi vòng rô-to**. Ở 1 photon mỗi vòng ($T_{\text{rep}}=T_0$), 37% số đợt bị mất thứ tự.
2. **Khoảng 36% cặp gãy→khóa kề nhau** cách nhau dưới 0.01 $T_0$. Đây là hiện tượng rung gãy/khóa liên tiếp bên trong một đợt, và không chuỗi xung thực tế nào tách được. Vì vậy học thuyết nên định nghĩa "Hiện tại" và "Tương lai" theo **đợt**, không theo từng lần gãy riêng lẻ. Đây cũng là cách đã làm trong [dong-ho-photon-kiem-chung.md](dong-ho-photon-kiem-chung.md).

### Kết luận câu 2

| Phương án | Đếm vòng? | Ghi chú |
|---|---|---|
| Một photon phân cực | Không | Chu kỳ π, 2 trạng thái |
| Một photon OAM | Không | Chu kỳ 2π, $2L+1$ trạng thái. Tương đương rô-to θ |
| Chuỗi photon (bộ đếm) | **Có** | Cần ≳ 25 photon mỗi vòng ở 1 MeV. Số đọc là số nguyên đơn điệu |
| OAM + bộ đếm vòng | **Có** | Pha trong vòng lấy từ OAM, số vòng lấy từ bộ đếm. Tương đương Φ tích lũy của $\mathcal E_{\text{dual}}$ |

Để giữ hình ảnh "photon làm đồng hồ", cần thay "**một** photon phân cực" bằng "**một chuỗi** photon". Mỗi photon là một tích tắc, và thời gian là số tích tắc đã đếm được. Mũi tên thời gian vẫn đến từ bản ghi tiêu tán trong môi trường (Câu 1), không đến từ bản thân đồng hồ.

## Ghi chú về hai chỗ sửa LaTeX
- **V.B.4:** khoảng "0.35–0.45 cho 0f/0g" nên sửa thành **0.38–0.45**. Với cả dải $r_{\text{rms}}=3.8$–4.8 fm, khoảng đúng là 0.22–0.45.
- **V.C.2:** đoạn đề xuất nhất quán với V.B.4. Cần bỏ câu cũ "…dominates".
