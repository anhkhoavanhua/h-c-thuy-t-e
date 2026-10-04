# Tách bước nhảy lượng tử MCWF (V.C.1): rô-to $\ell\to\ell-1$ hay hạt $n\to n-1$?

> **Cơ sở:** mục V.C.1 của [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf). Nối tiếp câu hỏi mở cuối [schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md).
> Số liệu tính bằng [giai_so_mcwf_tach_buoc.py](giai_so_mcwf_tach_buoc.py). Lệnh `python3 giai_so_mcwf_tach_buoc.py` chạy ≈ 17 s trên 8 lõi.
> **Mô hình y hệt V.C.1** (cấu hình "bão hòa" của [giai_so_bom_pha.py](giai_so_bom_pha.py)):
> - Rô-to $\ell\in[-8,8]$ ⊗ dao động tử $n\le6$, $H$ có số hạng $-g\,x\sin\theta$.
> - Bơm $L_p=\sqrt{0.2}\,e^{i\theta}f(\hat\ell)$, $\ell_0=4$.
> - Làm nguội Davies $L_c=\sqrt{0.5}\,\hat u_-$.
> - 16 quỹ đạo × 1500 đơn vị thời gian (≈ 240 $T_0$), bỏ 1/5 đầu. Tổng cộng 1582 bước làm nguội.

## 1. Năng lượng mất đi ở mỗi bước, theo từng khối

| Đại lượng (trước − sau) | Trung bình | Trung vị | 10–90% |
|---|---|---|---|
| $\Delta E$ toàn phần | **+2.051** $\hbar\omega_0$ | +2.268 | +0.79…+3.17 |
| Rô-to $\langle\ell^2/2I-\Delta E\cos\theta\rangle$ | **+2.044** | +2.262 | +0.78…+3.17 |
| Hạt $\langle n+\tfrac12\rangle$ | +0.008 | +0.005 | +0.001…+0.018 |
| Ghép $\langle-gx\sin\theta\rangle$ | −0.000 | | |
| $\Delta\langle\lvert\ell\rvert\rangle$ | +0.62 | +0.77 | |
| $\Delta\langle n\rangle$ | +0.008 | | |

- Con số **2.051** $\hbar\omega_0$ tái lập đúng giá trị "≈ 2.05 ħω0" của V.C.1. Mô hình được dựng lại chính xác.
- **99.7% năng lượng mất đi lấy từ rô-to.** Hạt chỉ đóng góp 0.4%.

## 2. Phân rã theo chuyển dời riêng $j\to i$

Mỗi trạng thái riêng của $H$ được gắn nhãn $(|\ell|,n)$ theo thành phần trội. Trọng số của chuyển dời $j\to i$ là $|A_{ij}c_j|^2$, trong đó $A$ là $L_c$ viết trong cơ sở riêng.

| Lớp chuyển dời | Tỉ phần số bước | Khe TB ($\hbar\omega_0$) | Tỉ phần năng lượng |
|---|---|---|---|
| **Rô-to $\ell\to\ell-1$** | **91.8%** | 2.84 | **93.5%** |
| Rô-to $\ell\to\ell-2$ | 2.9% | 3.81 | 4.0% |
| **Hạt $n\to n-1$** | **4.1%** | **1.00** | 1.5% |
| Trộn / khác | 1.2% | 2.41 | 1.0% |

**Phổ khe năng lượng (có trọng số):**

| Khe ($\hbar\omega_0$) | 0–0.5 | 0.5–1.5 | 1.5–2.5 | 2.5–3.5 | 3.5–5 | > 5 |
|---|---|---|---|---|---|---|
| Tỉ phần | **0.0%** | 10.1% | 43.3% | 43.7% | 2.4% | 0.5% |

**Các chuyển dời mạnh nhất.** Mỗi cặp ±ℓ là hai trạng thái riêng suy biến, nên xuất hiện thành hai dòng giống nhau:

| Chuyển dời $(\lvert\ell\rvert,n)$ | Khe | Trọng số |
|---|---|---|
| (4, 0) → (3, 0) | 3.487 | 2 × 21.2% |
| (3, 0) → (2, 0) | 2.447–2.467 | 2 × 19.5% |
| (2, 0) → (1, 0) | 1.21–1.75 | 9.5% |
| (1, 1) → (1, 0): bước của hạt | 0.991 | 1.9% |

## 3. Kết luận

1. **Làm nguội lượng tử trong V.C.1 thực chất là hãm rô-to theo bậc thang $\ell\to\ell-1$.** Bước này chiếm 92% số bước và 94% năng lượng. Các khe là $(2\ell-1)/2I$: ℓ 4→3 = 3.49, ℓ 3→2 = 2.45, khớp phổ rô-to tự do với $I=1$.
2. **Lượng tử $\hbar\omega_0$ của hạt gần như vắng mặt:** chỉ 4% số bước và 1.5% năng lượng. Hạt luôn ở gần trạng thái cơ bản ($\langle n\rangle\lesssim0.1$), vì bơm chỉ cấp năng lượng cho rô-to và ghép $g=0.1$ yếu.
3. **Không có bước nào dưới 0.5 $\hbar\omega_0$** (0.0%). Va chạm Carnot 0.4 $\hbar\omega_0$ của mô hình cổ điển không có đối ứng lượng tử. Nó chỉ là trung bình thống kê của các bước 1–3.5 $\hbar\omega_0$.
4. **So với Schrödinger trên đa tạp ràng buộc:** bước ℓ 3→2 ở đó là 2.32 (quán tính $\bar M=1.054$), ở đây là 2.45 (quán tính $I=1$). Lệch 5% đúng bằng phần khối lượng $m\kappa^2\cos^2\theta$ mà hạt góp vào rô-to khi bị ràng buộc. V.C.1 dùng ràng buộc mềm ($g=0.1$), nên rô-to gần như tự do.
5. **Hệ quả cho mục V.A (exciton):**
   - Bản thảo ánh xạ mỗi va chạm 4.0–4.6 MeV thành một bước tạo exciton $\Delta n=+2$, và cần trung vị **3 va chạm** để vượt $S_n=8$ MeV.
   - Trong mô hình lượng tử, mỗi bước mang **20–35 MeV**, vượt $S_n$ ngay từ **bước đầu tiên**.
   - Tính theo đợt: lượng tử 1.57 bước × 2.05 = **3.2 $\hbar\omega_0$ (32 MeV)**, còn cổ điển 5.12 × 0.43 ≈ 2.2 $\hbar\omega_0$ (22 MeV).
   - Phần Ericson của bản thảo vì vậy chỉ đúng cho giới hạn cổ điển. Nên ghi rõ điều này.

## 4. Câu sửa đề xuất cho bản thảo (tiếng Anh)

Mục V.B.2, thay câu "would be resolved into discrete operator transitions in multiples of ħω0" bằng:

> In the quantum model of Sec. V.C.1, 92% of the relaxation jumps are rotor-ladder transitions $\ell\to\ell-1$ (gaps $(2\ell-1)\hbar^2/2I\approx2.5$–$3.5\,\hbar\omega_0$), carrying 99.7% of the dissipated energy; particle transitions $n\to n-1$ (gap $\hbar\omega_0$) account for only 4% of the jumps, and no transition below $0.5\,\hbar\omega_0$ occurs. The classical sub-quantum losses ($\approx0.4\,\hbar\omega_0$) are therefore ensemble averages of rotor-ladder quanta.

## Câu hỏi còn mở
- Nếu muốn lượng tử $\hbar\omega_0$ của hạt thật sự tham gia, cần ghép mạnh hơn hoặc bơm cả vào hạt. Có thể quét $g$ từ 0.1 lên 1 để tìm điểm mà bước $n\to n-1$ bắt đầu chiếm phần đáng kể.
