# Part II: động lực học lượng tử mở của $\mathcal E_{\text{dual}}$, giải số và kiểm chứng kế hoạch

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf) (Part I), mô hình V.C.1. Dùng cùng cấu hình với [paper_part3.tex](paper_part3.tex).
> **Số liệu:** [giai_so_part2_master.py](giai_so_part2_master.py), gồm bốn module:
> - `M1`: ≈ 3 phút;
> - `M2`: ≈ 4 phút trên 8 lõi;
> - `M3`: vài giây;
> - `M4`: vài giây.
>
> **Mô hình:** $H_0=p_\theta^2/2I-\Delta E\cos\theta+(a^\dagger a+\tfrac12)-g\,x\sin\theta$ (Bảng I, dấu đúng), với ba kênh:
> - $L_p=\sqrt{0.2}\,e^{i\theta}f(\hat\ell)$, $f=\sqrt{\max(0,1-\ell^2/16)}$;
> - $L_c=\sqrt{0.5}\,\hat u_-$ (Davies);
> - $L_x=\sqrt{0.7}\,a$ là kênh **mới** mà kế hoạch đề xuất. Mình chạy cả có lẫn không có kênh này.

## 0. Kiểm kế hoạch: những chỗ cần sửa trước khi viết bài

| Điểm trong kế hoạch | Kết quả số | Cần sửa |
|---|---|---|
| $f(\hat\ell)=\Theta(\ell_0-\hat\ell)$ hoặc $\sqrt{\max(0,1-\hat\ell/\ell_0)}$ | Mô hình đã kiểm chứng (V.C.1, Part III) dùng $\sqrt{\max(0,1-\ell^2/\ell_0^2)}$ | Ghi đúng dạng $\ell^2$. Dạng tuyến tính là một mô hình khác, chưa kiểm |
| Module 1: trạng thái dừng có $\langle\ell\rangle\approx2.26$, $\operatorname{Tr}\rho^2\approx0.241$ | Đó là số của **mô hình va chạm** $\tau=T_0/25$ (Part III). Nghiệm **chính xác** của phương trình Lindblad là **$\langle\ell\rangle=2.325$, $\operatorname{Tr}\rho^2=0.248$** | Part II dùng 2.325 / 0.248. Độ lệch 3% là sai số $O(\tau)$ đã đo ở Part III |
| $\hat u=\kappa\dot{\hat\theta}\cos\hat\theta-\dot{\hat x}$ | Toán tử dùng trong số liệu là $\hat u=(\kappa/I)\tfrac12\{\cos\theta,p_\theta\}-p_x$ (dạng Hermit hóa, $m=1$) | Ghi dạng đối xứng hóa |
| Purity riêng hạt tăng 0.980 → 0.994 sau bước Carnot | ✓ **Xác nhận trên cơ sở đầy đủ V.C.1**: 0.9803 → 0.9943, tăng ở **99%** số bước | Không cần sửa |
| **$R_\delta\approx3.3$: "bó sóng ràng buộc là trạng thái con trỏ **duy nhất** nhân bản thông tin vào môi trường"** | ✗ **Không được chứng minh.** Con số 3.3 đến từ mô hình đồ chơi với bản ghi lý tưởng của các mức ℓ, không từ động lực học của bó sóng ràng buộc. Ở mô hình va chạm thật (M4): **không có bình nguyên**, $R_\delta\le1.75$ | Tách hai khẳng định: (a) *sieve*: bó sóng ràng buộc **bền nhất** (đã chứng minh, Part III bảng IV); (b) *Darwinism*: chỉ có ở mô hình bản ghi lý tưởng. Bỏ chữ "duy nhất" |
| HEOM ở Module 3 | Mình **không** chạy HEOM. Với dao động tử tuyến tính, nghiệm chuyển động Brown lượng tử chính xác (Feynman–Vernon / Ford–Lewis–O'Connell) **là** chuẩn mà HEOM hội tụ về | Ghi "exact QBM solution (the benchmark HEOM converges to)", không ghi "HEOM" |

## M1. Trạng thái dừng chính xác (siêu toán tử Liouville)

Giải trực tiếp $\mathcal L\,\mathrm{vec}(\rho)=0$ với điều kiện $\operatorname{Tr}\rho=1$:

| Cơ sở | Kênh | ⟨ℓ⟩ | ⟨H⟩ | ⟨n⟩ | Tr ρ² | S (nat) | Nhảy Carnot / $T_0$ | Nhảy bơm / $T_0$ | Nhảy $L_x$ / $T_0$ | Khe Liouvillian |
|---|---|---|---|---|---|---|---|---|---|---|
| ℓ≤5, n≤4 (3025²) | V.C.1 | **2.325** | 4.238 | 0.0143 | **0.2481** | 1.591 | 0.500 | 0.676 | — | **0.0488** ω₀ |
| ℓ≤5, n≤4 | + $L_x$ | 2.383 | 4.319 | 0.0040 | 0.2570 | 1.539 | 0.485 | 0.662 | 0.0174 | 0.0477 ω₀ |
| ℓ≤6, n≤5 (5929²) | V.C.1 | 2.326 | 4.238 | 0.0143 | 0.2482 | 1.590 | 0.501 | 0.676 | — | |
| ℓ≤6, n≤5 | + $L_x$ | 2.384 | 4.319 | 0.0040 | 0.2571 | 1.538 | 0.485 | 0.662 | 0.0174 | |

- **Hội tụ theo cỡ cơ sở:** sai khác ≤ 0.001. Trị riêng nhỏ nhất của ρ ~ 10⁻¹⁷, tức ρ dương, vật lý.
- **Khe Liouvillian 0.049 ω₀:** thời gian hồi phục chậm nhất là $1/0.049=20.5$ đơn vị thời gian ≈ **3.3 $T_0$**. Cùng bậc với thời gian chờ giữa các đợt (2.8 $T_0$, mục M2).
- **Tác dụng của $L_x$ (γx = 0.7) rất nhỏ:** kênh này chỉ nhảy 0.017 lần mỗi $T_0$, vì hạt hầu như đã ở trạng thái cơ bản ($\langle n\rangle$ giảm từ 0.014 xuống 0.004). ⟨ℓ⟩ đổi +0.06. Ma sát quang học không thay đổi bức tranh, vì rô-to mới là nơi tiêu tán.

## M2. Quỹ đạo lượng tử (MCWF): 100 quỹ đạo × 1500 đơn vị thời gian ≈ 239 $T_0$

| Đại lượng | V.C.1 | V.C.1 + $L_x$ |
|---|---|---|
| Tốc độ nhảy / $T_0$ | Carnot 0.502, bơm 0.682 | Carnot 0.490, bơm 0.666, $L_x$ 0.017 |
| Năng lượng mỗi bước Carnot | TB **2.013** $\hbar\omega_0$, trung vị 2.243 | TB 2.102, trung vị 2.297 |
| Phổ ΔE: 0–0.5 / 0.5–1.5 / 1.5–2.5 / 2.5–3.5 / > 3.5 | 2.6% / 14.5% / **48.4% / 28.4%** / 0.1% | 2.4% / 11.9% / 48.9% / 31.6% / 0.1% |
| **ℓ → ℓ−1** | **91.5% số bước / 93.2% năng lượng** (khe 2.80) | 93.1% / 93.9% (khe 2.83) |
| ℓ → ℓ−2 | 3.0% / 4.2% (khe 3.77) | 3.0% / 4.1% |
| **n → n−1** | **4.2% / 1.5%** (khe 1.00) | 2.8% / 1.0% |
| Khác | 1.2% / 1.1% | 1.2% / 1.0% |
| Đợt Carnot (khoảng cách < $T_0$) | 5984 đợt, **1.57** bước/đợt | 5949 đợt, 1.54 |
| Phân bố số bước mỗi đợt: 1 / 2 / 3 / 4 | 61% / 27% / 9% / 2% | 62% / 26% / 9% / 2% |
| Chờ giữa đợt | **2.83 $T_0$**, CV 0.73 | 2.86 $T_0$, CV 0.74 |
| Purity riêng hạt, trước → sau bước Carnot | **0.9803 → 0.9943** (tăng ở 99%) | 0.9939 → 0.9977 (100%) |

**Ghi chú:**
- "Phổ ΔE" ở đây là độ giảm $\langle H\rangle$ của trạng thái điều kiện. 2.6% số bước có ΔE < 0.5. Nhưng phân rã theo chuyển dời riêng **không có khe nào < 0.5** $\hbar\omega_0$: các ΔE nhỏ là do trạng thái trước bước nhảy là chồng chập nhiều mức.
- Phân bố số bước mỗi đợt **gần hình học**: phân bố hình học cùng trung bình 1.57 cho 64% / 23% / 8% / 3%, so với đo 61% / 27% / 9% / 2%. Đợt lượng tử ngắn, khác hẳn 5.12 va chạm cổ điển, nhưng giữ dạng hình học như Part I.
- Thời gian chờ có CV = 0.73 < 1. Đó **không phải** phân bố mũ thuần như pha Quá khứ cổ điển (CV ≈ 1–1.1).

**So với kế hoạch:** 92% / 99.7% / 4% / 1.51–1.57 / 2.84 $T_0$ đều được xác nhận. Lưu ý **99.7%** là tỉ phần năng lượng lấy từ khối rô-to (từ [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md)). Đó không phải tỉ phần của riêng lớp ℓ → ℓ−1, vốn là 93.2%. Viết bài cần tách hai con số này.

## M3. Phi Markov tại $k_BT=0.1\,\hbar\omega_0$ (1 MeV), γ = 0.7

| Phương trình | ⟨x²⟩ | ⟨p²⟩ | ⟨x²⟩⟨p²⟩ | Trị riêng nhỏ nhất của ρ |
|---|---|---|---|---|
| Lindblad $L=\sqrt\gamma\,a$ (+ $n_{\text{th}}=4.5\cdot10^{-5}$) | 0.5000 | 0.5000 | 0.2500 | 0 (thuần, trạng thái cơ bản) |
| **Caldeira–Leggett** (Markov nhiệt độ cao) | 0.1000 | 0.1000 | **0.0100 < 1/4** | **−1.111** |
| QBM chính xác, Drude $\hbar\omega_c=4.5$ MeV | 0.4726 | 0.5612 | 0.2652 | 0.015 |
| QBM chính xác, $\hbar\omega_c=50$ MeV | **0.4349** | 0.7942 | 0.3454 | 0.074 |
| QBM chính xác, $\hbar\omega_c=100$ MeV | **0.4288** | 0.9037 | 0.3875 | 0.097 |

- **Caldeira–Leggett:** trạng thái dừng giải tích là $\langle x^2\rangle=\langle p^2\rangle=k_BT$. Tích 0.01 vi phạm bất đẳng thức Robertson–Schrödinger. Trạng thái Gauss với ν = 0.1 có trị riêng $\lambda_1=\frac{1}{\nu+1/2}\cdot\frac{\nu-1/2}{\nu+1/2}=-1.111$. Đây đúng là con số "−1.11" trong Part I, nay đã có công thức.
- **Lindblad $\sqrt\gamma\,a$** giữ tính dương nhưng **bỏ sót** hiệu ứng liên kết mạnh: ⟨x²⟩ thật bị ép xuống 0.43, ⟨p²⟩ bị nâng lên 0.79–0.90.
- **QBM chính xác** (nghiệm Feynman–Vernon cho hệ tuyến tính) là chuẩn phi Markov. ⟨x²⟩ = 0.43–0.47 khớp V.C.3 của Part I.

## M4. Darwinism động lực trong mô hình va chạm thật

**Thiết lập:**
- Qubit tham chiếu R vướng với hai nhánh của S: $\tfrac1{\sqrt2}(\lvert0\rangle_R\lvert A\rangle+\lvert1\rangle_R\lvert B\rangle)$.
- 7 bản ghi 3 mức (không / Carnot / bơm), S 36 chiều.
- Đo $I(R{:}F)$ theo số bản ghi trong mảnh F. Darwinism đòi $I\approx1$ bit ngay ở mảnh nhỏ (bình nguyên).

| Cặp nhánh | τ | $I(R{:}F)$ (bit) với $\lvert F\rvert$ = 1, 2, 3, 4, 5, 6, 7 | $R_{0.1}$ |
|---|---|---|---|
| **C**: bó sóng ràng buộc θ₀ = 0 / π | $T_0/25$ (0.28 $T_0$) | 0.12, 0.25, 0.40, 0.53, 0.64, 0.72, **0.78** | < 1 |
| **E**: trạng thái cơ bản / mức \|ℓ\| = 3 | $T_0/25$ | 0.06, 0.11, 0.17, 0.24, 0.31, 0.38, **0.44** | < 1 |
| **C** | $T_0/5$ (1.4 $T_0$) | 0.47, 0.73, 0.89, **1.01**, 1.11, 1.18, 1.28 | 7/4 = **1.75** |
| **E** | $T_0/5$ | 0.21, 0.43, 0.63, 0.81, **0.95**, 1.07, 1.16 | 7/5 = 1.4 |

**Kết quả:**
1. Môi trường **ghi nhận nhánh C hiệu quả hơn nhánh E**: 0.78 so với 0.44 bit sau 0.28 $T_0$. Kênh bơm $e^{i\theta}$ mang pha của rô-to, nên bản ghi bơm phân biệt được θ₀ = 0 với θ₀ = π. Vậy **góc của rô-to là thứ được ghi lại**, phù hợp với hình ảnh quỹ đạo cổ điển.
2. **Không có bình nguyên.** Thông tin tích lũy dần theo từng bản ghi, không bão hòa ở mảnh nhỏ. Độ dư thừa thực tế $R_{0.1}\le1.75$, **thấp xa so với 3.3** của mô hình bản ghi lý tưởng.
3. **Kết luận đúng cho bài báo:** bó sóng ràng buộc là **trạng thái con trỏ** (bền nhất theo sieve) và **được ghi nhận nhiều nhất**. Nhưng **tính khách quan kiểu Darwinism** (nhiều mảnh độc lập cùng biết) **chưa đạt** trong thang thời gian ≤ 1.4 $T_0$ với 7 bản ghi. Muốn khẳng định, cần nhiều bản ghi hơn: cấu trúc MPS hoặc lấy mẫu nhánh, là hướng phát triển tiếp.

## Bảng số liệu chuẩn để đưa vào bài Part II

| Mục | Giá trị chuẩn | Nguồn |
|---|---|---|
| Trạng thái dừng (Lindblad chính xác) | ⟨ℓ⟩ = 2.325, Tr ρ² = 0.248, S = 1.59 nat, khe 0.049 ω₀ (3.3 $T_0$) | M1 |
| Bước nhảy | ℓ → ℓ−1: 91.5% số bước / 93.2% năng lượng; khối rô-to 99.7% năng lượng; n → n−1: 4.2% | M2, [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md) |
| Năng lượng mỗi bước | 2.01 $\hbar\omega_0$ (TB), khe ℓ 3→2 = 2.32 trên đa tạp, 2.45 không ràng buộc | M2, [schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md) |
| Đợt | 1.57 bước/đợt, chờ 2.83 $T_0$ (CV 0.73) | M2 |
| Purity riêng hạt | 0.980 → 0.994 (99% số bước) | M2 |
| Con trỏ (sieve, kênh Davies thật) | C 0.942 / E 0.710 tại 0.1 $T_0$; tốc độ nhảy 0.023 vs 0.182 | [nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md) mục D |
| Darwinism | Lý tưởng: R = 3.3. Động lực thật: không bình nguyên, R ≤ 1.75 | M4 |
| Phi Markov | CL: tích 0.010, λ = −1.111; QBM chính xác ⟨x²⟩ = 0.43–0.47; Lindblad 0.5 | M3 |

## Câu hỏi còn mở
- Có muốn mình viết **LaTeX Part II** (revtex, cùng kiểu Part I/III) từ bảng số liệu chuẩn trên không? Bài sẽ dùng các con số đã sửa: 2.325 thay 2.26, không ghi "duy nhất", không ghi "HEOM".
- Kênh $L_x$ gần như không đổi kết quả. Nên giữ trong Part II như một kiểm tra độ bền, hay bỏ cho gọn?
