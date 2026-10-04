# Bài toán người quan sát (observer / measurement problem) cho học thuyết $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I và mô hình lượng tử V.C.1.
> Số liệu tính bằng [giai_so_nguoi_quan_sat.py](giai_so_nguoi_quan_sat.py):
> - `A`: ≈ 3.5 phút;
> - `B`, `B2`: ≈ 20 s mỗi phần.
>
> Script chỉ import `basis()` của [giai_so_u_hermite_drude.py](giai_so_u_hermite_drude.py), không sửa file nào.

## Đặt bài toán

**Bài toán người quan sát** gồm ba câu hỏi:
1. **Cơ sở ưu tiên:** vì sao ta thấy trạng thái này mà không thấy trạng thái kia?
2. **Tính khách quan:** vì sao nhiều người quan sát độc lập cùng thấy một kết quả?
3. **Kết quả đơn nhất:** vì sao ra đúng một kết quả, với xác suất $|c|^2$?

**Trong $\mathcal E_{\text{dual}}$, người quan sát là môi trường nhận năng lượng Carnot.** Đó là các nucleon hấp thụ xung ở mỗi va chạm. Va chạm truyền xung $\Lambda=m_{\text{eff}}u$ dọc pháp tuyến $\nabla f$, nên môi trường "đo" **vận tốc tương đối**

$$\hat u=(\kappa/I)\,\tfrac12\{\cos\theta,\hat p_\theta\}-\hat p_x .$$

Mình giải số câu 1 bằng einselection (mục A), câu 2–3 bằng quantum Darwinism và quy tắc Born (mục B). Mục C nói rõ phần nào còn lại.

## A. Cơ sở ưu tiên: người quan sát "chọn" thấy gì?

**Mô hình:**
- Hệ: rô-to ⊗ hạt, $H$ Bảng I với $-gx\sin\theta$.
- Người quan sát đo liên tục $\hat u$ với cường độ Γ: $\dot\rho=-i[H,\rho]+\Gamma(\hat u\rho\hat u-\tfrac12\{\hat u^2,\rho\})$.

**Predictability sieve (Zurek):** so độ tinh khiết $\operatorname{Tr}\rho^2$ của ba họ trạng thái đầu thuần. Họ nào giữ thuần lâu nhất là **trạng thái con trỏ**, tức cái người quan sát "thấy":
- **E:** 10 trạng thái riêng năng lượng thấp nhất. Đây là **sóng dừng**.
- **C:** 8 bó sóng "cổ điển". Rô-to định xứ quanh $\theta_0$ với $\ell_0=2,3$, hạt ở trạng thái kết hợp tại $x=\kappa\sin\theta_0$, $v=\kappa\ell_0\cos\theta_0$. Đây là **hình ảnh hạt** của bản thảo.
- **U:** trạng thái riêng của $\hat u$.

**Theo thời gian, tại điểm làm việc** $\Gamma=0.5$ (= tốc độ làm nguội Carnot $\Gamma_c$ của V.C.1), cơ sở $\ell\le7$, $n\le6$:

| $t$ | E (sóng dừng) | C (bó sóng cổ điển) | U (riêng $\hat u$) | Thắng |
|---|---|---|---|---|
| 0.02 $T_0$ | 0.871 | 0.924 | 0.984 | U |
| 0.05 $T_0$ | 0.750 | 0.834 | 0.863 | U |
| **0.10 $T_0$** | 0.624 | **0.713** | 0.663 | **C** |
| **0.25 $T_0$** | 0.393 | **0.442** | 0.211 | **C** |

- Phương sai $\hat u$ trung bình: E 1.29, C 0.68, U 0.
- So với cơ sở nhỏ ($\ell\le5$, $n\le4$), các họ E và C khớp tới ±0.003. Riêng họ U lệch mạnh theo cỡ cơ sở (0.33 → 0.21 ở 0.25 $T_0$), vì trạng thái riêng của $\hat u$ trải đến biên cơ sở. Họ U chỉ có ý nghĩa ở thời gian rất ngắn.

**Theo cường độ Γ, tại $t=T_0/2$** (cơ sở $\ell\le5$, $n\le4$):

| Γ | E | C | U | Thắng |
|---|---|---|---|---|
| 0.02 | 0.861 | **0.905** | 0.830 | C |
| 0.2 | 0.409 | **0.444** | 0.312 | C |
| 0.5 | **0.237** | 0.207 | 0.146 | E (sát nút) |
| 2 | **0.095** | 0.056 | 0.050 | E |
| 20 | 0.025 | **0.040** | 0.027 | C |

**Kết quả:**
1. **Bó sóng cổ điển (hình ảnh hạt) bền nhất trong khoảng 0.1–0.25 $T_0$** ở điểm làm việc của động cơ. Theo nghĩa einselection, người quan sát Carnot "thấy" quỹ đạo cổ điển: rô-to có góc, hạt bám theo $x=\kappa\sin\theta$. Đây là lý lẽ vì sao mô tả cổ điển của bản thảo có cơ sở.
2. **Nhưng ưu thế rất mỏng và không bền.** Mọi họ đều mất hơn nửa độ tinh khiết trong 0.25 $T_0$. Thứ hạng đổi theo Γ và theo thời gian. Lý do: tốc độ giải kết hợp $2\Gamma\,\mathrm{Var}(\hat u)\approx0.7$–$1.3\,\omega_0$ **ngang** tần số động lực học ($\omega_0$, khoảng cách mức rô-to). E_dual nằm trong **vùng trung gian**: không có cơ sở con trỏ sạch, không thuần sóng, không thuần hạt.
3. Điều này khớp với [rung-dong-m-song-hat-kiem-chung.md](rung-dong-m-song-hat-kiem-chung.md), nơi 68–95% thời gian là "đồng thời sóng–hạt". Hai cách đo khác nhau cho cùng kết luận.

## B. Tính khách quan và quy tắc Born: quantum Darwinism

**Mô hình:**
- Hệ: 3 mức quay của pha Quá khứ ($|\ell|=1,2,3$) ở chồng chập, trọng số Born $p=(0.2,\ 0.3,\ 0.5)$.
- Người quan sát: $N=10$ nucleon-qubit. Mỗi qubit $k$ ghi mức $j$ bằng một phép quay **không hoàn hảo** góc $j\varphi_k$.
- Toàn bộ là **đơn nhất**: trạng thái toàn cục $\lvert\Psi\rangle\rangle=\sum_j c_j\lvert j\rangle\bigotimes_k\lvert e_{jk}\rangle$ thuần, $S=10^{-15}$ bit. **Không có sụp đổ nào được đưa vào.**

### B.1. Giải kết hợp: ρ_S trở nên chéo khi môi trường ghi

Bản ghi mạnh ($\varphi_k\in[0.6,1.0]$ rad):

| Số qubit đã ghi | $\lvert\rho_S(1,2)\rvert$ | $\lvert\rho_S(2,3)\rvert$ | $\lvert\rho_S(1,3)\rvert$ |
|---|---|---|---|
| 0 | 0.245 | 0.387 | 0.316 |
| 2 | 0.093 | 0.147 | 0.014 |
| 6 | 0.019 | 0.030 | 0.000 |
| 10 | **0.004** | **0.006** | **0.000** |

### B.2. Nhiều người quan sát có cùng thông tin không?

Thông tin tương hỗ $I(S:F)$ theo kích thước mảnh $F$, với $H(S)=1.486$ bit:

| $\lvert F\rvert$ | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| Bản ghi mạnh: $I/H(S)$ | 0.58 | 0.89 | **0.94** | **0.98** | **0.99** | **1.02** | **1.05** | 1.15 | 1.42 | 2.00 |
| Bản ghi yếu ($\varphi\in[0.25,0.75]$): $I/H(S)$ | 0.40 | 0.70 | 0.76 | 0.88 | 0.93 | 1.08 | 1.17 | 1.34 | 1.57 | 1.96 |

- **Có bình nguyên (plateau)** ở $I\approx H(S)$ từ 3 đến 7 qubit. Chỉ cần đọc **một mảnh nhỏ** môi trường là biết gần đủ $H(S)$, và đọc thêm không cho thêm thông tin cổ điển về hệ. Đến $|F|=N$ thì $I=2H$, phần thông tin lượng tử còn lại.
- Độ dư thừa: bản ghi mạnh $R_{0.1}=3.3$, bản ghi yếu $R_{0.1}=2.0$. Khoảng 3 người quan sát độc lập, mỗi người đọc một mảnh riêng, **cùng biết** mức ℓ. Đây là **tính khách quan** sinh ra từ đơn nhất, không cần tiên đề sụp đổ.

### B.3. Quy tắc Born

Một người quan sát đọc $m$ qubit. Xác suất mỗi chuỗi bit tính **trực tiếp từ $|\Psi\rangle\rangle$**, không giả định trước nhánh nào. Người quan sát kết luận mức bằng hợp lý cực đại:

| $m$ | Tần suất kết luận (mức 1, 2, 3) | Born $\lvert c_j\rvert^2$ | Kết luận khớp nhánh |
|---|---|---|---|
| 1 | 0.34, 0.00, 0.66 | 0.2, 0.3, 0.5 | 69% |
| 3 | 0.22, 0.22, 0.57 | | 84% |
| 5 | 0.21, 0.31, 0.49 | | 93% |
| **10** | **0.200, 0.294, 0.506** | | **97%** |

Khi người quan sát đọc đủ nhiều bản ghi, tần suất kết luận **hội tụ về $|c_j|^2$**.

## C. Phần nào được giải, phần nào không

| Câu hỏi | Giải bằng số? | Cơ chế |
|---|---|---|
| Vì sao không thấy chồng chập? | ✓ | Môi trường Carnot ghi thông tin, làm $\rho_S$ chéo (B.1) |
| Thấy cái gì? | ✓ | Với kênh Davies thật của V.C.1: bó sóng cổ điển thỏa ràng buộc (pha Quá khứ) là trạng thái con trỏ (D). Với phép đo $\hat u$ Hermit, ưu thế mỏng (A) |
| Vì sao mọi người quan sát đồng ý? | ✓ | Bình nguyên Darwinism, $R_\delta\approx2$–3 (B.2) |
| Xác suất là $\lvert c\rvert^2$? | ✓ (thống kê) | Tần suất kết luận hội tụ Born (B.3) |
| **Vì sao ra đúng một kết quả cụ thể?** | **✗** | Trạng thái toàn cục vẫn chứa **cả ba** nhánh. Câu này phụ thuộc diễn giải (Copenhagen, Everett, QBism…). **Không phép tính số nào phân xử được** |

**Ráp vào học thuyết:**
- **Người quan sát của E_dual là chính bể Carnot.** Mũi tên thời gian ([pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md)) và tính khách quan (B) **cùng một nguồn**: năng lượng Carnot mang thông tin ra môi trường, không thể thu hồi.
- **Quỹ đạo cổ điển của bản thảo là một nhánh có điều kiện.** Mỗi lần mô phỏng Langevin là một "bản ghi" mà người quan sát Carnot đọc được. Dọc bản ghi đó, hệ thuần.
- **Ba pha Quá khứ → Hiện tại → Tương lai được định nghĩa bởi bản ghi.** Quá khứ là cái đã được ghi dư thừa, Hiện tại là lúc đang ghi (đợt va chạm), Tương lai là cái chưa ghi.

**Câu có thể đưa vào bản thảo (tiếng Anh):**

> The observer of the engine is the Carnot bath itself: each impact transfers momentum $\Lambda=m_{\text{eff}}u$ and thereby records the relative velocity. Within a unitary system–environment model, these records diagonalise the rotor state, are redundantly stored (redundancy $R_{0.1}\approx3$ for ten nucleon records), and reproduce Born frequencies. A predictability sieve at the operating dissipation $\Gamma_c=0.5\,\omega_0$ selects localized classical wave packets for $0.1$–$0.25\,T_0$, but only marginally, since the decoherence rate is comparable to $\omega_0$. The classical trajectories of this work are thus to be read as conditional records; the selection of a single outcome remains interpretation-dependent.

## D. Kiểm lại bằng kênh thật của V.C.1: làm nguội Davies $\sqrt{\Gamma_c}\,\hat u_-$

Mục A dùng phép đo $\hat u$ liên tục (Lindblad Hermit). Mục này dùng **đúng kênh làm nguội Carnot** của mô hình V.C.1. Kênh này vừa đo vừa lấy năng lượng đi. Lệnh chạy: `python3 giai_so_nguoi_quan_sat.py D` (≈ 15 phút).

**Độ tinh khiết trung bình theo thời gian.** C = bó sóng cổ điển (hình ảnh hạt), E = trạng thái riêng năng lượng (sóng dừng), U = trạng thái riêng $\hat u$:

| Kênh | $t=0.1\,T_0$ | $0.25\,T_0$ | $0.5\,T_0$ | $1\,T_0$ |
|---|---|---|---|---|
| $\Gamma_c=0.05$ | **C 0.998** / E 0.978 / U 0.955 | **C 0.991** / 0.948 / 0.893 | **C 0.972** / 0.904 / 0.803 | **C 0.941** / 0.835 / 0.681 |
| $\Gamma_c=0.5$ (V.C.1) | **C 0.982** / 0.836 / 0.711 | **C 0.933** / 0.730 / 0.533 | **C 0.807** / 0.683 / 0.410 | C 0.656 / **E 0.718** / 0.368 |
| $\Gamma_c=0.5$ + bơm $\Gamma_p=0.2$ (**đúng V.C.1**) | **C 0.942** / 0.710 / 0.641 | **C 0.813** / 0.499 / 0.418 | **C 0.560** / 0.343 / 0.255 | **C 0.358** / 0.241 / 0.167 |
| $\Gamma_c=5$ | **C 0.928** / 0.721 / 0.491 | **C 0.853** / 0.817 / 0.517 | C 0.726 / **E 0.822** / 0.505 | C 0.583 / **E 0.811** / 0.501 |

- Số liệu ở cơ sở $\ell\le5$, $n\le4$. Ở cơ sở $\ell\le7$, $n\le6$, các họ C và E khớp tới ±0.004. Chỉ họ U đổi theo cỡ cơ sở, như ở mục A.
- **Lưu ý khi đọc thời gian dài:** làm nguội Davies ở $T=0$ kéo mọi trạng thái về trạng thái cơ bản, vốn thuần. Họ E có độ tinh khiết **tăng trở lại** ($\Gamma_c=5$: 0.72 → 0.82), phản ánh sự hồi phục về cơ bản chứ không phải einselection. Vì vậy chỉ so sánh trong cửa sổ trước khi hồi phục.

**Vì sao bó sóng cổ điển thắng?** So sánh công bằng ba họ:

| Họ | $\langle H\rangle$ TB | Tốc độ nhảy Carnot $\langle L_c^\dagger L_c\rangle$ ($\Gamma_c=0.5$) | $\lvert\langle\hat u\rangle\rvert$ TB |
|---|---|---|---|
| E (sóng dừng) | 1.72 $\hbar\omega_0$ | 0.182 | 0 |
| **C (bó sóng cổ điển)** | **4.19** $\hbar\omega_0$ | **0.023** | 0.094 |
| U (riêng $\hat u$) | 7.06 | 0.376 | 0.151 |

Bó sóng cổ điển có năng lượng **cao gấp 2.4 lần** nhưng bị bể Carnot "đá" **ít hơn 8 lần**. Lý do là chúng **thỏa ràng buộc**: $x=\kappa\sin\theta$ và $v=\kappa\dot\theta\cos\theta$, nên vận tốc tương đối $u\approx0$ và va chạm Carnot không có gì để ghi. Trạng thái riêng năng lượng có $\langle\hat u\rangle=0$ nhưng phương sai lớn. Chúng là chồng chập của các chuyển động vi phạm ràng buộc theo cả hai chiều, nên bị ghi liên tục.

### Kết luận của mục D (thay cho kết luận 1–2 của mục A)

1. **Với kênh thật của V.C.1, người quan sát Carnot chọn rõ ràng bó sóng cổ điển thỏa ràng buộc**, tức chính **pha Quá khứ (khóa đồng bộ)** của bản thảo. Mô hình đầy đủ (làm nguội + bơm): C thắng **ở mọi mốc tới 1 $T_0$**, hơn E 0.10–0.31.
2. **Pha Quá khứ là trạng thái con trỏ của người quan sát Carnot.** Đây là liên hệ mới giữa cơ học giải tích của bản thảo và cơ học lượng tử: ràng buộc holonom $f=x-\kappa\sin\theta=0$ chính là điều kiện để bể **không** ghi được gì. Vì vậy trạng thái khóa là trạng thái "khách quan", bền nhất trước sự quan sát.
3. **Pha Hiện tại (gãy) là lúc ràng buộc bị vi phạm**, $u\neq0$, và bể ghi lại bằng các va chạm. Đây chính là lúc quan sát xảy ra. **Pha Tương lai** là khi hệ đã bị đẩy về một trạng thái con trỏ mới.
4. Ưu thế chỉ đảo khi làm nguội rất mạnh ($\Gamma_c=5$) **và** ở thời gian dài, khi mọi thứ đã hồi phục về trạng thái cơ bản. Đó là vùng ngoài tham số của bản thảo.

**Câu có thể đưa vào bản thảo (tiếng Anh):**

> With the Davies relaxation channel of Sec. V.C.1 as the observer, a predictability sieve selects localized wave packets that satisfy the holonomic constraint: their Carnot jump rate $\langle L_c^\dagger L_c\rangle$ is eight times smaller than that of energy eigenstates despite a 2.4-fold higher mean energy, because the relative velocity $u$ vanishes on the constraint surface. The synchronized Past phase is therefore the pointer state of the Carnot bath; breakdown into the Present is the act of observation, in which $u\neq0$ is recorded by inelastic impacts.

## Câu hỏi còn mở
- Không còn câu hỏi số nào cho phần này. Câu "vì sao ra đúng một kết quả" (mục C) vẫn là vấn đề diễn giải, không giải được bằng số.
