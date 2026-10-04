# Schrödinger, Page–Wootters, Wheeler–DeWitt cho học thuyết $\mathcal E_{\text{dual}}$ (PDF 7): giải số tổng hợp

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I: $I=k=\Delta E=1$, $g=0.1$, $\kappa=1/3$, $R=3$, $\hbar\omega_0=10$ MeV. Mô hình lượng tử theo V.C.1, đồng hồ theo V.C.4.
> **Phần mới:** [giai_so_wdw_pw_e7.py](giai_so_wdw_pw_e7.py), chạy ≈ 2 s. Dựng **một** trạng thái Wheeler–DeWitt trọn vẹn đúng cấu hình PDF (7): đồng hồ chuỗi xung ⊗ $\mathcal E_{\text{dual}}$ ⊗ môi trường Carnot.
> **Phần tổng hợp:** gom kết quả đã giải ở các tài liệu trước. Tham số của chúng trùng Bảng I của PDF (7), nên dùng được nguyên vẹn.

## 1. Phần mới: một vũ trụ Wheeler–DeWitt trọn vẹn cho PDF (7)

### Cấu hình

| Thành phần | Mô hình | Nguồn trong PDF (7) |
|---|---|---|
| **Đồng hồ C** | Bộ đếm xung, mỗi tích tắc $\tau=T_0/25$ (25 xung mỗi chu kỳ) | V.C.4 |
| **Hệ S** | Rô-to $\ell\le4$ ⊗ hạt $n\le3$ (36 chiều), $H$ Bảng I với $-gx\sin\theta$ | Eq. (1), Bảng I |
| **Môi trường E** | Mỗi tích tắc một ancilla 3 mức {không, nhận lượng tử Carnot $L_c=\sqrt{0.5}\,\hat u_-$, nhận lượng tử bơm $L_p=\sqrt{0.2}\,e^{i\theta}f(\hat\ell)$} | V.C.1 |
| Va chạm S–E | Đơn nhất $U=\exp\{-i[\tau H\otimes1+\sqrt\tau\sum_j(L_j\otimes\lvert j\rangle\langle0\rvert+\text{h.c.})]\}$, $\lVert U^\dagger U-1\rVert=3\cdot10^{-14}$ | |
| Trạng thái đầu của S | Bó sóng thỏa ràng buộc (pha Quá khứ), $\ell_0=3$, $\theta_0=0$, hạt kết hợp tại $x=0$, $v=\kappa\ell_0$ | Trạng thái con trỏ ([nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md), mục D) |

**Trạng thái lịch sử**, gồm 7 tích tắc, $E$ có $3^7=2187$ chiều, tổng $8\times36\times2187$ chiều:

$$\lvert\Psi\rangle\rangle=\tfrac{1}{\sqrt8}\sum_{k=0}^{7}\lvert k\rangle_C\otimes U_k\cdots U_1\lvert\psi_0\rangle_S\lvert0\cdots0\rangle_E .$$

### (1) Wheeler–DeWitt: vũ trụ tĩnh

| Kiểm tra | Kết quả |
|---|---|
| Ràng buộc rời rạc $\hat J\lvert\Psi\rangle\rangle=\lvert\Psi\rangle\rangle$ (dạng WdW cho đồng hồ rời rạc, Giovannetti–Lloyd–Maccone) | $\lVert\hat J\Psi-\Psi\rVert=$ **4.0·10⁻¹⁵** |
| Phân bố số đọc đồng hồ $P(k)$ | **1/8 cho mọi k**. Người quan sát ngoài không thấy gì thay đổi |
| Entropy toàn vũ trụ | **0**, trạng thái thuần |

### (2) Page–Wootters: thời gian xuất hiện cho người quan sát bên trong

Trạng thái điều kiện $\rho_S(k)=\operatorname{Tr}_E\langle k\vert\Psi\rangle\rangle$:

| $k$ | $t/T_0$ | Lệch so với Lindblad liên tục | $\operatorname{Tr}\rho_S^2$ | $S(\rho_S)=S(E)$ | $P$(≥ 1 bản ghi Carnot) | $P$(≥ 1 bản ghi bơm) |
|---|---|---|---|---|---|---|
| 0 | 0 | 0 | 1.000 | 0.000 | 0 | 0 |
| 1 | 0.04 | 0.003 | 0.978 | 0.063 | 0.011 | 0.020 |
| 3 | 0.12 | 0.009 | 0.937 | 0.155 | 0.064 | 0.065 |
| 5 | 0.20 | 0.017 | 0.882 | 0.271 | 0.155 | 0.116 |
| 7 | 0.28 | 0.023 | **0.797** | **0.440** | **0.251** | 0.163 |

- Khi đọc đồng hồ ở xung $k$, hệ S **tiến hóa đúng phương trình chủ** (Lindblad của V.C.1). Sai lệch 0.3% mỗi xung là sai số $O(\tau)$ của mô hình va chạm, không phải lỗi của PW.
- **Mũi tên thời gian nằm trong vũ trụ tĩnh.** $S(E)$ và số bản ghi Carnot **tăng đơn điệu** theo số đọc đồng hồ, trong khi entropy toàn cục vẫn bằng 0. "Tương lai" là phía đồng hồ có nhiều bản ghi hơn.

### (3) Điều kiện hóa trên bản ghi môi trường, tại $k=7$

| Bản ghi (· không, C Carnot, B bơm) | Xác suất | $\operatorname{Tr}\rho_S^2$ | $\langle H\rangle$ |
|---|---|---|---|
| ······· (chưa ghi gì) | 0.627 | **1.000000** | 4.608 $\hbar\omega_0$ |
| ·····C· | 0.039 | **1.000000** | 4.343 |
| ······C | 0.038 | **1.000000** | 4.053 |

- Đọc thêm **bản ghi của môi trường** thì S **thuần** ở mọi nhánh. Mỗi nhánh là một quỹ đạo lượng tử, đúng hình ảnh "quỹ đạo cổ điển = bản ghi" của bản thảo.
- Nhánh có bản ghi Carnot có năng lượng **thấp hơn** (4.61 → 4.05–4.34), vì lượng tử Carnot đã mang năng lượng đi.

**Kết luận của phần mới.** Ba phương trình **không mâu thuẫn** nhau trong học thuyết, mà là ba tầng của cùng một đối tượng:
- **Wheeler–DeWitt** mô tả cả vũ trụ (tĩnh, thuần).
- **Page–Wootters** cho thời gian khi đọc đồng hồ chuỗi xung.
- **Schrödinger / phương trình chủ** là luật tiến hóa mà người bên trong thấy.
- Ba pha Quá khứ → Hiện tại → Tương lai được **ghi lại trong môi trường Carnot**, nên có hướng dù vũ trụ tĩnh.

## 2. Tổng hợp các kết quả đã giải

### 2.1. Phương trình Schrödinger

| Kết quả | Số liệu | Tài liệu |
|---|---|---|
| Phổ 2D $(\theta,x)$ không ràng buộc | Kích thích thấp nhất 0.866 và 1.011 $\hbar\omega_0$, so với mode cổ điển 0.949 và 1.049. Phi điều hòa 9–12% | [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md) |
| Liên kết $x=\kappa\sin\theta$ lượng tử | Giếng có lực ≤ R giữ hạt với độ nhòe **≥ 1.57κ**. Ràng buộc cứng không thực hiện được | như trên |
| "Gãy" lượng tử (TDSE 2D) | Xác suất thoát tăng **liên tục**: 1.4% ở ℓ = 3, 14% ở ℓ = 6. Không có ngưỡng sắc | như trên |
| Schrödinger **trên đa tạp ràng buộc** (Eq. 6) | Pha Quá khứ có 7 mức. **$\lvert\ell\rvert=3$ là mức khóa cao nhất**, $\lvert\ell\rvert=4$ đã gãy. Cách sắp toán tử lệch ≤ 0.1 MeV | [schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md) |
| Bước lượng tử quanh ngưỡng | ℓ 3→2 = **2.3 $\hbar\omega_0$**, gấp 5–6 lần Carnot cổ điển 0.4. V.B.2 ("multiples of ħω0") cần sửa | như trên |
| Bó sóng rô-to ở $\dot\theta_d$ | Pha nhòe sau 0.13–0.37 $T_0$ (cơ chế cổ điển). Hồi sinh lượng tử mỗi ≈ 1.05 $T_0$ | như trên |
| Bước nhảy MCWF (V.C.1) | **92% là rô-to ℓ → ℓ−1**, 99.7% năng lượng. Hạt $n\to n-1$ chỉ 4% | [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md) |
| Dirac (tương đối tính) | Lệch −0.2% đến −1.4%. Schrödinger là đủ | [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md) |
| Sóng–hạt theo năng lượng rung | $V=e^{-2E/\hbar\omega_0}$, $D=\sqrt{1-V^2}$. Dọc chu trình, 68–95% thời gian là "đồng thời" | [rung-dong-m-song-hat-kiem-chung.md](rung-dong-m-song-hat-kiem-chung.md) |

### 2.2. Page–Wootters

| Kết quả | Số liệu | Tài liệu |
|---|---|---|
| Nghiệm WdW "gần 0" của script đề xuất | Thực ra là trạng thái cơ bản: đồng hồ **không quay**. Đồng hồ đúng nghĩa cần cặp ±ℓ gần suy biến ($\langle\ell\rangle=2.98$) | [page-wootters-kiem-chung.md](page-wootters-kiem-chung.md) |
| Rô-to làm đồng hồ | Tuần hoàn, $\lvert\langle\chi(\theta)\vert\chi(\theta+2\pi)\rangle\rvert^2=1$. Không đếm vòng | [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md) |
| Đồng hồ 30ħ và ~3ħ (cỡ rô-to của mô hình) | Độ trùng với Schrödinger 0.986–0.999 (30ħ); **0.946–0.975** (~3ħ) | như trên |
| Đồng hồ photon (Moreva) | Thí nghiệm lý tưởng tái lập với độ trùng 1.000000. Động cơ xếp đúng Quá khứ → Hiện tại → Tương lai ở 100% số đợt khi $k_BT\le1$ MeV | [dong-ho-photon-kiem-chung.md](dong-ho-photon-kiem-chung.md) |
| OAM và chuỗi xung | OAM tuần hoàn 2π, không đếm vòng. Chuỗi xung đếm được, cần ≥ 25 xung mỗi vòng (đã đưa vào V.C.4) | [pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md) |
| PW có môi trường | $\operatorname{Tr}_E$ = Lindblad. Độ tinh khiết giảm **trơn**. Dọc một bản ghi, mỗi lần phát Carnot làm hạt **thuần hơn** | như trên |
| Người quan sát Carnot | Trạng thái con trỏ là **bó sóng thỏa ràng buộc** (pha Quá khứ), bị "đá" ít hơn 8 lần. Darwinism $R_\delta\approx3$; Born hội tụ | [nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md) |

### 2.3. Wheeler–DeWitt

| Kết quả | Số liệu | Tài liệu |
|---|---|---|
| Nghiệm riêng không suy biến của $\hat H$ | $\langle\ell\rangle=0$ (đối xứng đảo thời gian). Đồng hồ quay hai chiều như nhau | [page-wootters-kiem-chung.md](page-wootters-kiem-chung.md) |
| Giới hạn WdW với đồng hồ ~3ħ | Độ trùng giảm còn **0.69–0.91**. Đối xứng $t\leftrightarrow\tfrac12T_{\text{rot}}-t$, **không có mũi tên** | [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md) |
| Thời gian trong AdS | Hạt trong AdS lặp lại sau $2\pi L/c$. Δ không nguyên buộc thời gian ∈ ℝ, trùng đòi hỏi Φ ∈ ℝ của V.C.4 | [ads-cft-toan-anh-kiem-chung.md](ads-cft-toan-anh-kiem-chung.md) |
| **WdW trọn vẹn có môi trường (mới)** | Ràng buộc thỏa tới 4·10⁻¹⁵. Mũi tên thời gian **có**, nằm trong bản ghi Carnot | Mục 1 ở trên |

## 3. Những gì cần ghi vào bản thảo

1. **V.B.2:** thay "multiples of ħω0" bằng "rotor-ladder transitions ≈ 2.3–3.5 ħω0". Câu tiếng Anh có sẵn trong [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md).
2. **V.C.4:** bổ sung rằng trạng thái WdW trọn vẹn (đồng hồ chuỗi xung ⊗ động cơ ⊗ bể Carnot) thỏa $\hat J\lvert\Psi\rangle\rangle=\lvert\Psi\rangle\rangle$, và mũi tên thời gian chính là bản ghi Carnot tăng dần. Đề xuất câu:

> A complete timeless state, $\lvert\Psi\rangle\rangle=N^{-1/2}\sum_k\lvert k\rangle_C\otimes U_k\cdots U_1\lvert\psi_0\rangle_S\lvert0\rangle_E$, with a pulse-train clock ($\tau=T_0/25$), the engine and a Carnot record register, satisfies the discrete Wheeler–DeWitt constraint to $4\times10^{-15}$; conditioning on the clock reproduces the Lindblad dynamics of Sec. V.C.1, while the environmental entropy and the number of Carnot records grow monotonically with the clock reading, so the arrow of time is carried by the records within a globally static universe.

3. **V.B.1:** "effective classical analog" nên hiểu là **bản ghi có điều kiện** ([nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md)).

## Giới hạn của phần mới
- Chỉ 7 tích tắc (0.28 $T_0$), vì môi trường tăng theo $3^k$. Thời gian dài hơn cần chuyển sang quỹ đạo lượng tử (MCWF). Hai cách là tương đương, và các kết quả dài hạn đã có ở các tài liệu MCWF.
- Mô hình va chạm lệch Lindblad 0.3% mỗi xung. Đã kiểm đây đúng là sai số $O(\tau)$: trên cùng khoảng thời gian τ, chia thành 1, 2, 4 xung cho độ lệch 3.5·10⁻³, 1.7·10⁻³, 8.8·10⁻⁴, tức giảm đúng một nửa mỗi lần chia đôi. Muốn giảm thì cần nhiều xung hơn mỗi chu kỳ, tức càng nhiều chiều.
