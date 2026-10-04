# Phương trình Schrödinger trên đa tạp ràng buộc (pha Quá khứ) của $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** bản thảo [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Eq. (1)–(6) và Bảng I: $I=k=\Delta E=1$, $g=0.1$, $\kappa=1/3$, $\dot\theta_{\text{br}}=3.205$; đơn vị $\hbar=m=\omega_0=\ell_0=1$, $\hbar\omega_0=10$ MeV.
> Số liệu tính bằng [giai_so_schrodinger_rang_buoc.py](giai_so_schrodinger_rang_buoc.py). Lệnh `python3 giai_so_schrodinger_rang_buoc.py` chạy ≈ 70 s.
> **Bổ sung cho** [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md). Bài đó giải Schrödinger 2 chiều $(\theta,x)$ không ràng buộc. Bài này giải Schrödinger **trên chính đa tạp $x=\kappa\sin\theta$**, tức lượng tử hóa Eq. (6). Đây là nơi bản thảo đặt các khẳng định bán cổ điển (V.B.1–2).

## Mô hình

Thế $x=\kappa\sin\theta$ vào Eq. (1) cho Lagrangian một bậc tự do:

$$L=\tfrac12M(\theta)\dot\theta^2-U(\theta),\qquad M=I+m\kappa^2\cos^2\theta,\qquad U=-\Delta E\cos\theta+\left(\tfrac12k\kappa^2-g\kappa\right)\sin^2\theta.$$

Phương trình Euler–Lagrange của $L$ đúng là Eq. (6). Vì khối lượng phụ thuộc vị trí, cách sắp toán tử không duy nhất. Mình giải cả ba cách:
- **LB:** toán tử Laplace–Beltrami trên metric $M(\theta)$.
- **LB + da Costa:** thêm thế hình học $-\hbar^2K^2/8$. Đây là cách đúng khi ràng buộc được thực hiện bằng một giếng giam rất cứng. $K$ là độ cong của đường $x=\kappa\sin\theta$ trong tọa độ trọng-khối.
- **naive:** $-\tfrac12\partial_\theta M^{-1}\partial_\theta$.

Lưới tuần hoàn 512 điểm.

## S1. Phổ của pha Quá khứ

Năng lượng gãy (Eq. (10), $\dot\theta$ tại $\theta=0$): $E_{\text{br}}=\tfrac12M(0)\dot\theta_{\text{br}}^2+U(0)=\mathbf{4.706\,\hbar\omega_0}$ (47.1 MeV).

| Mức | LB | LB + da Costa | naive | $\sqrt{\langle p_\theta^2\rangle}$ | Rô-to tự do $\ell^2/2\bar M$ | Loại |
|---|---|---|---|---|---|---|
| 0 | −0.543 | −0.547 | −0.535 | 0.67 | | dao động quanh θ = 0 |
| 1 | 0.321 | 0.312 | 0.311 | 1.15 | | dao động |
| 2 | 0.842 | 0.835 | 0.841 | 1.01 | | dao động |
| 3, 4 | 1.968, 1.991 | 1.962, 1.986 | 1.967, 1.995 | 1.98, 1.94 | 1.897 (ℓ = 2) | quay, ±2 |
| **5, 6** | **4.306**, 4.306 | 4.300 | 4.306 | 2.99 | 4.269 (ℓ = 3) | **quay, ±3** |
| 7, 8 | 7.608 | 7.602 | 7.608 | 4.00 | 7.589 (ℓ = 4) | quay, ±4: **trên ngưỡng** |
| 9, 10 | 11.864 | 11.858 | 11.864 | 5.00 | 11.859 | quay, ±5 |

($\bar M=1/\langle1/M\rangle=1.0541$.)

**Kết quả:**
1. **Pha Quá khứ lượng tử chỉ có 7 trạng thái:** 3 mức dao động (libration) và hai cặp quay $\pm2$, $\pm3$. Đếm Bohr–Sommerfeld cho $N(E_{\text{br}})=6.3$, cùng bậc.
2. **Mức quay bền cao nhất là $|\ell|=3$. Mức $|\ell|=4$ đã nằm trên ngưỡng gãy.** Khẳng định "rotor angular momentum ∼ 3ħ" (Tóm tắt, V.B.1) được xác nhận theo nghĩa chính xác: $\ell=3$ là trạng thái quay cuối cùng còn khóa.
3. **Cách sắp toán tử gần như không quan trọng:** chênh ≤ 0.01 $\hbar\omega_0$ (0.1 MeV). Thế da Costa nhỏ vì độ cong của đường ràng buộc nhỏ ($\kappa=1/3$). Phổ trên không phụ thuộc lựa chọn lượng tử hóa.
4. Hai mức $\pm2$ tách nhau 0.023 do thế $\cos\theta$ trộn $+2$ với $-2$. Các cặp từ $\pm3$ trở lên gần như suy biến hoàn toàn (tách 10⁻⁴).

## S2. Bước lượng tử quanh ngưỡng và "va chạm dưới lượng tử"

| Bước | $\Delta E$ ($\hbar\omega_0$) | MeV |
|---|---|---|
| Carnot cổ điển $\langle\Delta E_{\text{Carnot}}\rangle$ (V.A, V.B.2) | 0.40–0.46 | 4.0–4.6 |
| **Lượng tử nhỏ nhất đi xuống từ $\ell=3$: $\ell=3\to2$** | **2.315** | 23.1 |
| $\ell=4\to3$ (từ mức gãy đầu tiên về mức khóa cao nhất) | 3.302 | 33.0 |
| MCWF của bản thảo (V.C.1), mỗi bước nhảy | ≈ 2.05 | 20.5 |
| Khoảng trống $E_{\text{br}}-E(\ell=3)$ | 0.400 | 4.0 |

**Kết quả:**
1. Trên đa tạp ràng buộc, không có mức nào cách mức $\ell=3$ chỉ 0.4 $\hbar\omega_0$. Bước nhỏ nhất là **2.3 $\hbar\omega_0$**, lớn gấp **5–6 lần** một va chạm Carnot cổ điển. Nghĩa là năng lượng mất đi phải đến theo **bậc thang của rô-to**, $(2\ell-1)/2\bar M$, chứ không phải theo $\hbar\omega_0$.
2. **Câu V.B.2 cần sửa.** Câu "resolved into discrete operator transitions in multiples of ħω0" không đúng cho pha khóa, vì ở đó hạt không có bậc tự do riêng nên không có lượng tử $\hbar\omega_0$ độc lập. Bước đúng là khoảng cách mức rô-to, ≈ 2.3 $\hbar\omega_0$ ở $\ell=3$. Con số này khớp cỡ với 2.05 $\hbar\omega_0$ mỗi bước nhảy trong mô phỏng MCWF của chính bản thảo (V.C.1). Như vậy hai mục V.B.2 và V.C.1 hiện đang mâu thuẫn nhau, và V.C.1 mới là cái đúng.
3. Mức $\ell=3$ nằm dưới ngưỡng đúng **0.400 $\hbar\omega_0$**, bằng một va chạm Carnot trung bình. Đây có thể chỉ là trùng hợp số học, vì $E_{\text{br}}$ do $R=3$ quyết định còn mức $\ell=3$ do $\bar M$ quyết định. Mình **không** đề xuất đưa nhận xét này vào bài.

## S3. Bó sóng rô-to ở điểm làm việc: pha còn xác định được bao lâu?

Bó sóng Gauss quanh $\theta=0$ với $p_\theta=M(0)\dot\theta_d=2.92\hbar$ ($\dot\theta_d=0.821\,\dot\theta_{\text{br}}$, V.B.4). Mình giải TDSE chính xác qua khai triển riêng, rồi so với **tập hợp cổ điển** (4000 quỹ đạo Eq. (6), cùng phân bố $(\theta,p_\theta)$ ban đầu).

Thước đo là độ kết hợp pha $\lvert\langle e^{i\theta}\rangle\rvert$. Bằng 1 nghĩa là pha hoàn toàn xác định, bằng 0 nghĩa là pha trải đều cả vòng.

| $\Delta p_\theta$ ($\hbar$) | $\Delta\theta$ (rad) | Mất pha, lượng tử ($T_0$) | Mất pha, cổ điển ($T_0$) | Lượng tử @ 4.5 $T_0$ | Cổ điển @ 4.5 $T_0$ | Hồi sinh lượng tử quanh 2 $T_0$ |
|---|---|---|---|---|---|---|
| 0.5 | 1.00 | 0.37 | 0.37 | 0.44 | 0.10 | 0.60 (1.93 $T_0$) |
| 1.0 | 0.50 | 0.24 | 0.24 | 0.19 | 0.02 | 0.52 (1.97 $T_0$) |
| 2.0 | 0.25 | 0.13 | 0.14 | 0.07 | 0.15 | 0.59 (2.10 $T_0$) |

("Mất pha" là lúc $\lvert\langle e^{i\theta}\rangle\rvert$ rơi xuống dưới 1/e giá trị đầu.)

**Kết quả:**
1. **Pha của rô-to nhòe hết trong 0.13–0.37 $T_0$**, ngắn hơn **12–46 lần** thời gian lưu 4.5–6 $T_0$ của pha Quá khứ.
2. **Sự nhòe này là cổ điển, không phải lượng tử.** Lượng tử và cổ điển cho cùng thời gian mất pha, sai khác ≤ 0.01 $T_0$. Lý do: $\Delta\ell\ge\tfrac12/\Delta\theta$ buộc tập hợp phải có độ trải vận tốc, và rô-to phi điều hòa tự làm nhòe pha. Mô hình cổ điển của bản thảo tránh được điều này chỉ vì nó theo dõi **một** quỹ đạo có $\theta$ và $\dot\theta$ xác định đồng thời. Lượng tử cấm điều đó với tác dụng ~3ħ.
3. **Khác biệt lượng tử duy nhất là hồi sinh.** Độ kết hợp quay lại 0.5–0.8 đều đặn mỗi ≈ 1.05 $T_0$, tức nửa chu kỳ hồi sinh $T_{\text{rev}}/2=2\pi\bar M/\hbar$, vì phổ $\ell^2/2\bar M$ gần như nguyên. Tập hợp cổ điển không có hồi sinh.
4. **Hệ quả cho bản thảo.** Ở mức ~3ħ, "góc pha θ của động cơ" chỉ có nghĩa thống kê trên khoảng thời gian dưới 0.4 $T_0$. Câu "the model operates as an effective classical analog" (V.B.1) nên ghi rõ: mỗi quỹ đạo cổ điển tương ứng với **một bản ghi đo liên tục**, không phải với trạng thái lượng tử của một rô-to cô lập. Kết luận này ăn khớp với kết quả Page–Wootters có môi trường ở [pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md): dọc theo một bản ghi của môi trường, hạt vẫn thuần.

## Tổng hợp cho bản thảo

| Khẳng định trong PDF (7) | Kết quả Schrödinger trên đa tạp | Đề xuất |
|---|---|---|
| Tóm tắt, V.B.1: rô-to ~3ħ | ✓ $\ell=3$ là mức quay khóa cao nhất, $\ell=4$ đã gãy | Có thể ghi "the highest locked rotor level is $\lvert\ell\rvert=3$; $\lvert\ell\rvert=4$ lies above $E_{\text{br}}$" |
| V.B.2: lượng tử hóa thành "multiples of ħω0" | ✗ Bước nhỏ nhất là 2.3 $\hbar\omega_0$ (thang rô-to) | Sửa thành "rotor-ladder transitions of ≈ 2.3 ħω0 near threshold, consistent with the 2.05 ħω0 MCWF jumps (V.C.1)" |
| V.B.1: mô hình là tương tự cổ điển hiệu dụng | ✓ có điều kiện. Pha nhòe sau ≤ 0.4 $T_0$ (cơ chế cổ điển); hồi sinh lượng tử mỗi ≈ 1.05 $T_0$ | Ghi rõ quỹ đạo cổ điển là một bản ghi đo, không phải trạng thái rô-to cô lập |
| Lượng tử hóa không phụ thuộc cách chọn | ✓ Ba cách sắp toán tử chênh ≤ 0.1 MeV | Không cần nhắc |

## Câu hỏi còn mở
- Muốn tính bước nhảy lượng tử có **cả** hạt (rô-to ⊗ dao động tử) với ràng buộc mềm, có thể dùng lại mô hình MCWF V.C.1. Phân bố bước nhảy của nó (2.05 $\hbar\omega_0$) đã khớp thang rô-to ở đây. Có cần chạy lại để tách phần $\ell\to\ell-1$ khỏi phần $n\to n-1$ không?
