# Bốn câu hỏi phản biện cho Part I: kiểm chứng các phương án trả lời

> Số liệu mới tính bằng [giai_so_phan_bien.py](giai_so_phan_bien.py) (chạy `python3 giai_so_phan_bien.py`, ≈ 1 phút trên 8 lõi). Các số còn lại dẫn từ các tài liệu kiểm chứng trước.
> **Bộ thông số:** $\tau=10$, $\gamma_x=0.7$, $f=1.05$, $g=0.1$, $\hbar\omega_0=10$ MeV.
> **Chưa kiểm được bản thảo:** vẫn không có file `paper_prc.tex` trong dự án. Các mục V.B.3 và V.C.2 được viện dẫn nhưng chưa đối chiếu được.

## Tóm tắt

| Câu phản biện | Phương án trả lời | Kiểm chứng | Kết luận |
|---|---|---|---|
| **1.** Nguồn năng lượng của mô-men lái Van der Pol | Hệ mở; nguồn lái là biểu diễn hiện tượng học; "công lái **cân bằng chính xác với tiêu tán Carnot**" | Trong 60 $T_0$ (8 quỹ đạo, nhiễu 1 MeV): công lái **700.2**, Carnot **142.7** (chỉ **20.4%**), phần còn lại 553.9 đi vào ma sát $\gamma_x$ (trừ phần nhiễu đưa vào). Không nhiễu: 100% công lái đi vào ma sát, Carnot = 0 | **Sửa câu cuối:** công lái cân bằng với **ma sát $\gamma_x$ cộng tổn hao Carnot**, trong đó Carnot chỉ chiếm khoảng 20% |
| **2.** "Tự duy trì" hay "do nhiễu kích hoạt" | "Động cơ tự duy trì có kích hoạt thăng giáng"; lề 13.7%; $\Delta_{\text{eff}}\approx0.4$ MeV; "sau khi vượt ngưỡng, toàn bộ… **hoàn toàn theo động lực học nội tại**" | $\dot\theta_d=2.575$ ✓, $M_{\det}=0.863$ ✓, lề 13.7% ✓, $\Delta_{\text{eff}}\approx0.4$ MeV ✓ (ước lượng từ 2 điểm), 0 lần gãy khi không nhiễu ✓. **Nhưng** tổng Carnot mỗi đợt (1 MeV): **19.4** MeV khi không ngắt, **7.6** khi ngắt lái cả đợt, **5.9** khi ngắt cả lái lẫn nhiễu | **Phần định danh đúng. Sửa câu cuối:** nếu nguồn lái vẫn hoạt động trong đợt, nó cấp khoảng **60%** năng lượng của đợt, nên đợt **không** thuần nội tại. Chỉ thuần nội tại khi ngắt lái (và nhiễu) trong đợt |
| **3.** Đuôi năng lượng không bị chặn | "Ngắt nguồn lái **trong pha tự do** tự động chặn trên $E_{\text{budget}}\le20.4$ MeV" | Chỉ ngắt ở **pha tự do**: max **80.6** MeV (1 MeV) và **321.5** MeV (5 MeV). Ngắt lái **suốt đợt** (kể cả các pha khóa ngắn): max **19.9** MeV (1 MeV) nhưng **271.5** MeV (5 MeV). Ngắt **lái + nhiễu** suốt đợt: max **11.6** và **20.4** MeV. Con số 20.4 là **tổng Carnot lớn nhất**, không phải ngân sách (ngân sách max là 24–52 MeV) | **Sai điều kiện.** Phải ghi: ngắt lái **trong toàn bộ đợt**, và ở $k_BT\gtrsim$ vài MeV thì ngắt cả **nhiễu**. Khi đó tổng tiêu tán mỗi đợt ≤ 20 MeV |
| **4.** Cơ sở của $\gamma_x=0.70$ | Hấp thụ mặt $W_D\approx3.0$–4.5 MeV, cộng với $W_V$, cho $\gamma_{\text{total}}\approx0.74$ | Koning–Delaroche (neutron): $\gamma_{\text{total}}=2(W_V+W_D)/\hbar\omega_0$ đạt **0.71** tại $\Delta E=7$ MeV và **0.84** tại 8 MeV (A=56), nhờ $W_D=3.45$–4.07 MeV. **Nhưng** $\gamma$ phụ thuộc mạnh vào năng lượng: 0.43 (5 MeV) → 1.07 (10 MeV) | **Đứng được về số**, với 3 lưu ý cần ghi vào bài (bên dưới) |

## 1. Ngân sách công của nguồn lái

| Nhiễu | Công lái $W$ | Carnot $C$ | $C/W$ | $\Delta E$ hệ | $W-C-\Delta E$ | Số lần gãy |
|---|---|---|---|---|---|---|
| 0 (tất định, 1 quỹ đạo) | 97.14 | 0.00 | 0% | +0.70 | 96.43 | 0 |
| 1 MeV (8 quỹ đạo) | 700.16 | 142.65 | **20.4%** | +3.61 | 553.91 | 362 |

$W-C-\Delta E$ bằng năng lượng mất do ma sát $\gamma_x$ trừ đi năng lượng nhiễu đưa vào.

## 2. Hấp thụ Koning–Delaroche (neutron)

$$W_D=d_1\frac{\Delta E^2e^{-d_2\Delta E}}{\Delta E^2+d_3^2},\quad d_1=16\Big(1-\frac{N-Z}{A}\Big),\quad d_3=11.5\text{ MeV};\qquad \gamma_{\text{total}}=\frac{2(W_V+W_D)}{\hbar\omega_0}$$

| $\Delta E=E-E_F$ (MeV) | A=56: $W_D$ / $\gamma_{\text{total}}$ | A=70: $\gamma_{\text{total}}$ | A=208: $\gamma_{\text{total}}$ |
|---|---|---|---|
| 2 | 0.42 / 0.085 | 0.084 | 0.073 |
| 4 | 1.47 / 0.301 | 0.296 | 0.260 |
| 5 | 2.12 / 0.434 | 0.428 | 0.376 |
| 6 | 2.79 / 0.573 | 0.565 | 0.498 |
| **7** | **3.45 / 0.711** | **0.700** | 0.620 |
| 8 | 4.07 / 0.841 | 0.828 | 0.737 |
| 10 | 5.14 / 1.071 | 1.055 | 0.945 |

**Ba lưu ý cần ghi vào bài khi dùng lập luận này:**

1. **$W_D$ là độ sâu của thừa số dạng đạo hàm Woods–Saxon**, có đỉnh ở bề mặt ($r\approx R$). Nucleon trong mô hình dao động quanh tâm, nên độ hấp thụ trung bình trên hàm sóng **nhỏ hơn** độ sâu này.
2. **Với A = 56, ΔE = 7 MeV ứng với $E=-2.8$ MeV (trạng thái liên kết).** KD03 được khớp cho tán xạ ($E>0$). Dùng nó ở $E<0$ là ngoại suy. Mô hình quang học tán sắc (dispersive optical model) mới cho câu trả lời ở vùng này.
3. **Một hằng số $\gamma_x=0.70$ chỉ đúng tại ΔE ≈ 7 MeV.** Trong khoảng 5–10 MeV, $\gamma$ thay đổi từ 0.43 đến 1.07 (gấp 2.5 lần). Nên ghi rõ đây là giá trị hiệu dụng, và nêu độ nhạy của kết quả theo $\gamma_x$. Ở các kiểm chứng trước, $\gamma_x$ = 0.5 → 0.7 đã làm $r$ đổi từ 0.23 xuống 0.16/$T_0$.

## 3. Câu trả lời phản biện đề xuất (tiếng Anh, đã sửa)

> **R1 (energy source).** The Van der Pol torque is a phenomenological representation of angular-momentum transfer from collective degrees of freedom into the single-particle phase coordinate; the system is open and non-Hamiltonian. In steady state the drive work balances the combined losses from nucleon friction $\gamma_x$ and Carnot impacts. In our simulations (1 MeV noise, 60 $T_0$) Carnot impacts account for about 20% of the drive work, and friction for the remainder.

> **R2 (self-sustained vs. noise-activated).** The engine is noise-activated. Deterministically the attractor lies 13.7% below the breakdown stress ($M_{\det}=0.863$) and no breakdown occurs. Fluctuations with $k_BT\gtrsim1$ MeV overcome an effective barrier of $\approx0.4$ MeV, and the rate saturates above that. Once breakdown occurs, the cascade is governed by the internal dynamics only if the drive is decoupled during the burst; otherwise the drive supplies about 60% of the burst energy.

> **R3 (energy tails).** Bounded bursts require decoupling the drive for the whole burst, including the short re-locked intervals; decoupling only during free flight leaves tails up to 80 MeV at 1 MeV noise. With the drive decoupled for the whole burst, the maximum dissipated energy per burst is 19.9 MeV at 1 MeV noise. At several-MeV noise the bath must also be decoupled during the burst (maximum 20.4 MeV).

> **R4 (friction parameter).** Using Koning–Delaroche absorption, surface absorption dominates below 10 MeV, and $2(W_V+W_D)/\hbar\omega_0$ reaches 0.70 at $E-E_F\approx7$ MeV. We note three caveats: (i) $W_D$ is a surface-peaked depth, so the average over an interior orbit is smaller; (ii) this energy lies in the bound region, outside the KD fit domain; (iii) the absorption varies by a factor of 2.5 between 5 and 10 MeV. $\gamma_x=0.70$ should therefore be read as an effective value, and results are given for $\gamma_x=0.5$–0.7.

## Câu hỏi còn mở

- **Bản thảo.** Để mình đối chiếu các mục V.B.3 và V.C.2 và xác nhận các lỗi cũ đã được sửa, cần đặt `paper_prc.tex` vào thư mục `hoc-thuyet-e/`.
- **Trung bình $W_D$ trên hàm sóng.** Có thể tính $\langle W_D\rangle$ trên quỹ đạo dao động tử có $\ell_0=2.04$ fm với bán kính hạt nhân $R=1.2A^{1/3}$ fm. Việc này định lượng lưu ý (i), và có thể làm tiếp nếu cần.
