# Phương trình Euler–Lagrange hệ $\mathcal{E}_{\text{kép}}$ với biến thời gian $t$: đưa $t$ của Học thuyết E vào học tăng cường

> Cùng nội dung với [rl_bien_t_he_e_kep.tex](rl_bien_t_he_e_kep.tex) → [rl_bien_t_he_e_kep.pdf](rl_bien_t_he_e_kep.pdf).
> Mọi số liệu là kết quả chạy thật của [giai_so_rl_bien_t.py](giai_so_rl_bien_t.py) (phần B, D, S) và [giai_so_rl_thu_du_lieu.py](giai_so_rl_thu_du_lieu.py). Kết quả thô nằm trong `rl_data/` (`ket_qua_BD.txt`, `ket_qua_S.txt`, `chu_trinh_f.json`).
> Nguồn mô hình: HocThuyetE (8).pdf.

## Tóm tắt

Biến $t$ của Học thuyết E được đưa vào phương trình RL theo ba cách:

1. **Chiết khấu.** Quỹ đạo dừng của $\int e^{-\rho t}L\,dt$ thỏa phương trình Euler–Lagrange có thêm số hạng $-\rho\,\partial L/\partial\dot q$. Đây là dạng Caldirola–Kanai với ma sát âm, và năng lượng tăng theo $dE/dt=2\rho T$.
   - Mô-men Van der Pol là một chiết khấu phụ thuộc trạng thái: $\rho_{\text{eff}}(\dot\theta)=\gamma_\theta(1-\dot\theta^2/\Omega_0^2)/I$.
   - Thời gian lưu Poisson cho hệ số chiết khấu $\gamma=e^{-\Delta t/\tau_{\text{dw}}}$.
2. **Đồng hồ.** $t=k\tau$, với $\tau=T_0/25$ là nhịp của đồng hồ xung Page–Wootters. Giải phương trình Bellman horizon hữu hạn $V_k(s)$ trên đa tạp ràng buộc; biên $|\lambda|\ge R$ là trạng thái hấp thụ.
3. **Bán-Markov (SMDP).** Mỗi chu trình $\mathcal E_{\text{past}}\to\mathcal E_{\text{present}}\to\mathcal E_{\text{future}}$ là một bước, với thời gian lưu ngẫu nhiên $\tau_c$. Hệ số chiết khấu là $\mathbb E[e^{-\rho\tau_c}]$; nó có dạng đóng suy từ luật $T_0+\mathrm{Exp}(\mu)$ và $\mathrm{Geom}(p)$.

---

## 1. Phương trình Euler–Lagrange hệ $\mathcal{E}_{\text{kép}}$

### 1.1 Hàm Lagrange và ràng buộc

Không gian cấu hình là $(\theta,x)\in S^1\times\mathbb R$:

$$L=\tfrac12I\dot\theta^2+\tfrac12m\dot x^2+\Delta E\cos\theta-\tfrac12kx^2+gx\sin\theta$$

$L$ không chứa $t$ tường minh. Pha đồng bộ $\mathcal E_{\text{past}}$ chịu ràng buộc holonom

$$f(\theta,x)=x-\kappa\sin\theta=0,\qquad \kappa=\ell_0\,\Delta E/E_{vt}=\ell_0/3.$$

### 1.2 Lagrange loại một có lực suy rộng và điều khiển $u$

Các lực không thế giữ nguyên như trong bài báo:

- tiêu tán Rayleigh $\mathcal F=\tfrac12\gamma_x\dot x^2$;
- lực Langevin với $\langle\xi_i(t)\xi_i(t')\rangle=2\gamma_ik_BT\,\delta(t-t')$.

Riêng mô-men lái Van der Pol được thay bằng mô-men điều khiển $u(t)$ do tác tử RL chọn:

$$\frac{d}{dt}\frac{\partial L}{\partial\dot\theta}-\frac{\partial L}{\partial\theta}=u-\gamma_\theta\dot\theta+\xi_\theta-\lambda\kappa\cos\theta$$

$$\frac{d}{dt}\frac{\partial L}{\partial\dot x}-\frac{\partial L}{\partial x}=-\gamma_x\dot x+\xi_x+\lambda$$

Khử $\lambda$ bằng $\ddot x=\kappa\ddot\theta\cos\theta-\kappa\dot\theta^2\sin\theta$ thì được phương trình trên đa tạp một bậc tự do:

$$M(\theta)\ddot\theta=u-\gamma_\theta\dot\theta-\gamma_x\kappa^2\cos^2\theta\,\dot\theta+m\kappa^2\dot\theta^2\sin\theta\cos\theta-\Delta E\sin\theta-(k\kappa^2-2g\kappa)\sin\theta\cos\theta+\eta(t)$$

Trong đó $M(\theta)=I+m\kappa^2\cos^2\theta$ và $\eta=\xi_\theta+\kappa\cos\theta\,\xi_x$. Khi $u=\gamma=\eta=0$, đây đúng là phương trình (6) của bài báo.

Phản lực liên kết:

$$\lambda=m\ddot x+kx-g\sin\theta+\gamma_x\dot x$$

Pha $\mathcal E_{\text{past}}$ gãy khi $|\lambda|\ge R=E_{vt}/\ell_0$.

### 1.3 Biến $t$ của Học thuyết E

Vì $L$ không phụ thuộc $t$, thời gian chỉ đi vào hệ qua năm cấu trúc sau:

| # | Cấu trúc | Giá trị | Dùng ở |
|---|---|---|---|
| 1 | Chu kỳ $T_0=2\pi/\omega_0$ | thang đo | toàn bài |
| 2 | Thời gian lưu $T_0+\mathrm{Exp}(\mu)$ | $\tau_{\text{dw}}\approx4.5$–$6.0\,T_0$ | chiết khấu $\gamma$ |
| 3 | Số cú va mỗi đợt $N\sim\mathrm{Geom}(p)$ | $p\approx0.195$ | SMDP |
| 4 | Đồng hồ xung Page–Wootters | $\tau=T_0/25$, $t=k\tau$ | Bellman $V_k$ |
| 5 | Mũi tên thời gian | số bản ghi Carnot tăng đơn điệu | (diễn giải) |

---

## 2. Từ Euler–Lagrange sang phương trình Bellman

### 2.1 Định lý (EL chiết khấu): chiết khấu là ma sát âm

Quỹ đạo làm dừng $S_\rho=\int_0^Te^{-\rho t}L\,dt$ thỏa

$$\frac{d}{dt}\frac{\partial L}{\partial\dot q}-\frac{\partial L}{\partial q}=\rho\,\frac{\partial L}{\partial\dot q}.$$

Với hệ $\mathcal E_{\text{kép}}$:

$$I\ddot\theta=\rho I\dot\theta-\Delta E\sin\theta+gx\cos\theta,\qquad m\ddot x=\rho m\dot x-kx+g\sin\theta,\qquad \frac{dE}{dt}=2\rho T.$$

**Chứng minh.** Lấy biến phân: $\delta S_\rho=\int e^{-\rho t}\big[\partial_qL-\tfrac{d}{dt}\partial_{\dot q}L+\rho\,\partial_{\dot q}L\big]\delta q\,dt$. Nhân phương trình với $\dot q$ và dùng $\partial_{\dot q}L\cdot\dot q=2T$ thì được $dE/dt=2\rho T$.

Hệ quả:

- Đặt $\rho=-\gamma$ thì được Lagrange Caldirola–Kanai $e^{\gamma t}L$ của dao động tắt dần. Như vậy chiết khấu trong RL đối ngẫu với ma sát.
- Mô-men Van der Pol $\gamma_\theta(1-\dot\theta^2/\Omega_0^2)\dot\theta$ có đúng dạng $\rho_{\text{eff}}I\dot\theta$, với
  $$\rho_{\text{eff}}(\dot\theta)=\frac{\gamma_\theta}{I}\Big(1-\frac{\dot\theta^2}{\Omega_0^2}\Big).$$
  Chiết khấu này bằng $0.1\,\omega_0$ khi $\dot\theta\to0$ và triệt tiêu trên chu trình giới hạn $\dot\theta=\Omega_0$.
- Đối ngẫu chỉ chính xác cho phương trình của $\theta$. Vì $\gamma_\theta/I\ne\gamma_x/m$, không gộp được hai bậc tự do thành một thừa số $e^{\gamma t}$ chung.

### 2.2 Hamilton–Jacobi là Bellman với phần thưởng $L$

Đặt $S(q,t)=\operatorname{ext}\int_t^TL\,dt'$ và chọn điều khiển $v=\dot q$:

$$-\partial_tS=\operatorname{ext}_v\big[L+\partial_qS\cdot v\big]\ \Longrightarrow\ \partial_tS=H(q,-\partial_qS).$$

Đặc trưng của phương trình này chính là phương trình EL.

**Lưu ý:** phải dùng $\operatorname{ext}$, không dùng $\max$. $L$ lồi theo $\dot q$, nên tác dụng là điểm yên ngựa và "tối đa hóa $\int L$" không có nghiệm hữu hạn. Vì vậy trong RL thực, phần thưởng $r$ phải khác $L$; phương trình EL chỉ đóng vai mô hình chuyển trạng thái.

### 2.3 HJB có $t$ cho hệ $\mathcal{E}_{\text{kép}}$

Trạng thái là $s=(\theta,\dot\theta)$. Biên gãy là $\mathcal B=\{|\lambda(s,u)|\ge R\}$:

$$\rho V-\partial_tV=\max_u\Big\{r+\dot\theta\,\partial_\theta V+a(\theta,\dot\theta,u)\,\partial_{\dot\theta}V+\tfrac12\sigma^2(\theta)\,\partial^2_{\dot\theta}V\Big\},\quad V|_{\mathcal B}=-P,\ V(\cdot,T)=0$$

Ở đây $\sigma^2(\theta)=2k_BT(\gamma_\theta+\gamma_x\kappa^2\cos^2\theta)/M(\theta)^2$.

Rời rạc theo đồng hồ $t=k\tau$, với $j$ chạy trên các nút Gauss–Hermite của cú hích nhiễu:

$$V_k(s)=\max_u\sum_j w_j\Big[\mathbb 1_{\mathcal B_j}(-P)+(1-\mathbb 1_{\mathcal B_j})\big(r_j+\gamma V_{k+1}(s'_j)\big)\Big]$$

Hai lựa chọn cho $t$:

- **(a) Horizon hữu hạn:** $\rho=0$, và $V_k$ phụ thuộc $k$ tường minh.
- **(b) Dừng:** $\rho=1/\tau_{\text{dw}}$, tức coi pha $\mathcal E_{\text{past}}$ kết thúc theo hazard Poisson:
  $$\gamma=e^{-\tau/\tau_{\text{dw}}}=e^{-1/(25\times5.25)}=\mathbf{0.99241},\qquad \frac{1}{1-\gamma}=131.8\ \text{nhịp}=5.25\,T_0.$$

### 2.4 Bellman bán-Markov trên chu trình ba pha

Hành động là tỉ số lái $f$ ($\Omega_0=f\dot\theta_{\text{br}}$), chọn mỗi lần tái khóa:

$$R_c=-c_WP_W\frac{1-e^{-\rho\tau_c}}{\rho}+E_ce^{-\rho\tau_c},\qquad V=\max_f\Big\{\mathbb E[R_c|f]+\mathbb E[e^{-\rho\tau_c}|f]\,V\Big\}$$

Nếu $\tau_c=T_0+\mathrm{Exp}(\mu)+N\delta$, với $N\sim\mathrm{Geom}(p)$ độc lập, thì

$$\mathbb E[e^{-\rho\tau_c}]=\frac{e^{-\rho T_0}}{1+\rho\mu}\cdot\frac{p\,e^{-\rho\delta}}{1-(1-p)e^{-\rho\delta}}.$$

Q-learning bán-Markov không cần mô hình:

$$Q(a)\leftarrow Q(a)+\alpha\big[R_c+e^{-\rho\tau_c}\max_{a'}Q(a')-Q(a)\big]$$

Khi $\rho\to0$, mục tiêu trở thành tốc độ thưởng $\mathbb E[R_c]/\mathbb E[\tau_c]$ (định lý renewal–reward).

---

## 3. Giải số

Đơn vị: $\hbar\omega_0=1$ (10 MeV), $\omega_0=\ell_0=1$, $T_0=2\pi$.

Thông số: $I=m=k=\Delta E=1$, $\kappa=1/3$, $g=0.1$, $R=3$, $\gamma_\theta=0.1$, $\gamma_x=0.70$, $k_BT=0.1$.

### 3.1 Phần B: kiểm Định lý EL chiết khấu

Dùng tích phân biến phân rời rạc với cập nhật $q_{k+1}=q_k+e^{\rho h}(q_k-q_{k-1})-h^2M^{-1}\nabla\Pi(q_k)$, so với RK4 tại $t=4T_0$:

| $\rho$ | $h=T/1250$ | $h=T/2500$ | $h=T/5000$ | bậc đo |
|---|---|---|---|---|
| 0.00 | 1.47e-4 | 3.68e-5 | 9.21e-6 | 2.00 |
| 0.05 | 3.81e-3 | 1.98e-3 | 1.01e-3 | 0.97 |
| 0.10 | 1.48e-2 | 7.54e-3 | 3.81e-3 | 0.99 |

- Sai số hội tụ về 0, nên quỹ đạo dừng của tác dụng chiết khấu đúng là nghiệm của phương trình EL có thêm $\rho$.
- Khi $\rho>0$, bậc chỉ là 1 vì trọng số $e^{-\rho t_k}$ lấy ở điểm trái.
- Đẳng thức năng lượng: với $\rho=0.1$ trên $[0,2T_0]$, $\Delta E=0.161515$ và $\int2\rho T\,dt=0.161515$, lệch tương đối $6\times10^{-9}$.

### 3.2 Phần D: Bellman có đồng hồ trên đa tạp $\mathcal E_{\text{past}}$

**Bài toán**

- Hành động: $u\in\{-0.2,0,0.1,0.2,0.3,0.4,0.6\}$, chọn mỗi nhịp $T_0/25$.
- Phần thưởng mỗi nhịp: $\Delta\theta/2\pi-c\max(\int u\,d\theta,0)$, với $c=0.15$.
- Gãy: phạt $-P$ và kết thúc pha.
- Horizon $20\,T_0$ (500 nhịp); trạng thái đầu $(0,\,2.15)$.

**Lưới**

- $72\times169$ điểm, nội suy song tuyến, nhiễu Gauss–Hermite 5 nút.
- Trong mỗi nhịp, nhiễu được áp **trước**, sau đó mới kiểm $|\lambda|$.

**Kiểm ngưỡng gãy**

| Định nghĩa ngưỡng | Giá trị |
|---|---|
| Tiệm cận $\sqrt{R/m\kappa}$ | 3.000 |
| Theo quỹ đạo, mã Part I | 3.2048 (bài báo 3.205) |
| Tĩnh, không cản | 3.107 |
| Tĩnh, có cản | 3.056 |
| Tĩnh, có cản, $u=0.4$ | 3.038 |

Ngưỡng tĩnh là $\min\{\dot\theta:\max_\theta|\lambda|\ge R\}$. Hai định nghĩa cho số khác nhau vì $\dot\theta$ dao động theo $\theta$ trong một vòng.

**Kết quả**

| | $V$ ($P=2$) | P(gãy) | $V$ ($P=10$) | P(gãy) |
|---|---|---|---|---|
| Lưới, phụ thuộc $t$ | 33.574 | 0.013 | 33.492 | 0.008 |
| Lưới, dừng $\gamma=0.99241$ | 33.310 | 0.016 | 33.263 | 0.007 |
| SDE 4000 quỹ đạo, phụ thuộc $t$ | 33.892 ± 0.019 | 0.006 | 33.842 ± 0.022 | 0.003 |
| SDE 4000 quỹ đạo, dừng | 33.665 ± 0.030 | 0.007 | 33.616 ± 0.025 | 0.002 |

- SDE gốc (bước $\tau/10$, kiểm $|\lambda|$ ở mọi bước) xác nhận mô hình lưới trong khoảng 1%.
- Chính sách phụ thuộc $t$ hơn chính sách dừng **+0.23 vòng (≈0.7%, khoảng 6σ)**. Lợi thế nhỏ và tập trung ở cuối horizon:

| Thời gian còn lại | 20T0 | 10T0 | 5T0 | 2T0 | 1T0 | 0.4T0 | dừng |
|---|---|---|---|---|---|---|---|
| $\dot\theta$ lớn nhất còn đẩy tại $\theta=0$ ($P=2$) | 2.90 | 2.90 | 2.93 | 2.98 | 2.95 | 2.95 | 2.95 |
| $\langle u\rangle$ trên $\dot\theta\in(1.5,3)$ | 0.404 | – | 0.424 | – | 0.437 | **0.206** | – |

Gần hạn chót, tác tử làm hai việc mà chính sách dừng không làm được:

- chạy sát ngưỡng hơn, vì rủi ro gãy còn ít thời gian để xảy ra;
- trong 0.4 T0 cuối, giảm mô-men đi một nửa, vì công bỏ ra không kịp biến thành vòng quay.

> **Lỗi đã gặp và đã sửa:** bản đầu áp nhiễu *sau* bước tất định, nên bỏ sót các lần vượt ngưỡng giữa nhịp. Kết quả là lưới báo P(gãy) = 0.9%, trong khi SDE thật cho **98%**. Sau khi chuyển nhiễu lên trước, lưới và SDE khớp nhau.

### 3.3 Phần S: Bellman bán-Markov trên chu trình ba pha

**Dữ liệu.** Lấy từ mô phỏng Langevin đầy đủ của Part I (gãy, va Carnot dẻo $e=0$, tái khóa): 3 giá trị $f$ × 4 quỹ đạo × 120 $T_0$.

**Thông số đo được** ($c_W=0.15$, $S_n=0.8$):

| $f$ | Số chu trình | $\langle\tau_c\rangle/T_0$ | $\mu/T_0$ | CV dư | $\langle N\rangle$ | $\langle E_c\rangle$ | P($E_c\ge S_n$) | $P_W$ | Thưởng/$T_0$ |
|---|---|---|---|---|---|---|---|---|---|
| 0.95 | 16 | 19.27 | 17.83 | 0.49 | 2.12 | 0.900 | 0.500 | 0.1748 | −0.118 |
| 1.05 | 69 | 6.06 | 3.90 | 1.11 | 4.51 | 1.675 | 0.623 | 0.2303 | +0.060 |
| 1.13 | 75 | 5.59 | 2.40 | 1.20 | 8.61 | 3.620 | 0.800 | 0.2661 | +0.396 |

So với bài báo tại $f=1.05$:

- $\mu=3.90\,T_0$ và $\langle N\rangle=4.51$ thấp hơn $4.5$–$6.0\,T_0$ và $5.12$ của bài báo.
- CV dư 1.11 khớp với khoảng 1.0–1.1 của bài báo.
- Với 69 chu trình, sai số thống kê của $\mu$ khoảng ±12%.

**Dạng đóng so với thực nghiệm** ($f=1.05$, $\mu=3.90\,T_0$, $p=0.222$, $\delta=0.257\,T_0$):

| $\rho$ | $\mathbb E[e^{-\rho\tau_c}]$ thực nghiệm | Dạng đóng | Lệch |
|---|---|---|---|
| 0.01 | 0.7089 | 0.7028 | 0.9% |
| 0.03 | 0.4134 | 0.3901 | 6% |
| 0.10 | 0.1148 | 0.0864 | 25% |

Khi $\rho$ lớn, kết quả nhạy với các chu trình ngắn, và giả thiết $T_0+\mathrm{Exp}$ thuần túy đánh giá thấp chúng.

**Chính sách tối ưu theo sở thích thời gian $\rho$:**

| $\rho$ | $V(0.95)$ | $V(1.05)$ | $V(1.13)$ | $f^*$ |
|---|---|---|---|---|
| 0 (thưởng/$T_0$) | −0.118 | +0.060 | +0.396 | 1.13 |
| 0.003 | −6.523 | +2.408 | +18.133 | 1.13 |
| 0.01 | −2.101 | +0.259 | +3.814 | 1.13 |
| 0.03 | −0.771 | −0.243 | +0.366 | 1.13 |
| 0.1 | −0.251 | −0.233 | −0.211 | 1.13 |
| 0.3 | −0.087 | −0.109 | −0.125 | **0.95** |

Công lái phải trả ngay, còn năng lượng Carnot chỉ đến ở cuối chu trình. Vì thế khi $\rho$ đủ lớn, lựa chọn tối ưu đảo sang lái yếu nhất.

**Horizon hữu hạn theo $t$** ($H=60\,T_0$, $\rho=0$): $V(0)=20.72$.

- Chính sách là $f=1.13$ cho đến $t\approx57\,T_0$.
- Từ 57.1 đến 57.4 $T_0$, chính sách dao động giữa 1.05 và 1.13 (do nhiễu mẫu).
- Từ $58.1\,T_0$ chuyển sang $f=0.95$: chu trình không kịp khép trước hạn thì chỉ còn chi phí công lái.

**Q-learning bán-Markov** (lấy mẫu lại chu trình thực nghiệm, $2\times10^5$ cập nhật, $\varepsilon=0.2$, $\alpha_n=n^{-0.6}$):

| $\rho$ | $f^*$ học được | $Q^*$ học được | $V^*$ chính xác |
|---|---|---|---|
| 0.01 | 1.13 | 3.814 | 3.814 |
| 0.1 | 1.13 | −0.210 | −0.211 |

---

## 4. Gắn biến $t$ vào RL của bạn

1. **Chiết khấu từ thời gian lưu.** Nếu bước RL dài $\Delta t$, đặt $\gamma=e^{-\Delta t/\tau_{\text{dw}}}$ với $\tau_{\text{dw}}\approx5.25\,T_0$. Cách này thay một siêu tham số tùy ý bằng thang thời gian vật lý của pha $\mathcal E_{\text{past}}$.
2. **Đưa đồng hồ vào trạng thái.**
   - Bài toán có hạn chót: dùng trạng thái $(s,k)$, với $k$ là số đọc đồng hồ nhịp $T_0/25$.
   - Bài toán tuần hoàn: dùng $(s,\,t\bmod T_0)$.
3. **Bước bán-Markov.** Khi một hành động kéo dài một pha có độ dài ngẫu nhiên $\tau$, cập nhật theo
   $$Q(s,a)\leftarrow Q(s,a)+\alpha\Big[\int_0^\tau e^{-\rho t}r\,dt+e^{-\rho\tau}\max_{a'}Q(s',a')-Q(s,a)\Big].$$
   Nếu chỉ biết luật thời gian của học thuyết, dùng dạng đóng của $\mathbb E[e^{-\rho\tau}]$ ở mục 2.4.
4. **Biên gãy là trạng thái hấp thụ.** Kiểm $|\lambda(s,u)|\ge R$ trong suốt nhịp, kể cả ngay sau cú hích nhiễu, không chỉ ở cuối nhịp.
5. **Mô hình chuyển.** Dùng phương trình đa tạp (pha khóa) và EL tự do (pha tự do) làm mô phỏng môi trường. Không đặt $r=L$, vì khi đó bài toán tối đa hóa không có nghiệm hữu hạn.

Mẫu Python (SMDP Q-learning có $t$):

```python
import math

def smdp_q_update(Q, s, a, cycle, s_next, rho, alpha, c_w, P_w):
    """cycle = (tau, E_carnot): một chu trình E_past → E_present → E_future."""
    tau, E = cycle
    disc = math.exp(-rho * tau)                          # chiết khấu theo thời gian lưu thật
    R = -c_w * P_w * (1 - disc) / rho + E * disc         # chi phí liên tục + thưởng cuối đợt
    target = R + disc * max(Q[s_next].values())
    Q[s][a] += alpha * (target - Q[s][a])
```

---

## 5. Giới hạn

- **Phần thưởng và các giá ($c$, $c_W$, $P$) là lựa chọn thiết kế**, không suy ra từ học thuyết. Kết luận định tính (chính sách đổi theo $t$ và $\rho$) có thể giữ khi đổi phần thưởng, nhưng các con số cụ thể thì không.
- **Phần D là RL dựa trên mô hình** (quy hoạch động chính xác trên lưới), không phải RL học sâu. Lợi thế của biến $t$ ở đây chỉ 0.7%.
- **Phần S giả thiết trạng thái sau tái khóa là Markov.** Dữ liệu nhỏ: 16–75 chu trình mỗi giá trị $f$, vì máy đang chịu tải nặng (load trên 100). Bộ thu dữ liệu chạy tiếp được nếu bị ngắt; muốn mở rộng thì sửa `F_LIST` và `N_TRAJ` trong `giai_so_rl_thu_du_lieu.py`.
- **Dạng đóng của $\mathbb E[e^{-\rho\tau}]$ chỉ chính xác khi $\rho\mu\ll1$.**
- **Tích phân biến phân chỉ có bậc 1 khi $\rho>0$.** Dùng trọng số điểm giữa sẽ được bậc 2.

## Cách chạy lại

```bash
cd hoc-thuyet-e
python3 giai_so_rl_thu_du_lieu.py      # dữ liệu Langevin → rl_data/chu_trinh_f.json (lâu, tiếp tục được)
python3 giai_so_rl_bien_t.py B         # cầu EL ↔ Bellman (< 5 s)
python3 giai_so_rl_bien_t.py D         # Bellman có đồng hồ (vài phút khi máy rảnh)
python3 giai_so_rl_bien_t.py S         # SMDP + Q-learning (≈ 3 s)
tectonic rl_bien_t_he_e_kep.tex        # biên dịch PDF
```

## Câu hỏi còn mở

- Phần thưởng thật trong RL của bạn là gì: sản lượng vòng quay, năng lượng Carnot, hay xác suất phát nucleon? Câu trả lời quyết định $f^*$ và mức phạt $P$.
- Có cần thêm dữ liệu Langevin (5 giá trị $f$, 8 quỹ đạo, 150 $T_0$) để $\mu$ và $\langle N\rangle$ khớp với bài báo hơn không?
