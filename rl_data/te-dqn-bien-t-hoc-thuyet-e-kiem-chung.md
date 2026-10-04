# TE-DQN: học tăng cường và Deep Q-Network có biến thời gian $t$ của Học thuyết E

> Tiếp nối [rl-bien-t-he-e-kep-kiem-chung.md](rl-bien-t-he-e-kep-kiem-chung.md). Bản đó chỉ có quy hoạch động và Q-learning dạng bảng; bản này thêm học sâu: Deep Q-Network (DQN) có biến $t$.
>
> **Nguồn đã đọc**
> - Cả 9 bản HocThuyetE (`HocThuyetE.pdf`, `(1)`…`(8)`). Đây là các bản sửa liên tiếp của cùng một bài báo, nên các số về thời gian lấy theo bản (8).
> - Thư mục `hoc-thuyet-e/`: mô phỏng Langevin, luật đợt Geom, tái khóa, đồng hồ Page–Wootters.
> - Sách *Cơ học giải tích* (Nguyễn Văn Đạo, NXB ĐHQG 2001, 348 trang). Bản PDF là ảnh quét, nên tôi tự OCR tiếng Việt bằng `tesseract -l vie`, đọc mục lục tra cứu rồi các trang 25–28, 72–75, 133–139, 189, 201–203, 213–215, 255–264.
>
> **Mã và kết quả**
> - Mã: [giai_so_te_dqn.py](giai_so_te_dqn.py) và [giai_so_e_dqn.py](giai_so_e_dqn.py). DQN viết tay bằng numpy vì máy không có torch.
> - Kết quả thô: `te_dqn_*.json`, `e_dqn_*.json`, `te_dqn_dp.json`, `*_train*.log` (cùng thư mục `rl_data/`).

---

## 1. Biến $t$ lấy từ đâu

### 1.1 Từ sách Cơ học giải tích (Nguyễn Văn Đạo)

| Chỗ trong sách | Nội dung | Vai trò của $t$ trong RL |
|---|---|---|
| tr. 25–28, §1.4: $\vec r_k=\vec r_k(q_1,\dots,q_n,t)$, liên kết **rêônôm** (không dừng) | $t$ xuất hiện tường minh trong tọa độ | MDP **không dừng**: $P(s'\mid s,a,t)$, $V(s,t)$ |
| (3.43)–(3.44): $\dfrac{dH}{dt}=\sum_i Q_i^*\dot q_i$ với hệ sclêrônôm | công suất của lực không thế là tốc độ đổi năng lượng | phần thưởng công $c\int u\,d\theta$; mô-men $u$ chính là $Q^*_\theta$ |
| (3.47): $\sum Q_i^*\dot q_i\le0$, hệ hao tán | ma sát làm năng lượng giảm | chiết khấu $\rho$ đối ngẫu với ma sát (Định lý EL chiết khấu, bản trước) |
| (5.26)–(5.27): $S=\int_{t_1}^{t_2}L(q,\dot q,t)\,dt$, $\delta S=0$ | nguyên lý Hamilton | tác dụng ↔ lợi tức; phải dùng ext thay cho max |
| (7.4): $\dfrac{\partial S}{\partial t}+H\big(q,\dfrac{\partial S}{\partial q},t\big)=0$ | phương trình Hamilton–Jacobi chuẩn | đây là HJB khi $\rho=0$; $V(s,t)$ là hàm sinh |
| (7.16): hệ bảo toàn $S=-ht+W(q)$ | tách biến thời gian | thế định hình $\Phi(s,t)\propto\dot\theta\,(H-t)$ (biến thể H) |

### 1.2 Từ Học thuyết E (bản 8 và thư mục)

| Cấu trúc thời gian | Giá trị | Trong TE-DQN |
|---|---|---|
| Nhịp đồng hồ xung Page–Wootters | $\tau=T_0/25$ | bước quyết định; $k$ là số đọc đồng hồ |
| Thời gian lưu pha $\mathcal E_{\text{past}}$: $T_0+\mathrm{Exp}(\mu)$ | $\tau_{\text{dw}}\approx5.25\,T_0$ | $\rho=1/\tau_{\text{dw}}$, nên $\gamma=e^{-\tau/\tau_{\text{dw}}}=0.99241$ |
| Đợt $\mathcal E_{\text{present}}$: $N\sim\mathrm{Geom}(p)$ cú va, khoảng cách $\delta$ | $p=0.195$, $\delta=0.36\,T_0$ = 9 nhịp | bước bán-Markov: đồng hồ nhảy $\Delta k=1+9N$ |
| Tái khóa $\mathcal E_{\text{future}}$ | $\dot\theta_{\text{end}}\approx2.15$ | trạng thái sau đợt |
| Ngưỡng gãy $|\lambda|\ge R$ | $R=3$ | sự kiện chuyển pha, kiểm ở từng bước con |

---

## 2. Môi trường: hệ $\mathcal E_{\text{kép}}$ đủ ba pha

**Trạng thái vật lý:** $(\theta,\dot\theta)$ trên đa tạp $x=\kappa\sin\theta$, cộng số đọc đồng hồ $k$.

**Pha $\mathcal E_{\text{past}}$.** Phương trình rút gọn (6) của bài báo, có thêm cản $\gamma_\theta$, cản $\gamma_x$ chiếu lên $\theta$, nhiễu Langevin $k_BT=0.1$ và mô-men điều khiển $u$:

$$M(\theta)\ddot\theta=u-\gamma_\theta\dot\theta-\gamma_x\kappa^2\cos^2\theta\,\dot\theta+m\kappa^2\dot\theta^2\sin\theta\cos\theta-\Delta E\sin\theta-(k\kappa^2-2g\kappa)\sin\theta\cos\theta+\eta(t)$$

Mỗi nhịp đồng hồ gồm 10 bước con Euler–Maruyama, kiểm $|\lambda|\ge R$ ở mọi bước con.

**Pha $\mathcal E_{\text{present}}$.** Khi gãy, đợt Carnot kéo dài $N\sim\mathrm{Geom}(0.195)$ cú va, mỗi cú 9 nhịp. Trong đợt động cơ mất đồng bộ nên không tính vòng quay.

**Pha $\mathcal E_{\text{future}}$.** Tái khóa với $\dot\theta=2.15$ và $\theta\sim U(-\pi,\pi)$. Đồng hồ đã tiến $1+9N$ nhịp.

**Hành động và phần thưởng.**
- $u\in\{-0.2,0,0.1,0.2,0.3,0.4,0.6\}$.
- $r=\Delta\theta/2\pi-0.15\max(u\,\Delta\theta,0)$: số vòng quay trừ giá công lái, theo (3.44).
- Horizon $H=20\,T_0=500$ nhịp. Chi phí của việc gãy tự sinh ra từ thời gian mất trong đợt và tốc độ bị kéo về 2.15, không cần đặt phạt $P$ tùy ý như bản trước.

Thông số lấy từ Bảng I của bài báo: $I=m=k=\Delta E=1$, $\kappa=1/3$, $g=0.1$, $R=3$, $\gamma_\theta=0.1$, $\gamma_x=0.7$, $k_BT=0.1$.

---

## 3. Thuật toán mới: TE-DQN (Temporal-E Deep Q-Network)

### 3.1 Phương trình Bellman có $t$

Trạng thái mở rộng là $\tilde s=(s,k)$. Với chuyển bán-Markov kéo dài $\Delta k$ nhịp:

$$Q^*(s,k,a)=\mathbb E\Big[r+\Gamma(\Delta k)\,\max_{a'}Q^*(s',k+\Delta k,a')\Big],\qquad Q^*(\cdot,k\ge K,\cdot)=0$$

Hệ số chiết khấu $\Gamma(\Delta k)$ lấy theo một trong ba cách:

| Cách | $\Gamma(\Delta k)$ | Nguồn |
|---|---|---|
| horizon hữu hạn đúng | $1$ | (7.4) với $\rho=0$; $t$ có trong trạng thái |
| chiết khấu Poisson bán-Markov | $e^{-\rho\,\tau\,\Delta k}$, $\rho=1/\tau_{\text{dw}}$ | thời gian lưu $T_0+\mathrm{Exp}(\mu)$ |
| DQN thường | $\gamma$ cố định | không dùng $t$ |

### 3.2 Năm thành phần của TE-DQN

1. **Đặc trưng thời gian.** Đầu vào mạng là $\big(\cos\theta,\sin\theta,\dot\theta/3,\lambda/R,(K-k)/K\big)$; thời gian còn lại lấy từ đồng hồ Page–Wootters. $\lambda/R$ là khoảng cách tới ngưỡng gãy, tính bằng phương trình Lagrange loại một.
2. **Đích bán-Markov.** $y=r+\Gamma(\Delta k)\,\bar Q\big(s',k',\arg\max_{a'}Q(s',k',a')\big)$ (Double DQN). $\Delta k$ là thời gian trôi thật, gồm cả đợt Carnot.
3. **Kết thúc đúng.** Chỉ khi $t$ có trong trạng thái thì hết horizon mới là trạng thái kết thúc thật. Nếu không có $t$, hết horizon chỉ là **cắt** và vẫn phải bootstrap (Pardo et al. 2018).
4. **Định hình thế Hamilton–Jacobi (tùy chọn).** Từ (7.16), $S=-ht+W$; với rô-to quay đều, $W=p\theta$, nên phần giá trị còn lại xấp xỉ $\dot\theta\,(H-t)/2\pi$ vòng. Đặt
   $$\Phi(s,t)=\zeta\,\dot\theta\,\frac{(K-k)\tau}{2\pi},\qquad F=\Phi(s',t')-\Phi(s,t),\qquad \Phi(\cdot,K)=0.$$
   Định hình theo thế không đổi chính sách tối ưu (Ng et al. 1999), vì tổng $F$ co lại thành $-\Phi(s_0,0)$.
5. **Huấn luyện.** Replay $2\times10^5$, batch 256, mạng đích mềm $\tau_{\text{soft}}=0.005$, Adam với lr giảm từ $5\times10^{-4}$ xuống $10^{-4}$, Huber, cắt chuẩn gradient 10, $\varepsilon$ giảm từ 1 xuống 0.05 trong nửa đầu. Chạy 32 môi trường song song, 30 000 vòng ($9.6\times10^5$ chuyển). Mạng $d\to128\to128\to7$, ReLU.

### 3.3 Mã giả

```text
TE-DQN(ρ, finite, use_t, ζ):
  khởi tạo Q_w, Q̄_w̄ ← Q_w, bộ nhớ D
  với mỗi vòng:
    a ← ε-tham-lam(Q_w(φ(s, k)))
    (r, s', Δk, gãy) ← EKepEnv.step(a)          # Δk = 1 + 9N nếu gãy (đợt Geom)
    k' ← k + Δk ; kết_thúc ← (k' ≥ K)
    nếu ζ > 0: r ← r + Φ(s', k')·[¬kết_thúc] − Φ(s, k)
    Γ ← 1 nếu finite, ngược lại e^{−ρ τ Δk}
    Γ ← 0 nếu kết_thúc và use_t                  # chỉ khi có t thì H mới là kết thúc thật
    lưu (φ(s,k), a, r, φ(s',k'), Γ) vào D
    lấy mẫu batch: y = r + Γ · Q̄(φ', argmax Q(φ'))
    giảm Huber(Q(φ,a) − y); Q̄ ← (1−τ_s)Q̄ + τ_s Q
```

### 3.4 Các biến thể được so sánh

| Mã | Tên | Có $t$ | $\Gamma$ | Hết horizon |
|---|---|---|---|---|
| V | DQN thường | không | $0.99241$ cố định | kết thúc |
| S | DQN bán-Markov | không | $e^{-\rho\tau\Delta k}$ | cắt (bootstrap) |
| T | TE-DQN | có | 1 | kết thúc |
| G | TE-DQN-γ | có | $e^{-\rho\tau\Delta k}$ | kết thúc |
| H | TE-DQN + Hamilton–Jacobi | có | 1 + định hình $\Phi$ | kết thúc |

Mục 4.3 giải số thêm bản thiết kế **E-DQN 4 trụ cột** (dueling, tanh, đồng hồ sin/cos, buffer bán-Markov, mặt nạ λ, chiết khấu Van der Pol).

---

## 4. Kết quả giải số

### 4.1 Chuẩn quy hoạch động và đường cơ sở

Mọi chính sách được đánh giá trên **cùng 2000 tập** của EKepEnv, dùng chung một chuỗi số ngẫu nhiên (seed 12345). Lợi tức tính bằng số vòng quay trừ giá công, trong 20 $T_0$.

| Chính sách | Lợi tức | Gãy/tập |
|---|---|---|
| **DP phụ thuộc $t$** (lưới $72\times169$, nhánh đợt Carnot) | **34.415 ± 0.025** | 0.037 |
| DP dừng, chiết khấu bán-Markov $\rho=1/\tau_{\text{dw}}$ | 34.218 ± 0.024 | 0.027 |
| $u=0.3$ hằng (tốt nhất trong các hằng số) | 28.627 ± 0.073 | 0.590 |
| Mô-men Van der Pol của bài báo ($\Omega_0=1.05\,\dot\theta_{\text{br}}$) | 21.753 ± 0.162 | 3.521 |
| $u=0$ | 1.101 ± 0.015 | 0 |

Biết $t$ chỉ đáng thêm **+0.197 vòng (0.6%)** khi lời giải là tối ưu chính xác. Đây là trần cho mọi lợi ích mà biến $t$ có thể đem lại trong bài toán này.

### 4.2 TE-DQN (mục 3): 30 000 vòng, $9.6\times10^5$ chuyển, 2 seed

| Mã | Biến thể | Lợi tức trung bình (min–max) | % DP | Gãy/tập |
|---|---|---|---|---|
| V | DQN thường ($\gamma$ cố định, không $t$) | 33.474 (33.34–33.61) | 97.3% | 0.150 |
| S | DQN bán-Markov $\gamma^{\Delta t}$ | 32.596 (32.22–32.97) | 94.7% | 0.043 |
| T | TE-DQN ($t$, $\rho=0$) | 32.325 (31.87–32.78) | 93.9% | 0.100 |
| H | TE-DQN + Hamilton–Jacobi | 30.572 (30.20–30.94) | 88.8% | 0.195 |

Seed 2 của bốn biến thể này và biến thể G không chạy được, vì file mã bị chuyển vào `rl_data/` trong lúc đang xếp hàng chạy.

**Kết quả ngược kỳ vọng, giữ nguyên:**
- Mục tiêu horizon hữu hạn **không chiết khấu** (T, H) học kém hơn DQN thường. Bootstrap phải truyền tín dụng qua 500 nhịp, và phần lợi ích tối đa mà $t$ có thể mang lại (0.6%) bị nhiễu học (khoảng 3%) lấn át.
- Định hình Hamilton–Jacobi làm **tệ hơn** (−1.75 so với T), dù về lý thuyết nó không đổi chính sách tối ưu. Thế $\Phi\propto\dot\theta\,(H-t)$ có độ lớn tới khoảng 20 vòng, nên phần dư mà mạng phải học bị nhiễu thêm.

### 4.3 Giải số bản thiết kế E-DQN 4 trụ cột (bản đề xuất của bạn)

Mã: [giai_so_e_dqn.py](giai_so_e_dqn.py), viết bằng numpy và giữ đúng thiết kế:
- tanh, đầu dueling $Q=V+A-\bar A$;
- Adam lr $10^{-3}$, Huber, cắt gradient 10;
- $\rho=0.03$, kênh đồng hồ $[T_{\text{rem}}/H,\sin,\cos]$;
- buffer $(s,a,R,\Delta t,s')$ và mặt nạ λ.

Mỗi biến thể chạy 30 000 vòng × 2 seed.

**Kiểm công thức λ trong bản đề xuất** (bằng phương trình (7) của bài báo) so với λ suy từ Lagrange loại một:

$$\lambda(s,u)=\frac{\sin\theta}{M}\big[-Im\kappa\dot\theta^2+I(k\kappa-g)-\Delta Em\kappa\cos\theta+gm\kappa^2\cos^2\theta\big]+\underbrace{\frac{m\kappa\cos\theta}{M(\theta)}\,u}_{\text{terms}(u)}$$

- Khi $u=0$: lệch tối đa $8.9\times10^{-16}$ tại $\dot\theta=1,\,2.15,\,3$.
- Khi $u=0.4$ mà bỏ số hạng $u$: lệch 0.12.
- Khi thêm $m\kappa\cos\theta\,u/M$: lệch còn $\le8.9\times10^{-16}$.

Vậy "terms(u)" chính là $m\kappa\cos\theta\,u/M(\theta)$.

**Bảng bỏ-từng-trụ:**

| Mã | Biến thể | Lợi tức trung bình (min–max) | % DP | Gãy/tập | $\langle u\rangle$ khi còn 20 / 5 / 1 / 0.4 $T_0$ |
|---|---|---|---|---|---|
| E | E-DQN đầy đủ (mặt nạ $+1\sigma$) | 33.319 (33.25–33.39) | 96.8% | 0.037 | 0.342 / 0.302 / 0.267 / **0.233** |
| E0 | mặt nạ λ tất định ($m=0$) | 33.139 (32.70–33.58) | 96.3% | 0.027 | 0.343 / 0.292 / 0.238 / 0.165 |
| Em | bỏ trụ 4 (không mặt nạ) | 32.943 (32.87–33.01) | 95.7% | 0.025 | 0.293 / 0.290 / 0.329 / 0.272 |
| Ec | bỏ sin/cos đồng hồ | 33.295 (33.15–33.44) | 96.7% | 0.022 | 0.298 / 0.352 / 0.318 / 0.216 |
| Ev | chiết khấu Van der Pol $\rho_{\text{eff}}(s)$ | 33.478 (33.00–33.95) | 97.3% | **0.231** | 0.410 / 0.403 / 0.304 / 0.169 |
| — | DP chính xác, để so | 34.415 | 100% | 0.037 | 0.40 / 0.42 / 0.44 / 0.21 (lưới, bản trước) |

**Đường học** (lợi tức kiểm định, trung bình 2 seed, theo số chuyển):

| | 80k | 240k | 400k | 560k | 720k | 880k | 959k |
|---|---|---|---|---|---|---|---|
| E | 27.9 | 29.0 | 30.4 | 31.3 | 32.2 | 33.4 | 33.0 |
| Em | 27.1 | 30.3 | 23.1 | 26.5 | 32.0 | 31.9 | 32.3 |
| Ec | 26.1 | 31.9 | 32.2 | 28.9 | 31.8 | 33.0 | 33.2 |
| Ev | 28.2 | 23.4 | 30.0 | 25.3 | 33.3 | 27.8 | 32.9 |

### 4.4 Kiểm từng tuyên bố của bản đề xuất bằng số

| Tuyên bố trong bản đề xuất | Số đo | Kết luận |
|---|---|---|
| T1: $\gamma(\Delta t)=e^{-\rho\Delta t}$, $\rho\approx1/\tau_{\text{dw}}$ thay cho $\gamma=0.99$ cảm tính | $\rho=0.03$ khớp $1/\tau_{\text{dw}}=0.0303$; DQN bán-Markov (S) 94.7%, DQN thường (V) 97.3% | **Có cơ sở vật lý**, nhưng không làm tăng điểm trong bài này |
| T1: chiết khấu Van der Pol — "tại $\dot\theta=\Omega_0$, $\gamma\to1$, nhìn xa" | $\Omega_0=3.365>\dot\theta_{\text{br}}=3.205$: pha khóa **gãy trước khi** chạm $\Omega_0$, nên $\rho_{\text{eff}}\ge0.0093$ trong $\mathcal E_{\text{past}}$. Với $\dot\theta>\Omega_0$ thì $\rho_{\text{eff}}<0$, $\gamma>1$ (1.0036 tại 3.6), Bellman phân kỳ, nên phải kẹp $\ge0$ | **Chế độ "nhìn xa" không đạt được.** Ev có lợi tức cao nhất nhưng gãy **6×** nhiều hơn E (0.231 so với 0.037) và đường học dao động mạnh nhất: liều lĩnh, kém ổn định |
| T2: kênh đồng hồ giúp liều hơn gần hạn chót và ngắt công ở $0.4\,T_0$ cuối | E: $\langle u\rangle$ 0.342 → 0.233; DP: 0.40 → 0.21 | **Đúng về hướng.** Mạng tự học giảm công cuối kỳ mà không cần lập trình tay |
| T2: thêm $\sin,\cos(2\pi t/T_0)$ | Ec (bỏ) 33.295 so với E (có) 33.319 | **Không có tác dụng**, đúng như lý thuyết: $L$ không chứa $t$ tường minh (hệ sclêrônôm, sách §1.4), động lực tự trị nên pha $t\bmod T_0$ không mang thông tin. Chỉ có ích với hệ rêônôm, $\vec r(q,t)$ |
| T3: SMDP gom chu trình, giảm số bước cập nhật | Tiết kiệm 24.0% bước với chính sách Van der Pol (3.5 gãy/tập), 5.6% với $u=0.3$; các chính sách học được gãy ≤ 0.23/tập nên tiết kiệm < 2% | **Chỉ đáng kể khi hệ gãy thường xuyên** |
| T4: mặt nạ λ là hàm chắn an toàn | E so với Em: +0.38 lợi tức, cả 2 seed của E đều cao hơn mọi seed của Em. Nhưng gãy/tập **không giảm** (0.037 so với 0.025) | Mặt nạ giúp **học nhanh hơn** (bỏ hành động tự sát khi khám phá; Em sụt còn 23.1 ở 400k), chứ **không** làm chính sách cuối an toàn hơn. Gãy còn lại do nhiễu trong nhịp, mặt nạ một nhịp không lường trước được |
| T4: áp nhiễu trước khi kiểm λ | E ($+1\sigma$) 33.319 so với E0 (tất định) 33.139; khoảng min–max chồng nhau | **Chưa phân biệt được** với 2 seed |
| Bảng 4: E-DQN tốt hơn DQN chuẩn | E 96.8% so với V 97.3% (khoảng chồng nhau) | **Không được xác nhận.** Mọi DQN đều nằm ở 95.7–97.3% DP; trần lợi ích của $t$ chỉ 0.6% nên nhiễu học quyết định thứ hạng |

### 4.5 Kết luận

1. Biến $t$ của Học thuyết E vào RL đúng chỗ qua ba đường: chiết khấu $e^{-\rho\Delta t}$ với $\rho=1/\tau_{\text{dw}}$, kênh $T_{\text{rem}}$, và chuyển bán-Markov cho đợt Carnot. Lời giải chính xác xác nhận cả ba, nhưng giá trị của thông tin $t$ chỉ **0.6%** trong bài toán này.
2. E-DQN tự học được hành vi theo hạn chót (giảm công ở $0.4\,T_0$ cuối), giống lời giải DP.
3. Trong 4 trụ cột, chỉ **mặt nạ λ** cho lợi ích đo được và nhất quán qua 2 seed (+0.38, ổn định đường học). Pha sin/cos không có tác dụng vì hệ tự trị. Chiết khấu Van der Pol làm tác tử liều lĩnh hơn (gãy 6×).
4. Thứ hạng giữa các DQN (95.7–97.3%) nằm trong nhiễu với 2 seed. Muốn tách được hiệu ứng 0.5% cần khoảng 10 seed, hoặc một bài toán mà $t$ đáng giá hơn (hệ rêônôm, hạn chót ngắn, hoặc phần thưởng phụ thuộc $t$).

### 4.6 Giới hạn

- Chỉ có 2 seed cho mỗi biến thể; chênh lệch dưới khoảng 0.5 không có ý nghĩa thống kê.
- Môi trường $\mathcal E_{\text{present}}$ dùng luật thống kê của bài báo ($\mathrm{Geom}(0.195)$, $\delta=0.36\,T_0$, tái khóa $\dot\theta=2.15$), không tích phân chi tiết va chạm Carnot. Mô phỏng đầy đủ (`giai_so_dong_co_langevin`) chậm khoảng 1000 lần, không đủ nhanh để huấn luyện DQN.
- Phần thưởng (vòng quay − 0.15·công) là lựa chọn thiết kế; DQN viết bằng numpy và chạy trên CPU.
- Sách Cơ học giải tích là bản quét: tôi OCR các trang đã trích, nên chữ trích có thể sai chính tả, nhưng số phương trình đã đối chiếu với ảnh trang.

## 5. Chạy lại

```bash
cd hoc-thuyet-e
PYTHONPATH=.:rl_data python3 rl_data/giai_so_te_dqn.py DP              # chuẩn DP + đường cơ sở (~1–2 phút)
PYTHONPATH=.:rl_data python3 rl_data/giai_so_te_dqn.py train T 0 30000 # TE-DQN: V|S|T|H|G
PYTHONPATH=.:rl_data python3 rl_data/giai_so_e_dqn.py train E 0 30000  # E-DQN: E|E0|Em|Ec|Ev (~15 phút/lượt có mặt nạ)
PYTHONPATH=.:rl_data python3 rl_data/giai_so_e_dqn.py eval             # bảng tổng hợp mọi lượt
```

## Câu hỏi còn mở

- Hệ RL thật của bạn có $t$ tường minh không (hệ rêônôm, hoặc phần thưởng hay ràng buộc đổi theo $t$)? Nếu có, kênh sin/cos và TE-DQN sẽ có tác dụng; nếu không, chỉ $T_{\text{rem}}$ và mặt nạ λ đáng giữ.
- Có cần chạy thêm seed (khoảng 10 mỗi biến thể, vài giờ CPU) để tách được các hiệu ứng nhỏ không?

