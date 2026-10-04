# Liên kết đồng bộ với nền $E_{\text{vũ trụ}}$ và hiện tượng "E gãy": giải và kiểm chứng số liệu

> Tiếp nối [he-e-kep-hat-vat-chat-lagrange.md](he-e-kep-hat-vat-chat-lagrange.md) (hệ (1)–(2)).
> Đại số kiểm bằng sympy. Số liệu tính bằng [giai_so_lien_ket_e_gay.py](giai_so_lien_ket_e_gay.py) (RK4, $dt=10^{-4}$, chỉ thư viện chuẩn, chạy `python3 giai_so_lien_ket_e_gay.py` ≈ 24 s).

## Tóm tắt kết luận

| Khẳng định trong bản gốc | Kết quả kiểm chứng |
|---|---|
| Phương trình Lagrăng loại một (4), (5) | **Đúng** (sympy) |
| Công thức $\lambda(t)$ (6) | **Đúng** nhưng chưa đóng, vì còn chứa $\ddot\theta$. Dạng đóng theo $(\theta,\dot\theta)$ ở mục 2.2. Số: lệch $4\times10^{-12}$ |
| $\lambda$ tăng theo $\dot\theta^2$ | **Đúng**, hệ số chính xác $-\dfrac{m\kappa I\sin\theta}{I+m\kappa^2\cos^2\theta}$; max\|λ\|/(mκθ̇²) → 1 khi quay nhanh |
| Gãy thì hệ chuyển từ $n=1$ sang $n=2$ bậc tự do | **Đúng** |
| Xung lực $\Delta p=\int_{t_g}^{t_g+\Delta t}\lambda\,dt\neq0$ khi $\Delta t\to0$ | **SAI**. Vì $\lvert\lambda\rvert\le R$ nên $\lvert\Delta p\rvert\le R\,\Delta t\to0$. Số: $\Delta p=-4.9\times10^{-4}$ ở $\Delta t=10^{-3}$, tỉ lệ tuyến tính với $\Delta t$ |
| "Năng lượng giải phóng" $\Delta E_{\text{gãy}}=\int\lambda\,dx$ | **SAI** về ý nghĩa. Năng lượng toàn phần liên tục qua điểm gãy (lệch $\sim10^{-14}$). $\int\lambda\,dx$ chỉ là năng lượng **chuyển** giữa pha và hạt, bằng đúng $-\int(-\lambda\kappa\cos\theta)\,d\theta$ |
| Gãy do cộng hưởng | **Không áp dụng** trong pha có liên kết: $x=\kappa\sin\theta$ không có bậc tự do riêng nên không có tần số $\omega_0$ để cộng hưởng |
| $\kappa\sim E_{\text{kép}}/E_{\text{vũ trụ}}$ | **Sai thứ nguyên**: $\kappa$ có thứ nguyên độ dài (vì $x=\kappa\sin\theta$), còn tỉ số năng lượng thì không thứ nguyên. Cần thêm một thang độ dài $\ell$: $\kappa=\ell\,E_{\text{kép}}/E_{\text{vũ trụ}}$ |
| Ngưỡng gãy lấy từ "liên kết không giữ" (1.3) | **Cần tách bạch**: liên kết không giữ trong sách là $f\ge0$, được giải phóng khi $\lambda$ **đổi dấu**. Ngưỡng $\lvert\lambda\rvert=R$ là một mô hình phá hủy bổ sung, không có trong (1.3) |

**Điều thực sự xảy ra khi gãy:** vị trí và xung lượng liên tục. Chỉ **gia tốc** nhảy bậc:

$$\Delta\ddot x=-\lambda/m,\qquad \Delta\ddot\theta=\lambda\kappa\cos\theta/I$$

Sau đó hạt dao động tự do ở tần số $\omega_0$ và lệch dần khỏi quỹ đạo đồng bộ.

## 1. Thiết lập

Liên kết holonom dừng:

$$f(\theta,x)=x-\kappa\sin\theta=0\qquad(3)$$

$$\frac{\partial f}{\partial\theta}=-\kappa\cos\theta,\qquad\frac{\partial f}{\partial x}=1$$

Phương trình Lagrăng loại một $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i}-\frac{\partial L}{\partial q_i}=\lambda\frac{\partial f}{\partial q_i}$, với $L$ từ tài liệu trước:

$$I_{\mathcal E}\ddot\theta+\Delta E\sin\theta-gx\cos\theta=-\lambda\kappa\cos\theta\qquad(4)$$

$$m\ddot x+kx-g\sin\theta=\lambda\qquad(5)$$

Cả hai được sympy xác nhận.

## 2. Khử liên kết

### 2.1. Phương trình rút gọn (1 bậc tự do)

Thế $x=\kappa\sin\theta$, $\dot x=\kappa\dot\theta\cos\theta$, $\ddot x=\kappa(\ddot\theta\cos\theta-\dot\theta^2\sin\theta)$ vào (5). Lấy (4) cộng với $\kappa\cos\theta\times$(5) để khử $\lambda$:

$$\big(I_{\mathcal E}+m\kappa^2\cos^2\theta\big)\ddot\theta-m\kappa^2\dot\theta^2\sin\theta\cos\theta+\Delta E\sin\theta+(k\kappa^2-2g\kappa)\sin\theta\cos\theta=0\qquad(7)$$

Phương trình này cũng là phương trình Euler–Lagrange của Lagrangian rút gọn:

$$L_{\text{rg}}=\tfrac12\big(I_{\mathcal E}+m\kappa^2\cos^2\theta\big)\dot\theta^2+\Delta E\cos\theta-\big(\tfrac12k\kappa^2-g\kappa\big)\sin^2\theta$$

Năng lượng bảo toàn trong pha có liên kết:

$$E=\tfrac12\big(I_{\mathcal E}+m\kappa^2\cos^2\theta\big)\dot\theta^2-\Delta E\cos\theta+\big(\tfrac12k\kappa^2-g\kappa\big)\sin^2\theta$$

### 2.2. Lực liên kết $\lambda$ dạng đóng

Công thức (6) của bản gốc, $\lambda=m\kappa(\ddot\theta\cos\theta-\dot\theta^2\sin\theta)+(k\kappa-g)\sin\theta$, là đúng. Thay $\ddot\theta$ từ (7) vào, ta được $\lambda$ chỉ phụ thuộc trạng thái $(\theta,\dot\theta)$:

$$\boxed{\lambda(\theta,\dot\theta)=\frac{\sin\theta\,\Big[-I_{\mathcal E}m\kappa\,\dot\theta^2+I_{\mathcal E}(k\kappa-g)-\Delta E\,m\kappa\cos\theta+g m\kappa^2\cos^2\theta\Big]}{I_{\mathcal E}+m\kappa^2\cos^2\theta}}\qquad(8)$$

- **Phần quán tính (ly tâm):** $-\dfrac{I_{\mathcal E}m\kappa\sin\theta}{I_{\mathcal E}+m\kappa^2\cos^2\theta}\,\dot\theta^2$. Đạt cực đại $m\kappa\dot\theta^2$ tại $\theta=\pm\pi/2$.
- **Phần thế:** các số hạng còn lại, không phụ thuộc $\dot\theta$.

Hệ quả: khi $\mathcal E_{\text{kép}}$ quay nhanh, $\max\lvert\lambda\rvert\approx m\kappa\dot\theta^2$. Ngưỡng gãy vì thế tương ứng với một tốc độ quay tới hạn:

$$\dot\theta_{\text{gãy}}\approx\sqrt{R/(m\kappa)}$$

## 3. Tại thời điểm gãy

Liên kết lý tưởng chỉ tác dụng lực hữu hạn $\lvert\lambda\rvert\le R$. Khi bỏ liên kết:

1. **Vị trí và vận tốc liên tục.** Do đó xung lượng không nhảy:
   $$\lvert\Delta p\rvert=\Big\lvert\int_{t_g}^{t_g+\Delta t}F\,dt\Big\rvert\le F_{\max}\,\Delta t\xrightarrow{\Delta t\to0}0$$
   Xung lực va chạm (phương trình va chạm Đalămbe) chỉ xuất hiện khi **đặt vào** một liên kết đột ngột (va đập), không xuất hiện khi **gỡ bỏ** liên kết.
2. **Gia tốc nhảy bậc** vì lực $\lambda$ biến mất:
   $$\Delta\ddot x=-\lambda/m,\qquad\Delta\ddot\theta=+\lambda\kappa\cos\theta/I_{\mathcal E}$$
3. **Năng lượng toàn phần liên tục.** Liên kết lý tưởng không sinh công tổng: $\lambda\big(\tfrac{\partial f}{\partial\theta}\dot\theta+\tfrac{\partial f}{\partial x}\dot x\big)=\lambda\dot f=0$. Vậy $\int\lambda\,dx$ (công lên hạt) $=-\int(-\lambda\kappa\cos\theta)\,d\theta$ (công lên pha). Đây là năng lượng **trao đổi**, không phải năng lượng giải phóng.
4. **Dấu hiệu quan sát được của "gãy":** hạt thoát khỏi quỹ đạo đồng bộ. Nó dao động tự do với biên độ $\sqrt{x_g^2+(\dot x_g/\omega_0)^2}$, nên $\lvert x-\kappa\sin\theta\rvert$ tăng lên cỡ biên độ đó cộng $\kappa$.

## 4. Kết quả tính số

**Thông số chung:** $I_{\mathcal E}=1,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.5,\ \kappa=0.2$ (không thứ nguyên).
**Điều kiện đầu:** $\theta(0)=0$, $\dot\theta(0)=\omega$, $x(0)=0$, $\dot x(0)=\kappa\omega$ (thỏa $f=\dot f=0$).

**Phương pháp:** pha có liên kết được giải như hệ DAE 4 biến. Mỗi bước giải hệ tuyến tính 3×3 gồm (4), (5) và $\ddot f=0$ để lấy $\ddot\theta,\ddot x,\lambda$, **độc lập** với (6) và (8). Khi $\lvert\lambda\rvert$ chạm $R$, thời điểm gãy được định vị bằng chia đôi bước (60 lần). Sau đó hệ chuyển sang (1)–(2).

### 4.1. Kiểm chứng pha có liên kết ($R=\infty$, $\omega=5$, $t\in[0,20]$)

| Kiểm tra | Kết quả |
|---|---|
| $\max\lvert\lambda_{\text{DAE}}-\lambda_{(6)}\rvert$ | $4.18\times10^{-12}$ |
| $\max\lvert\lambda_{\text{DAE}}-\lambda_{(8)}\rvert$ | $4.14\times10^{-12}$ |
| Trôi liên kết $\max\lvert x-\kappa\sin\theta\rvert$ | $1.05\times10^{-12}$ |
| Trôi năng lượng tương đối | $9.28\times10^{-14}$ |
| $\max\lvert\lambda\rvert$ | 4.28500 |
| $\int\lambda\,dx$ (lên hạt) / $\int(-\lambda\kappa\cos\theta)\,d\theta$ (lên pha) | −0.424708 / +0.424708, tổng $-8.9\times10^{-9}$ |

### 4.2. $\max\lvert\lambda\rvert$ theo tốc độ quay

Tính từ (8), với $\dot\theta(\theta)$ lấy từ bảo toàn năng lượng rút gọn:

| $\dot\theta(0)$ | $\max\lvert\lambda\rvert$ | $m\kappa\dot\theta(0)^2$ | Tỉ số |
|---|---|---|---|
| 2 | 0.2146 | 0.8 | 0.268 |
| 3 | 0.9853 | 1.8 | 0.547 |
| 5 | 4.2850 | 5.0 | 0.857 |
| 8 | 12.3911 | 12.8 | 0.968 |
| 12 | 29.0293 | 28.8 | 1.008 |

Giá trị ở $\dot\theta(0)=3$ và $5$ trùng với max|λ| đo trực tiếp từ mô phỏng (0.9853 và 4.28500). Tỉ số tiến về 1 khi quay nhanh, xác nhận $\lambda\propto\dot\theta^2$.

Với $R=3$, tốc độ tới hạn ước lượng $\sqrt{R/(m\kappa)}=3.87$.

### 4.3. Mô phỏng E gãy, ngưỡng $R=3$

| | $\omega=3$ | $\omega=5$ | $\omega=8$ |
|---|---|---|---|
| Có gãy? | **Không** (max\|λ\| = 0.9853 < 3) | **Có** | **Có** |
| $t_{\text{gãy}}$ | — | 0.153839 | 0.032187 |
| $\theta_{\text{gãy}}$ | — | 0.770071 rad (44.12°) | 0.257581 rad (14.76°) |
| $\dot\theta_{\text{gãy}}$ | — | 5.01413 | 8.00755 |
| $\lambda_{\text{gãy}}$ | — | −3.000000 | −3.000000 |
| $x_{\text{gãy}}$ (= $\kappa\sin\theta$) | — | +0.139237 | +0.050948 |
| $p_{\text{gãy}}=m\dot x$ (liên tục) | — | +0.719889 | +1.548674 |
| $\ddot x$: liên kết → tự do | — | −3.487330 → −0.487330 | −3.178319 → −0.178319 |
| $\Delta\ddot x$ / lý thuyết $-\lambda/m$ | — | +3.000000 / +3.000000 | +3.000000 / +3.000000 |
| $\ddot\theta$: liên kết → tự do | — | +0.092619 → −0.338098 | +0.457761 → −0.122444 |
| $\Delta\ddot\theta$ / lý thuyết $\lambda\kappa\cos\theta/I$ | — | −0.430717 / −0.430717 | −0.580205 / −0.580205 |
| $E$ đầu / lúc gãy | — | 12.5000000000 / 12.5000000000 | 32.7800000000 / 32.7800000000 |
| Trôi $E$ sau gãy (10 đv t) | — | $5.4\times10^{-14}$ | $2.8\times10^{-14}$ |
| $\int\lambda\,dx$ tới lúc gãy (lên hạt) | — | −0.206952 | −0.076261 |
| $\int(-\lambda\kappa\cos\theta)\,d\theta$ (lên pha) | — | +0.206952 | +0.076261 |
| $\tfrac12m\dot x^2+\tfrac12kx^2$ lúc gãy → max sau đó | — | 0.297894 → 0.331834 | 1.204388 → 1.245892 |
| $\max\lvert x-\kappa\sin\theta\rvert$ sau gãy | — | 0.598157 (≈3κ) | 0.957587 (≈4.8κ) |

**Xung lượng ngay sau gãy**, $\Delta p=p(t_g+\Delta t)-p(t_g)$:

| $\Delta t$ | $\omega=5$ | $\omega=8$ | Chặn trên $R\,\Delta t$ |
|---|---|---|---|
| 0.1 | $-6.129\times10^{-2}$ | $-4.519\times10^{-2}$ | $3\times10^{-1}$ |
| 0.01 | $-4.999\times10^{-3}$ | $-2.054\times10^{-3}$ | $3\times10^{-2}$ |
| 0.001 | $-4.886\times10^{-4}$ | $-1.810\times10^{-4}$ | $3\times10^{-3}$ |

$\Delta p$ giảm tuyến tính theo $\Delta t$, với hệ số bằng đúng gia tốc tự do $m\ddot x_{\text{tự do}}$ (−0.487 và −0.178). Nghĩa là **không có xung lực**: $\lim_{\Delta t\to0}\Delta p=0$.

### 4.4. Kiểm tra chéo biên độ sau gãy ($\omega=5$)

Biên độ dao động tự do dự đoán:
$$\sqrt{x_g^2+(\dot x_g/\omega_0)^2}=\sqrt{0.139237^2+(0.719889/2)^2}=0.3860$$

Độ lệch cực đại dự đoán ≈ 0.386 + κ = 0.586, còn đo được là 0.598. Phần chênh nhỏ đến từ lực cưỡng bức $g\sin\theta$ và sự trao đổi năng lượng qua $g$: năng lượng riêng của hạt tăng từ 0.2979 lên 0.3318.

## 5. Phát biểu đã sửa của mô hình "E gãy"

1. Trước khi gãy, hạt bị liên kết (3) giữ đồng bộ. Hệ có **1 bậc tự do**, chuyển động theo (7), và lực liên kết là (8).
2. $\lvert\lambda\rvert$ tăng theo $\dot\theta^2$. Với ngưỡng bền $R$, liên kết gãy khi $\mathcal E_{\text{kép}}$ quay nhanh hơn khoảng $\sqrt{R/(m\kappa)}$, tại góc mà $\lvert\lambda(\theta,\dot\theta)\rvert=R$ lần đầu.
3. Lúc gãy, **không có xung lực và không có năng lượng tỏa ra**. Vị trí, xung lượng và năng lượng toàn phần đều liên tục. Gia tốc hạt nhảy bậc đúng $-\lambda/m$ (độ lớn $R/m$).
4. Sau gãy, hệ có **2 bậc tự do** theo (1)–(2). Hạt mang theo động năng đã có sẵn ($\tfrac12m\dot x_g^2$), dao động tự do ở $\omega_0$ và lệch khỏi quỹ đạo đồng bộ nhiều lần $\kappa$. "Năng lượng của cú gãy" là phần năng lượng dao động riêng này. Nó đã tích lũy **dần** qua công $\int\lambda\,dx$ trong pha đồng bộ, không phát ra tức thời.
5. Nếu muốn mô hình có xung lực thật, cần một cơ chế **đặt vào** liên kết đột ngột (va chạm, tái đồng bộ). Khi đó dùng phương trình va chạm, với $\Delta p$ hữu hạn trong $\Delta t\to0$.

## Câu hỏi còn mở

- Thang độ dài $\ell$ trong $\kappa=\ell\,E_{\text{kép}}/E_{\text{vũ trụ}}$ lấy từ đâu?
- $R_{\text{ngưỡng}}(E_{\text{vũ trụ}})$ phụ thuộc $E_{\text{vũ trụ}}$ theo dạng hàm nào? Hiện đang cho $R$ là hằng số tự chọn ($R=3$).
