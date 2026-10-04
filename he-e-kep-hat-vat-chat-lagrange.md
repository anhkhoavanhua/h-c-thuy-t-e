# Hệ $\mathcal{E}_{\text{kép}}$ – Hạt vật chất: lập và giải phương trình Lagrăng loại hai

> Bản đã sửa. Phương trình (1)–(2), nghiệm cưỡng bức, nghiệm cộng hưởng và mode chuẩn đã kiểm bằng sympy (phần dư = 0).
> Phương pháp theo sách Cơ học giải tích, Chương III.

## Bước 1: Tọa độ suy rộng và động năng $T$

Hệ có hai bậc tự do:

- **Động cơ $\mathcal{E}_{\text{kép}}$** (quay hoặc dao động pha): tọa độ góc $q_1=\theta$, vận tốc suy rộng $\dot\theta$.
- **Hạt vật chất**: tọa độ tịnh tiến $q_2=x$, vận tốc suy rộng $\dot x$.

Theo (1.24) và (3.19), động năng là hàm toàn phương của các vận tốc suy rộng:

$$T=\tfrac12 I_{\mathcal E}\dot\theta^2+\tfrac12 m\dot x^2$$

- $I_{\mathcal E}$: quán tính pha của $\mathcal{E}_{\text{kép}}$, thứ nguyên [kg·m²].
- $m$: khối lượng của hạt.

## Bước 2: Thế năng $\Pi$ và hàm Lagrăng $L$

Thế năng gồm ba thành phần:

1. **Thế nội tại của $\mathcal{E}_{\text{kép}}$**, dạng con lắc (Mục 3.2–3.3):
   $$\Pi_{\mathcal E}(\theta)=-\Delta E\cos\theta,\qquad \Delta E=E_{\text{future}}-E_{\text{past}}$$
   Giả thiết $\Delta E>0$, khi đó cực tiểu nằm tại $\theta=0$. Nếu $\Delta E<0$ thì cực tiểu dời về $\theta=\pi$.
2. **Thế liên kết của hạt:** $V(x)=\tfrac12kx^2$
3. **Thế tương tác ghép cặp:** $\Pi_{\text{int}}=-gx\sin\theta$, trong đó $g$ có thứ nguyên [J/m] (tức là lực).

Tổng thế năng:

$$\Pi(\theta,x)=-\Delta E\cos\theta+\tfrac12kx^2-gx\sin\theta$$

Theo (3.23), $L=T-\Pi$:

$$L=\tfrac12 I_{\mathcal E}\dot\theta^2+\tfrac12 m\dot x^2+\Delta E\cos\theta-\tfrac12kx^2+gx\sin\theta$$

## Bước 3: Phương trình Lagrăng loại hai

$$\frac{d}{dt}\Big(\frac{\partial L}{\partial\dot q_i}\Big)-\frac{\partial L}{\partial q_i}=0\quad(i=1,2)$$

| | $\partial L/\partial\dot q$ | $\frac{d}{dt}(\cdot)$ | $\partial L/\partial q$ |
|---|---|---|---|
| $q_1=\theta$ | $I_{\mathcal E}\dot\theta$ | $I_{\mathcal E}\ddot\theta$ | $-\Delta E\sin\theta+gx\cos\theta$ |
| $q_2=x$ | $m\dot x$ | $m\ddot x$ | $-kx+g\sin\theta$ |

$$I_{\mathcal E}\ddot\theta+\Delta E\sin\theta-gx\cos\theta=0\qquad(1)$$

$$m\ddot x+kx=g\sin\theta\qquad(2)$$

**Tích phân năng lượng.** $L$ không chứa $t$ tường minh, nên năng lượng toàn phần được bảo toàn:

$$E=\tfrac12 I_{\mathcal E}\dot\theta^2+\tfrac12 m\dot x^2-\Delta E\cos\theta+\tfrac12kx^2-gx\sin\theta=\text{const}$$

**Vị trí cân bằng.** Cho $\partial\Pi/\partial\theta=\partial\Pi/\partial x=0$. Ta thấy $(\theta,x)=(0,0)$ là một vị trí cân bằng. Nghĩa là nếu $\mathcal E_{\text{kép}}$ không chuyển động ($\theta\equiv0$) thì hạt có thể đứng yên. Hạt chỉ bị kéo chuyển động khi pha $\theta(t)$ biến thiên.

## Bước 4: Khảo sát và giải

Phương trình (2) cho thấy hạt chịu một lực cưỡng bức $g\sin\theta(t)$ do pha của $\mathcal E_{\text{kép}}$ sinh ra. Với $\theta(t)$ bất kỳ, nghiệm tổng quát có dạng Duhamel:

$$x(t)=A\cos(\omega_0t+\alpha)+\frac{g}{m\omega_0}\int_0^t\sin\omega_0(t-\tau)\,\sin\theta(\tau)\,d\tau,\qquad \omega_0=\sqrt{k/m}$$

### Trường hợp A: $\mathcal E_{\text{kép}}$ quay gần đều, $\theta(t)\approx\Omega t$

**Điều kiện áp dụng** (so sánh các đại lượng cùng thứ nguyên năng lượng):

$$\tfrac12 I_{\mathcal E}\Omega^2\gg\Delta E\quad\text{và}\quad \tfrac12 I_{\mathcal E}\Omega^2\gg g|x|_{\max}$$

Khi đó $\theta$ quay tròn chứ không lắc, vì năng lượng đã vượt đường phân cách $2\Delta E$. Độ lệch của $\dot\theta$ so với $\Omega$ chỉ cỡ $\delta\Omega/\Omega\sim\Delta E/(I_{\mathcal E}\Omega^2)\ll1$.

**$\Omega$ được xác định bởi điều kiện ban đầu**, cụ thể $\Omega\approx\sqrt{2E/I_{\mathcal E}}$. Nó không được xác định bởi $\Delta E$. Nếu muốn gắn thêm $\Omega=\Delta E/\hbar$ thì đó là một **tiên đề lượng tử bổ sung** (quan hệ Bohr), không suy ra được từ $L$ cổ điển này.

Với $\theta\approx\Omega t$, phương trình (2) trở thành:

$$\ddot x+\omega_0^2x=\frac gm\sin\Omega t$$

**(a) Không cộng hưởng ($\Omega\neq\omega_0$):**

$$x(t)=A\cos(\omega_0t+\alpha)+\frac{g}{m(\omega_0^2-\Omega^2)}\sin\Omega t$$

$$v(t)=-A\omega_0\sin(\omega_0t+\alpha)+\frac{g\Omega}{m(\omega_0^2-\Omega^2)}\cos\Omega t$$

$$p(t)=-mA\omega_0\sin(\omega_0t+\alpha)+\frac{g\Omega}{\omega_0^2-\Omega^2}\cos\Omega t$$

**(b) Cộng hưởng ($\Omega=\omega_0$):** biên độ tăng tuyến tính theo thời gian:

$$x(t)=A\cos(\omega_0t+\alpha)-\frac{g}{2m\omega_0}\,t\cos\omega_0t$$

Nghiệm này chỉ đúng trong một khoảng thời gian hữu hạn. Khi $g|x|$ lớn dần, điều kiện của trường hợp A bị vi phạm, phương trình (1) phản ứng ngược lên $\theta$ và làm lệch cộng hưởng.

### Trường hợp B: dao động nhỏ quanh $(0,0)$

Tuyến tính hóa với $\sin\theta\approx\theta$, $\cos\theta\approx1$:

$$I_{\mathcal E}\ddot\theta+\Delta E\,\theta-gx=0,\qquad m\ddot x+kx-g\theta=0$$

Tìm nghiệm dạng $e^{i\omega t}$, ta được phương trình đặc trưng:

$$(\Delta E-I_{\mathcal E}\omega^2)(k-m\omega^2)=g^2$$

Đặt $\omega_\theta^2=\Delta E/I_{\mathcal E}$, hai tần số mode chuẩn là:

$$\omega_\pm^2=\tfrac12\Big[\omega_\theta^2+\omega_0^2\pm\sqrt{(\omega_\theta^2-\omega_0^2)^2+4g^2/(I_{\mathcal E}m)}\Big]$$

- Ghép cặp $g$ đẩy hai tần số ra xa nhau. Khi $\omega_\theta=\omega_0$, khoảng tách là $\omega_+^2-\omega_-^2=2g/\sqrt{I_{\mathcal E}m}$.
- Cân bằng $(0,0)$ chỉ **ổn định khi $\Delta E\cdot k>g^2$**, tức ma trận Hessian của $\Pi$ xác định dương. Nếu ghép cặp mạnh hơn mức này, hệ rời khỏi $(0,0)$.
- Năng lượng trao đổi qua lại tuần hoàn giữa pha $\theta$ và hạt $x$, với chu kỳ phách $2\pi/(\omega_+-\omega_-)$.

## Kết luận

1. Hạt chuyển động khi và chỉ khi pha $\theta(t)$ của $\mathcal E_{\text{kép}}$ biến thiên và $g\neq0$. Khi $\theta\equiv0$, hạt có thể đứng yên tại cân bằng $(0,0)$.
2. Xung lượng $p(t)$ **dao động** chứ không luôn khác 0. Nó triệt tiêu tại các thời điểm rời rạc, ví dụ với $A=0$ là tại $\Omega t=\pi/2+n\pi$, lúc hạt ở biên. Vì $x$ bị chặn nên trung bình theo thời gian $\langle p\rangle=0$.
3. Biên độ của phần xung lượng cưỡng bức là $\dfrac{g\Omega}{|\omega_0^2-\Omega^2|}$. Nó phụ thuộc vào $\Omega$ và tăng mạnh gần cộng hưởng $\Omega\to\omega_0$.
4. Việc đồng nhất $\Omega=\Delta E/\hbar$ là một giả thuyết cần kiểm chứng riêng, không phải hệ quả của phương trình Lagrăng.

## Bước 5: Thông số và kết quả tính số

Hệ (1)–(2) được tích phân trực tiếp bằng RK4 (bước $dt=10^{-3}$, trên $[0,20]$) rồi so với các công thức giải tích ở Bước 4. Mã nguồn: [giai_so_he_e_kep.py](giai_so_he_e_kep.py), chỉ dùng thư viện chuẩn, chạy bằng `python3 giai_so_he_e_kep.py` (< 1 giây). Đơn vị không thứ nguyên.

Trong mọi trường hợp, năng lượng $E$ được bảo toàn tới sai số tương đối $\le 3.4\times10^{-14}$, xác nhận tích phân năng lượng ở Bước 3.

### 5.1. Trường hợp A(a): quay gần đều, không cộng hưởng

**Thông số:** $I_{\mathcal E}=1,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.5$. Điều kiện đầu: $\theta(0)=0,\ \dot\theta(0)=5,\ x(0)=0,\ \dot x(0)=C\bar\Omega$ (chọn để $A=0$).

| Đại lượng | Công thức | Giá trị tính | Giá trị đo (số) |
|---|---|---|---|
| $\omega_0$ | $\sqrt{k/m}$ | 2.000000 | — |
| Năng lượng quay | $E_{\text{rot}}=\tfrac12I\dot\theta(0)^2-\Delta E$ | 12.0000 ($=24\,\Delta E$) | — |
| Tốc độ quay xấp xỉ | $\sqrt{2E_{\text{rot}}/I}$ | 4.898979 | — |
| Tốc độ quay trung bình | $\bar\Omega=2\pi/\oint d\theta/\dot\theta$ | 4.897384 | $\theta(20)/20=4.896809$ |
| Gợn tốc độ | $\delta\Omega/\Omega\sim\Delta E/(I\bar\Omega^2)$ | 0.0208 | 0.0208 ($\dot\theta\in[4.79583,\,5]$) |
| Biên độ $x$ cưỡng bức | $C=g/(m(\omega_0^2-\bar\Omega^2))$ | −0.0050039 | sai lệch max $6.73\times10^{-5}$ (1.35%) |
| Biên độ $p$ cưỡng bức | $g\bar\Omega/(\omega_0^2-\bar\Omega^2)$ | −0.0245061 | — |
| $\max\lvert\theta-\bar\Omega t\rvert$ | — | — | 0.0213 rad (gợn tuần hoàn, không trôi) |
| Trung bình $\langle p\rangle$ trên $[0,20]$ | $=0$ | 0 | $1.34\times10^{-4}$ (≈0.5% biên độ) |

**Các thời điểm $p(t)=0$** (kiểm chứng kết luận 2):

| $n$ | Giải tích $t=(\pi/2+n\pi)/\bar\Omega$ | Số |
|---|---|---|
| 0 | 0.3207 | 0.3200 |
| 1 | 0.9622 | 0.9630 |
| 2 | 1.6037 | 1.6040 |
| 3 | 2.2452 | 2.2460 |

> **Lưu ý khi dùng công thức:** $\Omega$ trong nghiệm phải là tốc độ quay **trung bình** $\bar\Omega$, không phải $\dot\theta(0)$. Nếu thế $\Omega=\dot\theta(0)=5$, pha trôi $\approx2.06$ rad sau $t=20$ và sai số của $x$ lên tới 173%. Dùng $\bar\Omega=4.897384$ thì sai số giảm còn 1.35%. Phần còn lại đến từ gợn 2% của $\dot\theta$ do thế $-\Delta E\cos\theta$.

### 5.2. Trường hợp A(b): cộng hưởng $\Omega=\omega_0$

**Thông số:** $I_{\mathcal E}=100,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.01$ (quán tính lớn để $\theta$ quay đều), $\dot\theta(0)=\omega_0=2$. Điều kiện đầu $x(0)=0,\ \dot x(0)=-g/(2m\omega_0)=-0.025$.

Tốc độ tăng biên độ: $g/(2m\omega_0)=0.025$ trên một đơn vị $t$.

| Đỉnh $t_n=n\pi/\omega_0$ | $\lvert x\rvert$ số | $g\,t_n/(2m\omega_0)$ |
|---|---|---|
| 4.7124 | 0.11780 | 0.11781 |
| 9.4248 | 0.23564 | 0.23562 |
| 14.1372 | 0.35346 | 0.35343 |
| 18.8496 | 0.47134 | 0.47124 |

Sai lệch max so với $-\frac{g}{2m\omega_0}t\cos\omega_0t$ trên $[0,20]$: $2.17\times10^{-3}$ (≈0.5% biên độ cuối).

### 5.3. Trường hợp B: dao động nhỏ, ổn định

**Thông số:** $I_{\mathcal E}=1,\ m=1,\ k=1,\ \Delta E=1,\ g=0.2$, suy ra $\omega_\theta=\omega_0=1$. Điều kiện đầu: $\theta(0)=0.01$, các biến còn lại bằng 0.

| Đại lượng | Công thức | Giá trị |
|---|---|---|
| Điều kiện ổn định | $\Delta E\cdot k=1 > g^2=0.04$ | ổn định |
| $\omega_+^2,\ \omega_-^2$ | $1\pm g$ (vì $\omega_\theta=\omega_0$) | 1.200000, 0.800000 |
| $\omega_+,\ \omega_-$ | — | 1.095445, 0.894427 |
| Chu kỳ phách | $2\pi/(\omega_+-\omega_-)$ | 31.2568 |
| Thời điểm năng lượng chuyển hết sang $x$ | dự đoán $T_{\text{phách}}/2=15.628$ | đo: 15.628 |
| Phần năng lượng chuyển sang $x$ | — | 99.0% |

### 5.4. Trường hợp B: ghép cặp mạnh, mất ổn định

**Thông số:** như 5.3 nhưng $g=1.2$. Điều kiện đầu $\theta(0)=10^{-6}$.

| Đại lượng | Công thức | Giá trị |
|---|---|---|
| Điều kiện ổn định | $\Delta E\cdot k=1 < g^2=1.44$ | **không** ổn định |
| $\omega_+^2,\ \omega_-^2$ | — | 2.2, **−0.2** |
| Tốc độ tăng | $\lambda=\sqrt{-\omega_-^2}$ | 0.44721 |
| Tốc độ tăng đo được | $\ln(\theta(20)/\theta(10))/10$ | 0.44863 (lệch 0.3%) |
| $\max\lvert\theta\rvert$ trên $[0,20]$ | — | 0.002 rad (tăng ~2000 lần từ $10^{-6}$) |

## Phụ lục: những chỗ đã sửa so với bản gốc

| Bản gốc | Vấn đề | Đã sửa thành |
|---|---|---|
| Điều kiện $I_{\mathcal E}\gg g$ | Sai thứ nguyên ([kg·m²] so với [J/m]) | $\tfrac12 I_{\mathcal E}\Omega^2\gg\Delta E,\ g\lvert x\rvert_{\max}$ |
| "$p(t)$ luôn khác 0 tại mọi $t$" | Sai: $p$ là tổng các sin, có nghiệm $p=0$; $(0,0)$ là cân bằng | $p(t)$ dao động, $\langle p\rangle=0$ |
| $\Omega=\Delta E/\hbar$ | $\hbar$ không có trong $L$; $\Omega$ do điều kiện ban đầu | Ghi rõ là tiên đề lượng tử bổ sung |
| Thiếu $\Omega=\omega_0$ | Công thức phân kỳ | Thêm nghiệm cộng hưởng $-\frac{g}{2m\omega_0}t\cos\omega_0t$ |
| Nhắc "dao động nhỏ" nhưng không giải | Thiếu | Thêm mode chuẩn $\omega_\pm$ và điều kiện ổn định $\Delta E\,k>g^2$ |
| — | Thiếu | Thêm tích phân năng lượng, vị trí cân bằng, nghiệm Duhamel |
| $\Omega$ trong nghiệm cưỡng bức | Không nói rõ là $\dot\theta(0)$ hay trung bình | $\Omega=\bar\Omega$ (trung bình); dùng $\dot\theta(0)$ gây trôi pha (Bước 5.1) |
| — | Không có số liệu kiểm chứng | Thêm Bước 5: giải số RK4 cho 4 trường hợp |
