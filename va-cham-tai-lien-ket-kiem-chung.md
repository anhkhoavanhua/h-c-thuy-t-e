# Thang $\ell$, ngưỡng $R(E_{\text{vũ trụ}})$ và va chạm tái liên kết: giải và kiểm chứng số liệu

> Tiếp nối [lien-ket-dong-bo-e-gay-kiem-chung.md](lien-ket-dong-bo-e-gay-kiem-chung.md).
> Đại số kiểm bằng sympy. Số liệu tính bằng [giai_so_va_cham_tai_lien_ket.py](giai_so_va_cham_tai_lien_ket.py) (RK4, $dt=10^{-4}$, sự kiện định vị bằng chia đôi 60 lần, chỉ thư viện chuẩn, chạy `python3 giai_so_va_cham_tai_lien_ket.py` ≈ 35 s).
> Ký hiệu: $E_{\text{vt}}\equiv E_{\text{vũ trụ}}$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| Cần thang $\ell$; $\ell_0=\sqrt{\hbar/(m\omega_0)}$ | **Đúng thứ nguyên.** Proton với $\hbar\omega_0=1$ MeV cho $\ell_0=6.44\times10^{-15}$ m (cỡ fm) |
| $\ell=L_{\text{AdS}}$, $E_{\text{vt}}\sim\hbar c/L_{\text{AdS}}$ | **Đúng thứ nguyên**, nhưng là giả thuyết thêm vào, không suy ra từ cơ học. Hai lựa chọn $\ell_0$ và $L_{\text{AdS}}$ cho $\kappa$ khác nhau, nên phải **chọn một** |
| $R=E_{\text{vt}}/\ell$ | **Đúng thứ nguyên.** Hệ quả mới: $R\kappa=E_{\text{kép}}$ **chính xác**, không phụ thuộc $\ell$ và $E_{\text{vt}}$ |
| $\dot\theta_{\text{gãy}}\approx\dfrac{1}{\ell\sqrt m}\dfrac{E_{\text{vt}}}{\sqrt{E_{\text{kép}}}}$ | **Đúng đại số**, nhưng chỉ là **tiệm cận**. Sai +97% ở $E_{\text{vt}}=1$, −3.3% ở 5, −0.03% ở 20. Điều kiện dùng được: $E_{\text{vt}}\gg\ell\sqrt{kE_{\text{kép}}}$ |
| $E_{\text{vt}}$ lớn thì nền bền, nhỏ thì dễ gãy | **Đúng** (số: 0.717 → 28.29 khi $E_{\text{vt}}$ tăng 1 → 20). Với $E_{\text{vt}}$ nhỏ, gãy ngay cả khi quay rất chậm |
| Tái khóa pha tạo xung lực hữu hạn $\Delta p\neq0$ | **Đúng.** Đây là cơ chế đúng để có xung lực |
| $\Delta p=m(\kappa\dot\theta\cos\theta-v^-)$ | **Thiếu**: bỏ qua phản xung lên pha $\theta$. Đúng: $\Lambda=m_{\text{eff}}(\kappa\dot\theta^-\cos\theta-v^-)$ với $m_{\text{eff}}=\dfrac{mI}{I+m\kappa^2\cos^2\theta}$. Bản gốc lớn hơn 2.7–4.5% trong các lần mô phỏng |
| Năng lượng phát xạ $\Delta T=\tfrac12m(v^-)^2-\tfrac12m(v^+)^2$ | **SAI.** Số ra **âm** (−0.043, −0.159, −0.241) ở 3/3 lần va chạm với $I=1$, tức là "tạo ra năng lượng". Đúng là **định lý Carnot**: $\Delta E_{\text{mất}}=\tfrac12m_{\text{eff}}(\kappa\dot\theta^-\cos\theta-v^-)^2\ge0$. Số: khớp mô phỏng tuyệt đối (1.257884 = 1.257884) |
| Điều kiện tái khóa | **Bổ sung**: liên kết holonom chỉ khóa lại được khi quỹ đạo tự do **cắt lại mặt $f=0$** (vị trí khớp). Chỉ vận tốc nhảy, vị trí không nhảy |

## 1. Thang độ dài $\ell$

**Thang dao động lượng tử:**

$$\ell_0=\sqrt{\frac{\hbar}{m\omega_0}},\qquad \kappa=\ell_0\frac{E_{\text{kép}}}{E_{\text{vt}}}$$

Ví dụ số với hạt là proton:

| $\hbar\omega_0$ | $\omega_0$ (rad/s) | $\ell_0$ (m) |
|---|---|---|
| 0.1 MeV | $1.519\times10^{20}$ | $2.037\times10^{-14}$ |
| 1 MeV | $1.519\times10^{21}$ | $6.442\times10^{-15}$ |
| 10 MeV | $1.519\times10^{22}$ | $2.037\times10^{-15}$ |

**Thang AdS:** nếu $E_{\text{vt}}\sim\hbar c/L_{\text{AdS}}$ thì $\kappa\sim L_{\text{AdS}}^2E_{\text{kép}}/(\hbar c)$. Công thức này đúng thứ nguyên, nhưng việc gán $\ell=L_{\text{AdS}}$ là giả thuyết của học thuyết, không phải kết quả của cơ học giải tích.

## 2. Ngưỡng $R(E_{\text{vt}})$ và tốc độ gãy tới hạn

Với $R=E_{\text{vt}}/\ell$ và $\kappa=\ell E_{\text{kép}}/E_{\text{vt}}$ (sympy xác nhận):

$$\boxed{R\,\kappa=E_{\text{kép}}}$$

Nghĩa vật lý: công cực đại mà lực liên kết sinh ra trên biên độ đồng bộ đúng bằng $E_{\text{kép}}$.

Tốc độ gãy tiệm cận là $\dot\theta_{\text{gãy}}\approx\sqrt{R/(m\kappa)}$, tương đương:

$$\dot\theta_{\text{gãy}}\approx\frac{E_{\text{vt}}}{\ell\sqrt{mE_{\text{kép}}}}\quad\Longleftrightarrow\quad \tfrac12m\kappa^2\dot\theta_{\text{gãy}}^2\approx\tfrac12E_{\text{kép}}$$

Nói cách khác, gãy xảy ra khi động năng đồng bộ của hạt đạt cỡ $E_{\text{kép}}$.

Với $\ell=\ell_0$, công thức có dạng không thứ nguyên gọn:

$$\frac{\dot\theta_{\text{gãy}}}{\omega_0}\approx\frac{E_{\text{vt}}}{\sqrt{\hbar\omega_0\,E_{\text{kép}}}}$$

**Giới hạn hiệu lực.** Công thức bỏ qua phần thế của $\lambda$ (công thức (8) ở tài liệu trước). Phần thế có độ lớn cỡ $(k\kappa-g)$, nên chỉ bỏ qua được khi $k\kappa\ll R$, tức là:

$$E_{\text{vt}}\gg\ell\sqrt{k\,E_{\text{kép}}}$$

Nếu $k\kappa-g>R$ thì liên kết gãy **ngay cả khi quay chậm**, chỉ do thế đàn hồi.

### Số liệu

**Thông số:** $\ell=1,\ m=1,\ I=1,\ k=4,\ g=0.1,\ E_{\text{kép}}=\Delta E=0.5$. Khi đó $\ell\sqrt{kE_{\text{kép}}}=1.414$.

Tốc độ tới hạn chính xác được tìm bằng chia đôi trên $\max\lvert\lambda\rvert$ dọc quỹ đạo đồng bộ:

| $E_{\text{vt}}$ | $\kappa$ | $R$ | $R\kappa$ | $\dot\theta_{\text{gãy}}$ công thức | $\dot\theta_{\text{gãy}}$ chính xác | Sai lệch |
|---|---|---|---|---|---|---|
| 1 | 0.5000 | 1 | 0.50 | 1.41421 | 0.71749 | +97.11% |
| 2 | 0.2500 | 2 | 0.50 | 2.82843 | 3.45596 | −18.16% |
| 5 | 0.1000 | 5 | 0.50 | 7.07107 | 7.31189 | −3.29% |
| 10 | 0.0500 | 10 | 0.50 | 14.14214 | 14.22984 | −0.62% |
| 20 | 0.0250 | 20 | 0.50 | 28.28427 | 28.29303 | −0.03% |

**Kiểm chứng bằng mô phỏng đầy đủ** ($E_{\text{vt}}=5$):

| $\dot\theta(0)$ | Kết quả |
|---|---|
| 7.2388 ($0.99\times$ tới hạn) | Không gãy trên $[0,20]$ |
| 7.3850 ($1.01\times$ tới hạn) | Gãy tại $t=0.1817$ |

## 3. Va chạm tái liên kết

### 3.1. Phương trình va chạm đúng

Khi quỹ đạo tự do cắt lại mặt $f=x-\kappa\sin\theta=0$, liên kết khóa lại tức thời. Xung lực $\Lambda$ tác dụng theo $\nabla f=(-\kappa\cos\theta,\,1)$, tức là tác dụng lên **cả hai** tọa độ:

$$m(v^+-v^-)=\Lambda,\qquad I(\dot\theta^+-\dot\theta^-)=-\Lambda\kappa\cos\theta,\qquad v^+=\kappa\cos\theta\,\dot\theta^+$$

Giải hệ (sympy xác nhận):

$$\Lambda=m_{\text{eff}}\,u,\qquad u=\kappa\cos\theta\,\dot\theta^--v^-,\qquad m_{\text{eff}}=\frac{mI}{I+m\kappa^2\cos^2\theta}$$

$$\Delta E_{\text{mất}}=\tfrac12m_{\text{eff}}u^2\ \ge0\qquad\text{(định lý Carnot)}$$

**Liên hệ với bản gốc:**
- $\Lambda_{\text{bản gốc}}/\Lambda=1+m\kappa^2\cos^2\theta/I$. Hai công thức chỉ trùng nhau khi $I\to\infty$.
- Công thức $\tfrac12m(v^-)^2-\tfrac12m(v^+)^2$ sai **ngay cả khi $I\to\infty$**. Lý do: khi đó $\Delta\dot\theta\to0$, nhưng động năng pha thay đổi một lượng $I\dot\theta\,\Delta\dot\theta=-\Lambda v^+$, vẫn hữu hạn. Cộng phần này vào thì năng lượng mất đúng là $\tfrac12m(v^+-v^-)^2$.

### 3.2. Chu kỳ gãy → tự do → tái khóa, trường hợp $I=1$

**Thông số:** $I=1,\ m=1,\ k=4,\ g=0.1,\ \Delta E=0.5,\ \kappa=0.2,\ R=3$. Điều kiện đầu $\dot\theta(0)=5$. Thời gian $t\in[0,30]$. Năng lượng ban đầu $E_0=12.5$. Tốc độ tới hạn chính xác là 4.33634.

| # | Sự kiện | $t$ | $\theta$ | Chi tiết | $E$ |
|---|---|---|---|---|---|
| 1 | Gãy | 0.15384 | 44.12° | $\dot\theta=5.01413$, $\lambda=-3$ | 12.500000 |
| 2 | Khóa | 1.39258 | 34.13° | $v$: −0.77687 → +0.78795; $\dot\theta$: 5.01861 → 4.75955 | 11.242116 |
| 3 | Gãy | 1.45596 | 51.42° | $\dot\theta=4.76253$ | 11.242116 |
| 4 | Khóa | 2.68159 | 20.06° | $v$: −0.70001 → +0.84487; $\dot\theta$: 4.78742 → 4.49719 | 10.006663 |
| 5 | Gãy | 2.84937 | 63.31° | $\dot\theta=4.49485$ | 10.006663 |
| 6 | Khóa | 4.03528 | 3.53° | $v$: −0.58438 → +0.85038; $\dot\theta$: 4.54637 → 4.25996 | 8.936376 |
| — | Từ đó đến $t=30$ | | | $\dot\theta=4.19426<4.33634$ (dưới tới hạn), **đồng bộ vĩnh viễn** | 8.936376 |

Ở mọi lần khóa, $\lvert f\rvert\le2.4\times10^{-16}$, tức vị trí khớp đúng.

**So sánh xung lực và năng lượng tại các lần khóa:**

| Lần khóa | $\Lambda$ đúng | $\Lambda$ bản gốc | Tỉ số | $E$ mất (mô phỏng) | Carnot $\tfrac12m_{\text{eff}}u^2$ | $\Delta T$ bản gốc |
|---|---|---|---|---|---|---|
| 1 | +1.564818 | +1.607706 | 1.02741 | 1.257884 | 1.257884 | **−0.043381** |
| 2 | +1.544886 | +1.599410 | 1.03529 | 1.235453 | 1.235453 | **−0.159447** |
| 3 | +1.434763 | +1.491937 | 1.03985 | 1.070288 | 1.070288 | **−0.241074** |

**Tổng năng lượng mất:** 3.563624 (28.51% của $E_0$), chia thành 3 xung phát xạ rời rạc.

### 3.3. Trường hợp quán tính pha nhỏ, $I=0.05$

Thông số như 3.2 nhưng $I=0.05$, $t\in[0,10]$. Năng lượng ban đầu $E_0=0.625$.

| Sự kiện | Chi tiết |
|---|---|
| Gãy | $t=0.15142$, $\theta=44.01°$, $\dot\theta=5.19975$ |
| Khóa | $t=1.84198$, $\theta=256.35°$. Pha đã **đảo chiều** ($\dot\theta^-=-2.68562$) vì $E_0$ chỉ nhỉnh hơn đỉnh thế $\Delta E=0.5$ |
| Xung lực | $\Lambda=0.824902$ so với bản gốc $0.861648$ (tỉ số 1.04455) |
| Năng lượng mất | Mô phỏng = Carnot = **0.355388**; $\Delta T$ bản gốc = 0.262007 (thiếu 26%) |
| Cuối $t=10$ | Đồng bộ, $E=0.269612$ (mất 56.86%) |

## 4. Phát biểu đã sửa của chu kỳ Quá khứ → Hiện tại → Tương lai

1. **Pha đồng bộ** ($n=1$): hạt khóa theo $\mathcal E_{\text{kép}}$ với $\kappa=\ell E_{\text{kép}}/E_{\text{vt}}$, và lực giữ là $\lambda(\theta,\dot\theta)$.
2. **Gãy** khi $\lvert\lambda\rvert=R=E_{\text{vt}}/\ell$. Khi quay nhanh, điều này tương ứng với $\tfrac12m\kappa^2\dot\theta^2\approx\tfrac12E_{\text{kép}}$. Lúc gãy không có xung lực và không mất năng lượng.
3. **Pha tự do** ($n=2$): hạt dao động ở $\omega_0$, lệch khỏi quỹ đạo đồng bộ.
4. **Tái khóa** khi quỹ đạo cắt lại $f=0$. Xung lực hữu hạn $\Lambda=m_{\text{eff}}u$ tác dụng lên cả hạt lẫn pha. Năng lượng phát xạ là $\tfrac12m_{\text{eff}}u^2$ (Carnot), luôn dương.
5. **Hệ quả đo được:** mỗi chu kỳ mất một lượng năng lượng rời rạc, nên $\dot\theta$ giảm dần theo từng nấc (5.01 → 4.76 → 4.50 → 4.26 → 4.19). Khi $\dot\theta$ xuống dưới tốc độ tới hạn, hệ **tự ổn định ở trạng thái đồng bộ**. Chuỗi phát xạ vì vậy **hữu hạn**: 3 xung với $I=1$, 1 xung với $I=0.05$.

## Câu hỏi còn mở

- Chọn $\ell=\ell_0$ (lượng tử) hay $\ell=L_{\text{AdS}}$ (vũ trụ học)? Hai lựa chọn cho $\kappa$ khác nhau nhiều bậc độ lớn.
- Va chạm đang giả thiết khóa **hoàn toàn dẻo** (sau va chạm $\dot f=0$). Nếu có hệ số hồi phục $e>0$ thì năng lượng mất nhân thêm $(1-e^2)$, và cần thêm một quy tắc để xác định khi nào khóa lại.
- $E_{\text{kép}}$ đang lấy bằng $\Delta E$ của thế $-\Delta E\cos\theta$. Cần xác nhận đây có phải cùng một đại lượng trong học thuyết hay không.
