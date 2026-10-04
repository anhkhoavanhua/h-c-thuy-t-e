# Bộ thông số chuẩn hóa: $I=\chi\hbar/\omega_0$, thang QCD $L\approx1$ fm, cửa sổ $T\sim10T_0$ — giải và kiểm chứng số liệu

> Tiếp nối [chi-va-dai-hoi-phuc-kiem-chung.md](chi-va-dai-hoi-phuc-kiem-chung.md).
> Số liệu tính bằng [giai_so_cong_huong_qcd.py](giai_so_cong_huong_qcd.py). Chạy `python3 giai_so_cong_huong_qcd.py A B` (≈ 8 s) và `… C` (vài phút).
> Hằng số: $\hbar c=197.327$ MeV·fm, $m_pc^2=938.272$ MeV.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $I=m\Delta E/k=\Delta E/\omega_0^2=\chi\hbar/\omega_0$ | **Đúng đại số và thứ nguyên**. Với $\hbar\omega_0=1$ MeV, $\chi=1$: $I=6.94\times10^{-56}$ kg·m² |
| Cộng hưởng khiến liên kết "bền vững nhất", "chỉ gãy bởi ly tâm" | **Chỉ đúng một phần.** Góc giữ tĩnh tăng khi gần biên (21.8° → 62.1° ở $\chi=0.1$). Nhưng ở $\varepsilon/\sqrt\chi=0.5$ liên kết **vẫn gãy tĩnh ở 62°**, và có trường hợp $\dot\theta_{\text{gãy}}$ **thấp hơn** so với $I=1$ (1.225 so với 1.436). Ở thang QCD ($\varepsilon/\sqrt\chi\approx197$), chọn $I$ nào cũng **không ảnh hưởng**: phần tĩnh đã nhỏ cỡ $10^{-5}$ |
| Hệ quả mới của cộng hưởng | Với $g=0$, động lực học trở nên **không phụ thuộc $\chi$** (các hàng $\chi=0.1,1,10$ trùng nhau tới 5 chữ số) |
| Loại hướng B, chọn hướng A (giam hãm QCD) | **Hợp lý** nếu $\mathcal E_{\text{kép}}$ ghép với vật chất: lúc đó chiều phụ không thể chỉ chứa hấp dẫn. $L=\hbar c/200\text{ MeV}=0.987$ fm **đúng** |
| $L=1$ fm, $\chi=1$: $k\kappa/R=2.6\times10^{-5}$, $\dot\theta_{\text{gãy}}/\omega_0=197$ | **Đúng số**, nhưng kéo theo hệ quả quan trọng: **năng lượng quay cần để gãy là 19.5 GeV = 20.8 $m_pc^2$**. Với $\hbar\omega_0$ thực của nucleon (7–16 MeV), vẫn cần 1.2–2.9 GeV = 1.3–3.1 $m_pc^2$. Mức này vượt xa năng lượng kích thích hạt nhân (cỡ MeV) và nằm ngoài miền phi tương đối tính của mô hình |
| $T_0=2\pi/\omega_0\approx4\times10^{-21}$ s; $T_{\text{QCD}}=L/c\approx3\times10^{-24}$ s | **Đúng số** ($4.14\times10^{-21}$ s và $3.3\times10^{-24}$ s). Nhưng nếu decoherence theo $T_{\text{QCD}}$ thì cửa sổ chỉ là $\approx8\times10^{-4}\,T_0$, ngắn hơn cả một dao động, nên không $e$ nào kịp khóa. Lựa chọn $(3$–$10)T_0$ chưa có căn cứ |
| "Đơn vị $m=k=1$, $T_0=2\pi$" áp cho bảng 3.1/3.2 | **Sai.** Bảng đó dùng $k=4$, tức $\omega_0=2$ và $T_0=\pi$. Vậy $10T_0=31.4$, không phải 62.8 |
| Trong $T\sim10T_0$, dải khóa là $0\le e\le0.3$ | **Sai với chính bảng cũ:** $e=0.3$ khóa ở $t=46.6$ và $54.6$, đều $>31.4$. Đúng là $e^*(10T_0)\approx0.2$. Ở $3T_0=9.4$ chỉ còn $e=0$ (và $e=0.1$ khi $\dot\theta_0=6$) |
| "$e\ge0.4$: $t>200$ tức hơn 30 chu kỳ" | **Sai hệ số:** $200/\pi=63.7$ chu kỳ |
| Kết luận chuẩn hóa "$e\le0.3$ tái khóa trong $10T_0$" áp cho bộ thông số QCD cộng hưởng | **Bị bác bỏ bằng mô phỏng.** Ở đúng bộ thông số đề xuất ($\chi=1$, $I=\chi$, $\varepsilon=197.3$), trong $10T_0$ **không $e$ nào tái đồng bộ, kể cả $e=0$**: 20 lần khóa đều gãy lại ngay. Cần ≈3 120 va chạm ($\approx1\,560\,T_0$) mới về được đồng bộ |

## 1. Đơn vị thật (proton, $\chi=1$)

Các đại lượng trong bảng:
- $\hbar\omega_0^{\text{thực}}\approx41A^{-1/3}$ MeV (mô hình vỏ hạt nhân)
- $E_{\text{quay}}^{\text{tới hạn}}=\tfrac12I\dot\theta_{\text{gãy}}^2\approx E_{\text{vt}}^2/(2\hbar\omega_0)$, không phụ thuộc $\chi$
- $v/c$ của hạt lúc gãy $=\sqrt{\chi\hbar\omega_0/m_pc^2}$

| Nguồn $\hbar\omega_0$ | $\hbar\omega_0$ (MeV) | $E_{\text{vt}}$ (MeV) | $I$ (kg·m²) | $k\kappa/R$ | $\dot\theta_{\text{gãy}}/\omega_0$ | $E_{\text{quay}}^{\text{tới hạn}}$ | $/m_pc^2$ | $v/c$ |
|---|---|---|---|---|---|---|---|---|
| Giả thiết cũ | 1.000 | 197.33 | $6.94\times10^{-56}$ | $2.57\times10^{-5}$ | 197.33 | 19.47 GeV | **20.75** | 0.033 |
| $A=16$ | 16.271 | 197.33 | $4.27\times10^{-57}$ | $6.80\times10^{-3}$ | 12.13 | 1.197 GeV | **1.28** | 0.132 |
| $A=40$ | 11.988 | 197.33 | $5.79\times10^{-57}$ | $3.69\times10^{-3}$ | 16.46 | 1.624 GeV | **1.73** | 0.113 |
| $A=208$ | 6.920 | 197.33 | $1.00\times10^{-56}$ | $1.23\times10^{-3}$ | 28.52 | 2.813 GeV | **3.00** | 0.086 |
| Giả thiết cũ | 1.000 | 200.00 | $6.94\times10^{-56}$ | $2.50\times10^{-5}$ | 200.00 | 20.00 GeV | 21.32 | 0.033 |
| $A=16$ | 16.271 | 200.00 | $4.27\times10^{-57}$ | $6.62\times10^{-3}$ | 12.29 | 1.229 GeV | 1.31 | 0.132 |
| $A=40$ | 11.988 | 200.00 | $5.79\times10^{-57}$ | $3.59\times10^{-3}$ | 16.68 | 1.668 GeV | 1.78 | 0.113 |
| $A=208$ | 6.920 | 200.00 | $1.00\times10^{-56}$ | $1.20\times10^{-3}$ | 28.90 | 2.890 GeV | 3.08 | 0.086 |

**Hệ quả:** ở thang QCD, liên kết đồng bộ **thực tế không thể gãy** bởi các quá trình năng lượng hạt nhân.
- Muốn gãy, động cơ $\mathcal E_{\text{kép}}$ phải quay với năng lượng lớn hơn năng lượng nghỉ của proton.
- Ở mức năng lượng đó, mô hình Lagrange cổ điển phi tương đối tính không còn áp dụng được.
- Hạt thì vẫn chậm ($v/c\le0.13$). Phần vượt giới hạn là **động cơ**, không phải hạt.

## 2. Cộng hưởng $I=\chi$ so với $I=1$

**Đơn vị không thứ nguyên:** $m=\omega_0=\hbar=\ell_0=1$.

**Thông số:** $\Delta E=\chi$, $E_{\text{vt}}=\varepsilon$, $\kappa=\chi/\varepsilon$, $R=\varepsilon$.

| $g$ | $\chi$ | $\varepsilon/\sqrt\chi$ | $I$ | Góc giữ tĩnh | $\dot\theta_{\text{gãy}}$ | Tiệm cận | Sai lệch |
|---|---|---|---|---|---|---|---|
| 0 | 0.1 | 0.5 | $=\chi$ | 62.05° | 0.91611 | 0.5 | −45.42% |
| 0 | 0.1 | 0.5 | $=1$ | 21.75° | 0.22256 | 0.5 | +124.66% |
| 0 | 0.1 | 1 | $=\chi$ | 89.99° | **1.22539** | 1 | −18.39% |
| 0 | 0.1 | 1 | $=1$ | 89.99° | **1.43613** | 1 | −30.37% |
| 0 | 0.1 | 5 | $=\chi$ | 180° | 5.17817 | 5 | −3.44% |
| 0 | 0.1 | 5 | $=1$ | 180° | 5.10860 | 5 | −2.13% |
| 0 | 1 | 0.5 / 1 / 5 | $=\chi=1$ | 62.05° / 90.00° / 180° | 0.91611 / 1.22539 / 5.17817 | | −45.42% / −18.39% / −3.44% |
| 0 | 10 | 0.5 | $=\chi$ | 62.05° | 0.91611 | 0.5 | −45.42% |
| 0 | 10 | 0.5 | $=1$ | 51.80° | 0.88858 | 0.5 | −43.73% |
| 0 | 10 | 1 | $=\chi$ | 90.00° | 1.22539 | 1 | −18.39% |
| 0 | 10 | 1 | $=1$ | 60.68° | 1.25264 | 1 | −20.17% |
| 0 | 10 | 5 | $=\chi$ | 180° | 5.17817 | 5 | −3.44% |
| 0 | 10 | 5 | $=1$ | 180° | 5.34633 | 5 | −6.48% |
| 0.1 | 0.1 | 0.5 | $=\chi$ / $=1$ | 62.76° / 24.47° | 0.80633 / 1.04100 | 0.5 | −37.99% / −51.97% |
| 0.1 | 1 | 1 | $=\chi$ | 96.92° | 1.23250 | 1 | −18.86% |
| 0.1 | 10 | 1 | $=\chi$ / $=1$ | 91.90° / 61.45° | 1.22633 / 1.24545 | 1 | −18.46% / −19.71% |

(Toàn bộ 30 dòng, gồm cả $g=0.1$, có trong output của script.)

**Đọc bảng:**

1. **Bất biến theo $\chi$:** khi $I=\chi$ và $g=0$, các dòng $\chi=0.1,1,10$ cho **cùng** góc giữ và cùng $\dot\theta_{\text{gãy}}$. Lý do: cộng hưởng làm $L$ tỉ lệ đều theo $\chi$, nên chỉ còn một tham số $\varepsilon/\sqrt\chi$.
2. **Cộng hưởng mở rộng góc giữ tĩnh** gần biên, nhưng không đưa góc giữ lên 180°. Ở $\varepsilon/\sqrt\chi=0.5$, liên kết vẫn gãy tĩnh ở 62°, do các số hạng phi tuyến bậc cao trong $\lambda$.
3. **Cộng hưởng không luôn làm liên kết bền hơn.** Ở $\chi=0.1$, $\varepsilon/\sqrt\chi=1$: $\dot\theta_{\text{gãy}}$ với $I=\chi$ là 1.225, **thấp hơn** 1.436 của $I=1$.
4. **Với $\varepsilon/\sqrt\chi\gtrsim5$**, mọi lựa chọn $I$ đều giữ tĩnh 180° và $\dot\theta_{\text{gãy}}$ chỉ lệch vài %. Ở thang QCD ($\approx197$), lựa chọn $I$ không còn vai trò.

## 3. Cửa sổ $T$ và dải $e^*$

### 3.1. Đọc lại bảng quét cũ với đúng $T_0=\pi$ ($k=4$)

Thời điểm khóa **đầu tiên**, lấy từ [chi-va-dai-hoi-phuc-kiem-chung.md](chi-va-dai-hoi-phuc-kiem-chung.md):

| $e$ | $t_{\text{khóa}}$ ($\dot\theta_0=5$) | $t_{\text{khóa}}$ ($\dot\theta_0=6$) | Trong $3T_0=9.42$? | Trong $10T_0=31.42$? |
|---|---|---|---|---|
| 0.0 | 1.41 | 1.41 | có / có | có / có |
| 0.1 | 13.87 | 2.31 | **không** / có | có / có |
| 0.2 | 24.33 | 30.03 | không / không | có / có |
| 0.3 | 46.62 | 54.59 | không / không | **không / không** |
| ≥0.4 | >200 (trừ $e=0.9$, $\dot\theta_0=6$: 175.65) | | không | không |

Vậy $e^*(10T_0)\approx0.2$, còn $e^*(3T_0)\approx0$–$0.1$. **Không phải 0.3.**

### 3.2. Mô phỏng tại đúng bộ thông số đề xuất (QCD + cộng hưởng)

**Thông số:** $\chi=1$, $I=1$, $\Delta E=1$, $k=m=1$ ($T_0=2\pi$), $\varepsilon=197.327$, $\kappa=0.005068$, $R=197.327$, $g=0$.

**Điều kiện đầu:** $\dot\theta_0=1.05\,\dot\theta_{\text{gãy}}$. Cửa sổ $T=10T_0=62.83$.

**Ngưỡng:** $\dot\theta_{\text{gãy}}$ chính xác là 197.3320 (tiệm cận 197.3270). Năng lượng tại ngưỡng là 19 469.469. Năng lượng ban đầu $E_0=21\,465.192$, nên phần năng lượng dư phải tiêu tán là 1 995.723.

| $e$ | Gãy | Va chạm (khóa) | $E$ mất trong $10T_0$ | % phần dư | Mất TB/va chạm | Va chạm cần để về đồng bộ | Cuối |
|---|---|---|---|---|---|---|---|
| 0.0 | 21 | 20 (20) | 12.793 | 0.64% | 0.640 | ≈ 3 120 | **Tự do** |
| 0.3 | 1 | 116 (0) | 16.396 | 0.82% | 0.141 | ≈ 14 120 | **Tự do** |
| 0.6 | 1 | 37 (0) | 6.915 | 0.35% | 0.187 | ≈ 10 678 | **Tự do** |

**Cơ chế:**
- Với $e=0$, mỗi lần khóa lại, $\dot\theta$ vẫn ở trên ngưỡng, nên liên kết **gãy lại ngay lập tức**.
- Mỗi va chạm chỉ lấy đi cỡ $\tfrac12m_{\text{eff}}u^2\sim\tfrac12m(\kappa\dot\theta)^2\approx0.6$. Con số này rất nhỏ so với phần dư ~2000 của động cơ.
- Mỗi chu kỳ gãy/khóa mất khoảng nửa chu kỳ dao động ($\approx\pi$). Do đó thời gian về đồng bộ vào cỡ $3120\pi\approx9\,800$, tức **≈ 1 560 $T_0$**, vượt xa cửa sổ $10T_0$.

**Tổng quát:** tỉ số giữa năng lượng dư của động cơ và năng lượng mất mỗi va chạm tăng theo $(\varepsilon/\sqrt\chi)^2$. Ở thang QCD, cơ chế tự ổn định (Carnot) quá chậm so với mọi cửa sổ vài $T_0$.

## 4. Kết luận chuẩn hóa đã sửa

1. $I=\chi\hbar/\omega_0$ là một cách chuẩn hóa hợp lệ. Nó làm động lực học bất biến theo $\chi$, nhưng chỉ có tác dụng khi $E_{\text{vt}}\lesssim5\sqrt\chi\,\hbar\omega_0$. Ở thang QCD, nó không thay đổi kết quả.
2. Với $L\approx1$ fm ($E_{\text{vt}}\approx200$ MeV), liên kết đồng bộ **rất bền**. Muốn gãy cần năng lượng quay $E_{\text{vt}}^2/(2\hbar\omega_0)$, tức 1.2–20 GeV, lớn hơn $m_pc^2$. "E gãy" vì vậy **không xảy ra** trong vật lý hạt nhân thông thường.
3. Nếu cưỡng bức làm gãy (ví dụ vượt ngưỡng 5%), hệ **không** tái đồng bộ trong $10T_0$ với bất kỳ $e$ nào. Việc tái đồng bộ cần cỡ $10^3\,T_0$.
4. Dải "$0\le e\le0.3$ trong $10T_0$" không đứng được với bộ thông số nào đã thử:
   - Bộ cũ ($k=4$): $e^*(10T_0)\approx0.2$.
   - Bộ QCD: không có $e$ nào khóa được trong cửa sổ.

## Câu hỏi còn mở

- **Mâu thuẫn cốt lõi** giữa hai yêu cầu:
  - (a) Liên kết phải tồn tại, đòi $E_{\text{vt}}\gg\sqrt\chi\hbar\omega_0$.
  - (b) "E gãy" và tái đồng bộ phải xảy ra ở năng lượng hạt nhân và trong vài $T_0$, đòi $E_{\text{vt}}$ chỉ cỡ vài lần $\sqrt\chi\hbar\omega_0$.

  Thang QCD thỏa (a) nhưng phá (b). Học thuyết cần chọn: có muốn $E_{\text{vt}}/(\sqrt\chi\hbar\omega_0)\sim2$–$5$ hay không? Nếu có, $L$ phải cỡ 40–100 fm (với $\chi=1$, $\hbar\omega_0=1$ MeV), không phải 1 fm.
- **Cửa sổ decoherence** là $T_{\text{QCD}}$ ($\approx8\times10^{-4}\,T_0$) hay vài $T_0$? Với $T_{\text{QCD}}$ thì không có tái đồng bộ.
- **$\hbar\omega_0$** dùng giá trị mô hình vỏ (7–16 MeV) hay 1 MeV? Giá trị này làm năng lượng gãy thay đổi 7–16 lần.
