# Bộ thông số hạt nhân ($\hbar\omega_0=10$ MeV, $\chi=1$, $E_{\text{vt}}=30$ MeV): giải và kiểm chứng số liệu

> Tiếp nối [cong-huong-qcd-kiem-chung.md](cong-huong-qcd-kiem-chung.md).
> Số liệu tính bằng [giai_so_bo_thong_so_hat_nhan.py](giai_so_bo_thong_so_hat_nhan.py) (chạy `python3 giai_so_bo_thong_so_hat_nhan.py`, ≈ 5.5 phút).
> Hằng số CODATA: $\hbar c=197.327$ MeV·fm, $\hbar=6.582\times10^{-22}$ MeV·s, $m_pc^2=938.272$ MeV.

## Tóm tắt kết luận

**Về động lực học, bộ thông số này đạt yêu cầu:**
- Giữ tĩnh 180°.
- Gãy ở năng lượng quay 53 MeV, tức phi tương đối tính.
- Tái đồng bộ trong $10T_0$ ở **12/12** lần chạy với $e\le0.2$.

**Về cửa sổ thời gian, lập luận không đứng được:** $\tau_{\text{relax}}$ mà chính đề xuất tính ra chỉ bằng $0.20\,T_0$, không phải $10$–$30\,T_0$. Trong cửa sổ đó, **0/12** lần chạy khóa lại được.

Ngoài ra, 4 con số trong bảng đề xuất bị tính sai ($\lambda$, $\tau/T_0$, $I$, $R$).

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| Vùng $E_{\text{vt}}/(\sqrt\chi\hbar\omega_0)\approx2$–$5$ | **Hợp lý**, và là vùng mô phỏng xác nhận hoạt động được |
| Tiêu đề: "$L\sim20$–$80$ fm" | **Không khớp** với chính bộ số chọn. Với $\hbar\omega_0=10$ MeV, tỉ số 2–5 cho $E_{\text{vt}}=20$–$50$ MeV, tức $L=$ **3.95–9.87 fm**. Giá trị $L=6.58$ fm được chọn nằm đúng trong khoảng này |
| Bước sóng de Broglie nucleon 5 MeV là $\approx20$–$40$ fm | **Sai.** $pc=96.86$ MeV, nên $ƛ=\hbar c/pc=$ **2.04 fm** và $\lambda=2\piƛ=$ **12.8 fm** |
| $\tau_{\text{relax}}=5\text{ fm}/0.2c\approx10^{-22}$ s | **Đúng số:** $8.34\times10^{-23}$ s |
| $\tau_{\text{relax}}\approx(10$–$30)\,T_0$, hoặc "$\sim5T_0$" trong bảng | **Sai.** Với $\hbar\omega_0=10$ MeV, $T_0=4.14\times10^{-22}$ s, nên $\tau_{\text{relax}}=$ **0.20 $T_0$**. (Với $\hbar\omega_0=1$ MeV thì chỉ còn $0.02\,T_0$.) Hơn nữa $v_F$ thực là $0.275c$ chứ không phải $0.2c$, làm $\tau$ còn ngắn hơn |
| Cửa sổ $T\sim(5$–$10)T_0$ "có căn cứ" | **Không có căn cứ** từ $\tau_{\text{relax}}$. $10T_0=4.1\times10^{-21}$ s dài gấp **50 lần** $\tau_{\text{relax}}$ |
| $I=\chi\hbar/\omega_0=6.94\times10^{-58}$ kg·m² | **Sai 10 lần:** đúng là **$6.94\times10^{-57}$** kg·m² |
| $\kappa=\ell_0E_{\text{kép}}/E_{\text{vt}}\approx0.68$ fm | **Đúng:** $\ell_0=2.037$ fm, $\kappa=0.679$ fm |
| $R=E_{\text{vt}}/\ell_0\approx1.5\times10^4$ N | **Sai ~6.4 lần:** $R=14.73$ MeV/fm = **2 359 N** (1 MeV/fm = 160.2 N) |
| $k\kappa/R=0.11$, giữ tĩnh 180° | **Đúng:** 0.1111, góc giữ 180.00° |
| $\dot\theta_{\text{gãy}}/\omega_0\approx3.2$ | **Đúng:** giá trị chính xác là 3.255 (tiệm cận 3.0) |
| $E_{\text{quay}}=E_{\text{vt}}^2/(2\hbar\omega_0)=45$ MeV $\approx0.048\,m_pc^2$ | **Đúng tiệm cận.** Giá trị chính xác là **53.0 MeV** $=0.056\,m_pc^2$. Vẫn phi tương đối tính, $v/c=0.103$ |
| "Cần 4–8 va chạm, $t\approx12$–$25$" | **Gần đúng với $e=0$:** 2–9 va chạm, tái đồng bộ ở $t=5.1$–$19.0$. Với $e=0.1$–$0.2$ cần 6–32 va chạm, $t$ lên tới 53.2 |
| Tái khóa "trong đúng khoảng thời gian sống vi mô" | **Sai.** Với cửa sổ $\tau_{\text{relax}}=1.27$ (đơn vị $1/\omega_0$), không lần chạy nào khóa kịp. Lần sớm nhất là $t=3.72$, gấp gần 3 lần $\tau_{\text{relax}}$ |

## 1. Các đại lượng SI

| Đại lượng | Giá trị tính | Đề xuất |
|---|---|---|
| $\omega_0$ | $1.519\times10^{22}$ rad/s | — |
| $T_0=2\pi\hbar/\hbar\omega_0$ | $4.136\times10^{-22}$ s | — |
| $\tau_{\text{relax}}=5\text{ fm}/0.2c$ | $8.339\times10^{-23}$ s $=0.202\,T_0$ | "$10$–$30\,T_0$", "$\sim5T_0$" |
| $v_F$ thực ($k_F=1.36$ fm⁻¹) | $0.275c$ | $0.2c$ |
| $I=\chi\hbar/\omega_0$ | $6.941\times10^{-57}$ kg·m² | $6.94\times10^{-58}$ |
| $\ell_0=\hbar c/\sqrt{mc^2\hbar\omega_0}$ | 2.0371 fm | — |
| $\kappa$ | 0.6790 fm | 0.68 fm ✓ |
| $R=E_{\text{vt}}/\ell_0$ | 14.726 MeV/fm = 2 359.4 N | $1.5\times10^4$ N |
| $L=\hbar c/E_{\text{vt}}$ | 6.578 fm | 6.5 fm ✓ |
| $k\kappa/R$ | 0.1111 | 0.11 ✓ |
| $v/c$ lúc gãy | 0.1032 | 0.1 ✓ |
| de Broglie 5 MeV | $ƛ=2.037$ fm, $\lambda=12.80$ fm | 20–40 fm |

## 2. Động lực học đầy đủ

**Đơn vị:** $m=\omega_0=\hbar=\ell_0=1$. Một đơn vị năng lượng bằng $\hbar\omega_0=10$ MeV. $T_0=2\pi$.

**Thông số:** $I=\chi=1$, $\Delta E=1$, $\kappa=1/3$, $R=3$, $g=0$.

**Kết quả tĩnh:**
- Góc giữ tĩnh: **180.00°**.
- $\dot\theta_{\text{gãy}}$ chính xác: **3.25544** (tiệm cận 3.0).
- Năng lượng quay tới hạn: $\tfrac12I\dot\theta_{\text{gãy}}^2=5.2989$ đơn vị, tức **52.99 MeV**.

**Các cửa sổ thời gian** (đơn vị $1/\omega_0$):

| Cửa sổ | Giá trị |
|---|---|
| $\tau_{\text{relax}}$ | 1.267 |
| $5T_0$ | 31.416 |
| $10T_0$ | 62.832 |

### Chu kỳ gãy → va chạm → tái đồng bộ

Điều kiện đầu: $\dot\theta_0=f\cdot\dot\theta_{\text{gãy}}$. "E dư" là năng lượng vượt ngưỡng. "$t$ đồng bộ cuối" là lần khóa cuối cùng mà sau đó không gãy lại.

| $f$ | E dư (MeV) | $e$ | Gãy | Va chạm | $t$ khóa đầu | $t$ đồng bộ cuối | Trong $\tau_{\text{relax}}$? | Trong $5T_0$? | Trong $10T_0$? | E mất / E dư |
|---|---|---|---|---|---|---|---|---|---|---|
| 1.05 | 6.04 | 0.0 | 2 | 2 | 2.573 | 5.065 | không | có | có | 111% |
| 1.05 | 6.04 | 0.1 | 1 | 6 | 3.725 | 3.725 | không | có | có | 106% |
| 1.05 | 6.04 | 0.2 | 1 | 22 | 36.299 | 36.299 | không | **không** | có | 845% |
| 1.20 | 25.91 | 0.0 | 3 | 3 | 3.264 | 7.501 | không | có | có | 130% |
| 1.20 | 25.91 | 0.1 | 2 | 20 | 6.813 | 38.706 | không | **không** | có | 251% |
| 1.20 | 25.91 | 0.2 | 1 | 8 | 7.378 | 7.378 | không | có | có | 103% |
| 1.50 | 73.60 | 0.0 | 5 | 5 | 2.924 | 13.761 | không | có | có | 102% |
| 1.50 | 73.60 | 0.1 | 1 | 14 | 25.646 | 25.646 | không | có | có | 153% |
| 1.50 | 73.60 | 0.2 | 1 | 12 | 20.618 | 20.618 | không | có | có | 147% |
| 2.00 | 176.63 | 0.0 | 9 | 9 | 3.051 | 19.007 | không | có | có | 103% |
| 2.00 | 176.63 | 0.1 | 2 | 19 | 7.141 | 20.304 | không | có | có | 101% |
| 2.00 | 176.63 | 0.2 | 3 | 32 | 6.353 | 53.227 | không | **không** | có | 124% |

**Đọc kết quả:**

1. **Trong $10T_0$:** tái đồng bộ **12/12**. **Trong $5T_0$:** 9/12. **Trong $\tau_{\text{relax}}$:** **0/12**.
2. **Tỉ số E mất / E dư luôn ≥ 100%.** Nghĩa là hệ không chỉ tiêu tán phần dư mà còn tụt xuống dưới ngưỡng rồi mới khóa. Trường hợp $e=0.2$, $f=1.05$ mất tới 8.5 lần phần dư, vì phải nảy 22 lần trước khi khóa.
3. **$e=0$ ổn định nhất:** 2–9 va chạm, đồng bộ ở $t=5.1$–$19.0$ ($0.8$–$3.0\,T_0$). Gần với con số "4–8 va chạm, $t\approx12$–$25$" của đề xuất.
4. **So với thang QCD** ([cong-huong-qcd-kiem-chung.md](cong-huong-qcd-kiem-chung.md)), bộ số này giải được mâu thuẫn năng lượng: gãy ở 53 MeV thay vì 1.2–20 GeV, và va chạm cần là 2–32 thay vì khoảng 3 000.

## 3. Mâu thuẫn còn lại: thang thời gian

Mâu thuẫn đã chuyển từ **năng lượng** sang **thời gian**. Động lực học cần $\approx1$–$8\,T_0$ để tái đồng bộ, trong khi thời gian hồi phục hạt nhân ước tính chỉ $\approx0.2\,T_0$. Bộ thông số chỉ khép kín được nếu thời gian giữ kết hợp pha dài hơn $\tau_{\text{relax}}$ khoảng **25–50 lần**, tức $\gtrsim2$–$4\times10^{-21}$ s.

Hai cách để mô hình nhất quán:
- **(a)** Chứng minh cơ chế decoherence của cặp $\mathcal E_{\text{kép}}$–hạt yếu hơn va chạm nucleon–nucleon. Cần một đại lượng độc lập, không lấy $\lambda_{\text{mfp}}/v_F$.
- **(b)** Chấp nhận rằng trong hạt nhân thực, tái đồng bộ **không** hoàn tất trước khi mất pha. Khi đó chu trình Quá khứ → Hiện tại → Tương lai chỉ xảy ra một phần.

## Câu hỏi còn mở

- Thời gian giữ kết hợp pha của $\mathcal E_{\text{kép}}$ được xác định từ đâu, nếu không phải $\tau_{\text{relax}}$ nucleon–nucleon?
- $g$ (hằng số ghép cặp của hệ (1)–(2)) đang đặt bằng 0 trong mô phỏng. Giá trị thực của nó trong bộ thông số hạt nhân là bao nhiêu?
