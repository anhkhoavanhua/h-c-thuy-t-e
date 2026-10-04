# Bộ đếm vòng (số cuộn, năng lượng tiêu tán) và mở rộng Lindblad: giải và kiểm chứng số liệu

> Tiếp nối [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md).
> Số liệu tính bằng [giai_so_winding_lindblad.py](giai_so_winding_lindblad.py) (cần numpy). Chạy `python3 giai_so_winding_lindblad.py`: phần B mất vài giây, phần A ≈ 2 phút trên 8 lõi. Script **không sửa** các file mô phỏng cũ.
> **Phần A (Langevin cổ điển):** bộ thông số hạt nhân, $f=1.05$, $\tau=10$, $\gamma_x=0.7$, $g=0.1$; 8 quỹ đạo × $60T_0$.
> **Phần B (Lindblad):** $\hbar=\omega_0=1$, $\gamma=\gamma_x=0.7$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| **Cách 1.** Nâng $\theta\in S^1$ lên $\Phi\in\mathbb R$; số cuộn $N=\lfloor\Phi/2\pi\rfloor$ phân định quá khứ → hiện tại → tương lai | **Không đơn điệu, nên không phải mũi tên thời gian.** Với nhiễu 1 MeV: **8.4%** số bước đếm là **lùi**, 13 lần đảo chiều quay, **1/8 quỹ đạo quay ngược ròng**. Với 3 MeV: 38% bước lùi, 3/8 quỹ đạo lùi ròng. Với 5 MeV: **53%** bước lùi, và tốc độ quay trung bình **âm** (−0.15 vòng/$T_0$). $N$ là bộ đếm **có hướng theo chiều quay**, không phải theo chiều thời gian |
| $N(t)$ làm đồng hồ | **Kém đều:** CV gộp của $\Delta N$ mỗi $T_0$ là **0.72** (1 MeV), 4.2 (3 MeV), 15.9 (5 MeV) |
| Nâng lên $\mathbb R$ về mặt lượng tử | **Đây là một lý thuyết khác, không chỉ là cách ghi.** Trên $S^1$, $p_\theta=\hbar\ell$ với $\ell$ nguyên (phổ rô-to rời rạc, đúng như đã dùng ở S1/PW). Trên $\mathbb R$ với thế tuần hoàn $-\Delta E\cos\theta$, trạng thái được gắn nhãn theo **chuẩn xung lượng Bloch** và cho phổ dải. Cần chọn một trong hai và nói rõ |
| **Cách 2.** $E_{\text{diss}}(t)$ tăng đơn điệu, làm "con trỏ entropy" | **Đơn điệu ✓** theo định nghĩa. **Nhưng là đồng hồ rất thô:** riêng phần Carnot tăng theo từng cục, **49–68%** số cửa sổ dài $T_0$ **không tăng chút nào**, CV gộp 1.5–1.8. Phần ma sát $\gamma_x\dot x^2$ làm $E_{\text{diss}}$ tăng **nghiêm ngặt**, nhưng với tốc độ thay đổi theo pha. Độ phân giải thời gian thực tế chỉ vài $T_0$ |
| Lindblad "tạo mũi tên thời gian: mất kết hợp diễn ra theo một chiều" | **Cần phát biểu chính xác hơn.** Entropy von Neumann $S(\rho)$ **không đơn điệu**: xuất phát từ $\lvert2\rangle$, $S$ tăng tới **1.04** (1 MeV) hoặc **1.23** (5 MeV), rồi **giảm** về 0.005 hoặc 0.459. Đại lượng **đơn điệu** là **entropy tương đối** $D(\rho\Vert\rho_{\text{ss}})$ tới trạng thái dừng (định lý Spohn): 19.9 → $9\times10^{-4}$ và 4.07 → $6\times10^{-7}$, **giảm đơn điệu** ✓ |
| Toán tử nhảy $L_j=\sqrt{\Gamma_j}\,a^j$ | **$j\ge2$ gần như không đóng góp.** Với năng lượng hạt cổ điển 0.571 $\hbar\omega_0$ (tức $\langle n\rangle\approx0.071$): tốc độ $j=2$ chỉ bằng **7.1%** ($j=1$) nếu trạng thái kết hợp, **14.2%** nếu trạng thái nhiệt; $j=3$ còn nhỏ hơn 20–30 lần. Chỉ cần **$L\propto a$** (mỗi bước nhảy đúng $1\hbar\omega_0$) |
| "Tái lập hoàn hảo" công suất tiêu tán cổ điển | **Chỉ đúng gần cộng hưởng.** Hạt bị lái bởi $-g\,x\sin\Omega t$: tại $\Omega=1$ (cộng hưởng), Lindblad/cổ điển = **1.036**. Tại $\Omega=\dot\theta_d=2.63$ (tần số của mô hình), tỉ số = **0.609** (lệch −39%). Nguyên nhân: tiêu tán Lindblad dạng $a$ (xấp xỉ sóng quay) khác ma sát vận tốc khi lái xa cộng hưởng. Cần bộ tiêu tán không-RWA (kiểu Caldeira–Leggett), hoặc hiệu chỉnh $\Gamma$ theo tần số |
| $\langle\Delta E\rangle=\sum j\hbar\omega_0P(j)\approx0.40\,\hbar\omega_0$ mỗi va chạm | **Nhất quán về cân bằng năng lượng (cỡ bậc).** Tốc độ khóa ≈ 0.152 đợt/$T_0$ × 5.12 ≈ 0.78 lần/$T_0$, nhân với 0.40 cho ≈ **0.050 $\hbar\omega_0$ mỗi đơn vị thời gian**. Tốc độ nhảy Lindblad $\gamma\langle n\rangle=0.7\times0.071\approx$ **0.050**. Hai con số cùng cỡ, nhưng $\langle n\rangle$ ở đây là ước lượng thô từ năng lượng đồng bộ cực đại. **Lưu ý vật lý:** trong mô hình cổ điển, năng lượng Carnot phần lớn đến từ **rô-to**, nên toán tử nhảy nên tác dụng trên **không gian chung rô-to–hạt** (ví dụ chứa $e^{-i\theta}$), không chỉ trên hạt |

## 1. Bộ đếm vòng trong Langevin cổ điển

| $k_BT$ | Tốc độ TB (vòng/$T_0$) | Bước $N$ lùi | Số lần đảo chiều quay | Lùi xa nhất TB | Quỹ đạo lùi ròng | CV gộp $\Delta N$/$T_0$ | CV gộp $\Delta E_{\text{Carnot}}$/$T_0$ | Cửa sổ $T_0$ không có Carnot |
|---|---|---|---|---|---|---|---|---|
| 1 MeV | 2.042 | 99/1184 (**8.4%**) | 13 | 12.0 vòng | **1/8** | 0.72 | 1.83 | **67.7%** |
| 3 MeV | 0.556 | 431/1124 (**38.3%**) | 151 | 47.3 vòng | 3/8 | 4.16 | 1.64 | 59.8% |
| 5 MeV | **−0.152** | 604/1135 (**53.2%**) | 163 | 63.0 vòng | **5/8** | 15.9 | 1.53 | 49.0% |

**Đọc bảng:**
- Số cuộn đo chiều quay của $\mathcal E_{\text{kép}}$. Nhiễu làm rô-to đảo chiều, nên $N$ lùi.
- Nếu học thuyết đồng nhất "chiều thời gian" với chiều quay, thì ở nhiễu ≥ 3 MeV thời gian sẽ chạy ngược trong nhiều quỹ đạo. Điều đó mâu thuẫn với chính khái niệm mũi tên thời gian.

## 2. Lindblad

### B1. Mũi tên thời gian: quá trình từ $\lvert2\rangle$ về trạng thái dừng

| $k_BT$ | $\bar n$ | $S(\rho)$ max (thời điểm) | $S$ cuối / $S_{\text{ss}}$ | $D(\rho\Vert\rho_{\text{ss}})$ đầu → cuối | $D$ đơn điệu? |
|---|---|---|---|---|---|
| 1 MeV | $4.5\times10^{-5}$ | **1.040** ($t=1.01$) | 0.0048 / 0.0005 | 19.9 → $9.0\times10^{-4}$ | **có** |
| 5 MeV | 0.157 | **1.228** ($t=0.81$) | 0.459 / 0.458 | 4.07 → $6.2\times10^{-7}$ | **có** |

### B2. Công suất tiêu tán: Lindblad ($L=\sqrt\gamma a$) so với cổ điển ($\gamma_x\dot x^2$)

| Tần số lái | Lindblad $\gamma\langle n\rangle$ | Cổ điển $\gamma_x\langle\dot x^2\rangle$ | Tỉ số |
|---|---|---|---|
| $\Omega=1$ (cộng hưởng) | $7.40\times10^{-3}$ | $7.14\times10^{-3}$ | 1.036 |
| $\Omega=2.63$ (tốc độ quay dừng của mô hình) | $3.84\times10^{-4}$ | $6.31\times10^{-4}$ | **0.609** |

### B3. Bước nhảy nhiều lượng tử $a^j$ ($\langle n\rangle=0.071$)

| Trạng thái hạt | $\langle a^\dagger a\rangle$ | $\langle a^{\dagger2}a^2\rangle$ | $\langle a^{\dagger3}a^3\rangle$ | $j{=}2$ : $j{=}1$ |
|---|---|---|---|---|
| Kết hợp | 0.0710 | 0.00504 | 0.000358 | 0.071 |
| Nhiệt | 0.0710 | 0.01008 | 0.002147 | 0.142 |

## 3. Đề xuất sửa đoạn văn cho bài báo

**Về chiều thời gian:** câu "chuỗi ba pha được định hướng nhờ tính đơn điệu của năng lượng tiêu tán tích lũy" là **đúng**, và nên dùng thay cho số cuộn. Cần thêm hai ý:
- $E_{\text{diss}}$ tăng đơn điệu nhưng không đều, nên độ phân giải thời gian chỉ cỡ vài $T_0$.
- Số cuộn $N$ **không** đơn điệu khi có nhiễu, nên **không** dùng làm mũi tên.

**Đoạn Lindblad (tiếng Anh), sửa ba chỗ:**

> "Within the classical model, the nonlinear breakdown threshold and the avalanche statistics are obtained numerically. The sub-quantum Carnot dissipation ($\approx0.4\,\hbar\omega_0$ per impact) indicates, however, that the system operates in the strong quantum regime. The natural quantum extension is an open-system (Lindblad) description in which plastic impacts are replaced by **single-quantum jump operators acting on the joint rotor–particle space**, with a **non-rotating-wave dissipator** so that off-resonant driving at $\Omega\approx2.6\,\omega_0$ reproduces the classical dissipated power. In this description irreversibility appears as the **monotonic decrease of the relative entropy to the stationary state**, providing an arrow of time that the unitary Page–Wootters and Wheeler–DeWitt formulations lack."

Ba chỗ đã sửa so với bản đề xuất:
1. Bỏ "accurately captures", vì chưa có số liệu thực nghiệm để so.
2. Thay $a^j$ bằng bước nhảy đơn lượng tử trên không gian chung, kèm bộ tiêu tán không-RWA.
3. Nói rõ mũi tên thời gian là **entropy tương đối giảm đơn điệu**, vì entropy von Neumann **không** đơn điệu.

## Câu hỏi còn mở

- **Chiều quay của $\mathcal E_{\text{kép}}$ có ý nghĩa vật lý không?** Nếu có, việc rô-to đảo chiều dưới nhiễu (53% bước lùi ở 5 MeV) cần được giải thích hoặc bị loại, ví dụ bằng cách giới hạn học thuyết ở $k_BT\lesssim1$ MeV.
- **Bài Part II (Lindblad):** toán tử nhảy nào trên không gian chung (rô-to–hạt) tương ứng với va chạm dẻo cổ điển? Một ứng viên là $L\propto e^{-i\theta}a$ (hạ đồng thời mô-men rô-to và lượng tử hạt). Cần kiểm xem nó có tái lập được bội số hình học và thời gian lưu dạng mũ không.
