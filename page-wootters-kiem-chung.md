# Page–Wootters cho hệ $(\theta,x)$: chạy script đề xuất và kiểm chứng

> Số liệu tính bằng [giai_so_page_wootters.py](giai_so_page_wootters.py) (cần numpy; chạy vài giây).
> Phần (A) là script đề xuất, **giữ nguyên**. Chỉ thay `scipy.linalg.eigh` bằng `numpy.linalg.eigh` (máy không có scipy; hai hàm cho cùng kết quả với ma trận Hermite) và bỏ matplotlib (script không vẽ gì).
> **Đơn vị:** $\hbar=m=\omega_0=1$, $I=1$, $\Delta E=1$, $g=0.1$. **Cơ sở tính:** $\ell\in[-8,8]$, $n\le6$ (tổng 119 chiều).

## Tóm tắt

| Khẳng định / dự báo | Kết quả |
|---|---|
| Tìm $\lvert\Psi\rangle\rangle$ với $\hat{\mathbb H}\lvert\Psi\rangle\rangle=0$ bằng "trị riêng gần 0 nhất" | **Chọn nhầm trạng thái cơ bản.** Trị riêng gần 0 nhất là $E=-0.0359$, chính là **mức thấp nhất** của phổ. $\hat{\mathbb H}$ không có trị riêng bằng 0, và "0" chỉ là một cách chọn gốc năng lượng. Ở trạng thái cơ bản, con lắc nằm ở đáy: $p(\theta)$ chênh **340 lần** giữa giá trị lớn nhất và nhỏ nhất. **Đồng hồ không quay** |
| Kết quả script: $\langle x\rangle\in[-0.0628,+0.0628]$ | **Biên độ chỉ bằng 19% của κ** (không phải ±0.333). $\langle H_S\rangle=0.501$, tức hạt ở trạng thái cơ bản |
| "Hạt dao động với pha bám theo $\sin\theta_0$, mô phỏng lại liên kết đồng bộ cổ điển $x\approx\kappa\sin\theta_0$" | **Không đúng.** $\hat{\mathbb H}$ **không chứa** liên kết $x=\kappa\sin\theta$, chỉ có số hạng ghép $-g\hat x\sin\theta$. Kết quả ứng với **phản ứng tuyến tính với lực ghép**, không phải liên kết. Với đồng hồ đứng yên, biên độ $\propto g/k=0.1$. Với đồng hồ **quay một chiều** ($\langle\ell\rangle=2.98$): $\langle x\rangle\approx$ **−0.023 sinθ**, tức **7% của κ và ngược pha**. Kết quả khớp dấu với nghiệm cưỡng bức cổ điển $g/(k-\Omega^2)=-0.013$ |
| Quét $\theta_0\in[0,4\pi]$ ("2 chu kỳ") | Hai chu kỳ trùng nhau **theo định nghĩa**, vì θ tuần hoàn. Không có thông tin mới (xem [luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md)) |
| Dấu trong script | **Hai lỗi dấu triệt tiêu nhau.** Ma trận `sin_theta` của script là **−sinθ** (sai lệch 0 khi so với −sinθ đúng). Phép chiếu `clock_state.conj() @ psi` cho $\psi(-\theta)$ thay vì $\psi(\theta)$. Đổi θ → −θ cũng đổi dấu sinθ, nên tương quan $\langle x\rangle(\theta)$ ra đúng. Bản sửa nên dùng $\sin\theta=(e^{i\theta}-e^{-i\theta})/2i$ và $\psi_S(\theta)=\sum_\ell e^{i\ell\theta}\Psi_\ell$ |
| Mọi nghiệm WdW đều cho đồng hồ chạy một chiều | **Không.** $\hat{\mathbb H}$ đối xứng đảo thời gian, nên **mọi nghiệm riêng không suy biến có $\langle\ell\rangle=0$** (đo được max $3\times10^{-10}$ trên 16 nghiệm trong $E\in(3,6)$). Đồng hồ quay cả hai chiều như nhau. Muốn có một chiều, phải chồng chập **cặp gần suy biến** $\pm\ell$ |
| "PW đơn nhất chỉ mô tả được pha đồng bộ tĩnh; gãy và Carnot không sinh ra từ $\hat{\mathbb H}\lvert\Psi\rangle\rangle=0$" | **✓ Đúng hướng**, khớp các kiểm chứng trước. Thêm một điều chỉnh: ngay cả **pha đồng bộ** cũng không được mô tả, vì liên kết không có trong $\hat{\mathbb H}$. Nếu thêm liên kết dạng giếng, cơ học lượng tử chỉ giữ được hạt với độ nhòe ≥ 1.57κ (S2 trong tài liệu lượng tử) |

## Số liệu chi tiết

| Nghiệm | $E$ | $\langle\ell\rangle$ | $\sqrt{\langle\ell^2\rangle}$ | $p(\theta)$ max/min | $A$ (hệ số sinθ) | So sánh |
|---|---|---|---|---|---|---|
| "Gần 0" (= cơ bản) | −0.036 | 0.000 | 0.647 | **340** | +0.052 | $g/k=0.100$ |
| Riêng đơn | 3.334 | 0.000 | 1.083 | 15 | +0.732 | sóng đứng hai chiều, gần cộng hưởng |
| Riêng đơn | 3.568 | 0.000 | 1.940 | 54 | −0.375 | sóng đứng |
| Riêng đơn | 3.582 | 0.000 | 1.946 | 677 | −0.056 | sóng đứng |
| **Cặp ±ℓ, quay một chiều** | 5.029 | **2.980** | 2.990 | **1.26** (đồng hồ đều) | **−0.023** | $g/(k-\Omega^2)=-0.013$; κ = 0.333 |

$A$ là hệ số khớp $\langle x\rangle(\theta)\approx A\sin\theta+B\cos\theta+C$, lấy trọng số theo $p(\theta)$. Ở mọi dòng $B\approx0$.

**Đọc bảng:**
- Chỉ trạng thái **cặp ±ℓ** là đồng hồ đúng nghĩa: phân bố góc gần đều, và có chiều quay.
- Với đồng hồ đó, hạt dao động **ngược pha** với biên độ nhỏ. Đây là đáp ứng cưỡng bức trên tần số riêng ($\Omega\approx3>\omega_0$), đúng như nghiệm (2) của tài liệu đầu tiên. Nó **không** phải liên kết đồng bộ.
- Các nghiệm riêng đơn lẻ là **sóng đứng**. Hệ số $A$ lớn ở mức $E=3.334$ phản ánh sự trộn gần cộng hưởng giữa rô-to và hạt, không phải chuyển động theo thời gian.

## Script đã sửa dấu (thay các dòng tương ứng trong đề xuất)

```python
up = np.diag(np.ones(dim_C - 1), -1)          # e^{iθ}: |ℓ⟩ → |ℓ+1⟩
sin_theta = (up - up.T) / 2j                  # sinθ đúng dấu
cos_theta = (up + up.T) / 2
...
# chọn nghiệm: KHÔNG lấy trị riêng gần 0; lấy cặp gần suy biến ±ℓ ở năng lượng rô-to quay
# và chồng chập để ⟨ℓ⟩ > 0 (xem hàm checks trong giai_so_page_wootters.py)
...
psi_S = np.exp(1j * l_vals * th) @ psi_mat    # ψ_S(θ) = Σ_ℓ e^{iℓθ} Ψ_ℓ
```

## Kết luận cho học thuyết

1. **Page–Wootters với rô-to làm đồng hồ cho ra "thời gian" tuần hoàn**, và chỉ khi chọn đúng trạng thái: cặp ±ℓ quay một chiều. Nghiệm "gần 0" của đề xuất là trạng thái mà đồng hồ đứng yên.
2. **Hạt không bị khóa theo κ sinθ,** vì liên kết không có trong Hamiltonian. Mối tương quan thu được là đáp ứng tuyến tính: biên độ 0.02–0.06, và ngược pha khi đồng hồ quay nhanh.
3. **Gãy và Carnot không xuất hiện** trong mô tả đơn nhất này, như đề xuất đã nhận định. Chúng cần kênh mở (xem các tài liệu Lindblad).

## Câu hỏi còn mở

- **Muốn PW mô tả pha đồng bộ,** cần đưa liên kết vào $\hat{\mathbb H}$, ví dụ giếng bond $-D e^{-(x-\kappa\sin\theta)^2/2w^2}$. Mình có thể chạy lại PW với giếng này để xem biên độ $\langle x\rangle(\theta)$ tiến gần κ đến đâu. Theo S2, độ nhòe tối thiểu là 1.57κ.
- **Chọn nghiệm WdW nào?** Phương trình $\hat{\mathbb H}\Psi=0$ không tự chọn năng lượng. Học thuyết cần một nguyên tắc chọn, ví dụ năng lượng ứng với $\dot\theta\approx\dot\theta_d$, hoặc một trạng thái chồng chập.
