# Đồng hồ photon Page–Wootters trên học thuyết $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** bản thảo [HocThuyetE (5).pdf](../../HocThuyetE%20(5).pdf) (Part I, 5 trang). Bộ tham số lấy từ Bảng I: $\tau=10$, $\gamma_x=0.70$, $g=0.1$, $f=1.05$, $k_BT=1$ MeV, $\dot\theta_{\text{br}}=3.205$.
> **Thí nghiệm tham chiếu:** E. Moreva, G. Brida, M. Gramegna, V. Giovannetti, L. Maccone, M. Genovese, *Time from quantum entanglement: an experimental illustration*, Phys. Rev. A **89**, 052122 (2014).
> Số liệu tính bằng [giai_so_dong_ho_photon.py](giai_so_dong_ho_photon.py) (cần numpy). Chạy `python3 giai_so_dong_ho_photon.py`: phần A mất vài giây, phần B ≈ 2.5 phút trên 8 lõi. Phần B dùng lại mô phỏng Langevin của Part I, **không sửa** file cũ nào.

## 1. Tái lập thí nghiệm Moreva (dạng lý tưởng)

**Thiết lập:**
- Photon đồng hồ $C$ và photon hệ $S$ mang phân cực H/V, ở trạng thái vướng víu singlet $\lvert\Psi\rangle\rangle=(\lvert HV\rangle-\lvert VH\rangle)/\sqrt2$.
- Bản làm trễ quay phân cực cả hai photon, tức $\hat{\mathbb H}=\omega(\sigma_y^C+\sigma_y^S)$.

| Kiểm tra | Kết quả |
|---|---|
| $\lVert\hat{\mathbb H}\lvert\Psi\rangle\rangle\rVert$ | **0** (nghiệm Wheeler–DeWitt chính xác) |
| Người quan sát ngoài: $\max\big(1-\lvert\langle\Psi\vert e^{-i\hat{\mathbb H}t}\vert\Psi\rangle\rvert^2\big)$ | **0**, tức vũ trụ **tĩnh** |
| Người quan sát trong: đọc $C$ ở phân cực góc φ, so trạng thái điều kiện của $S$ với $e^{-i\omega\sigma_yt}\lvert\psi_0\rangle$ tại $t=\varphi/\omega$ | Độ trùng khớp **1.000000** với mọi φ ∈ [0, π]. Hệ **tiến hóa đúng Schrödinger** |
| Độ phân giải của đồng hồ qubit | $\langle\varphi\vert\varphi'\rangle=\cos(\varphi-\varphi')$. Chỉ φ và φ+π/2 trực giao, nên đồng hồ phân biệt **tối đa 2 thời điểm**. φ và φ+π cho **cùng** trạng thái (chu kỳ π) |

**Kết luận phần 1:** thời gian nội tại xuất hiện đúng như cơ chế Page–Wootters. Nhưng một photon phân cực là đồng hồ **2 mức, tuần hoàn chu kỳ π**. Nó chỉ phân biệt được "trước" và "sau" trong vòng nửa chu kỳ.

## 2. Động cơ $\mathcal E_{\text{dual}}$ làm đồng hồ photon

**Ánh xạ:**
- Rô-to θ đóng vai photon đồng hồ.
- **Số đọc** là pha tích lũy Φ(t), tương tự góc phân cực đã quay qua.
- Hạt $x$ đóng vai photon hệ.

**Ba pha:**
- **Quá khứ:** khóa đồng bộ.
- **Hiện tại:** lần gãy mở đầu một đợt.
- **Tương lai:** lần tái khóa kết thúc đợt.

Với mỗi đợt, mình đọc Φ tại lúc gãy và lúc tái khóa, rồi kiểm xem đồng hồ có xếp đúng thứ tự không.

Mô phỏng: 8 quỹ đạo × $60T_0$ cho mỗi mức nhiễu, tham số theo Bảng I.

| $k_BT$ | Số đợt | Đồng hồ đọc "Tương lai **trước** Hiện tại" ($\Phi_{\text{khóa}}\le\Phi_{\text{gãy}}$) | Đợt kéo dài ≥ π (nửa vòng) | Đợt ≥ 2π (một vòng) | Sai số suy thời gian $\hat t=\Phi/\langle\dot\Phi\rangle$: RMS / max |
|---|---|---|---|---|---|
| 0.3 MeV | 37 | **0/37 (0%)** | 100% | 45.9% | 0.30 / 0.98 $T_0$ |
| **1 MeV** (Bảng I) | 66 | **0/66 (0%)** | **100%** | 66.7% | **0.58** / 1.82 $T_0$ |
| 5 MeV | 66 | **36/66 (54.5%)** | 97.0% | 78.8% | 14.05 / 65.08 $T_0$ |

**Đọc bảng:**

1. **Ở $k_BT\le1$ MeV (vùng hiệu lực của bản thảo):** đồng hồ $\mathcal E_{\text{dual}}$ xếp **đúng thứ tự** Quá khứ → Hiện tại → Tương lai ở **100%** số đợt. Sai số suy ra thời gian từ số đọc là ~0.6 $T_0$ (RMS). Con số này nhỏ hơn thời lượng một đợt (~1.2 $T_0$) và nhỏ hơn thời gian lưu (~4.5 $T_0$).
2. **Điều kiện bắt buộc là đồng hồ phải đếm vòng.**
   - 100% số đợt quay quá nửa vòng. Một đồng hồ phân cực kiểu Moreva (chu kỳ π) đọc cùng một giá trị ở hai thời điểm khác nhau trong đợt, nên **không xếp được thứ tự**.
   - 2/3 số đợt quay quá một vòng. Ngay cả số đọc góc modulo 2π cũng không đủ.
   - Vậy "đồng hồ photon" của học thuyết phải là **đồng hồ tích lũy** (Φ trên đường thẳng thực, hoặc kèm bộ đếm vòng), khác với photon phân cực đơn lẻ.
3. **Ở 5 MeV:** rô-to đảo chiều thường xuyên (đúng như kiểm chứng số cuộn trước đó). Đồng hồ đọc **Tương lai trước Hiện tại** ở **54.5%** số đợt, và sai số thời gian lên tới 14 $T_0$. Ở mức nhiễu này, đồng hồ $\mathcal E_{\text{dual}}$ **mất khả năng định hướng thời gian**. Điều này khớp với giới hạn $k_BT\lesssim1$ MeV trong bản thảo.

## 3. So sánh hai đồng hồ

| | Đồng hồ photon Moreva | Đồng hồ $\mathcal E_{\text{dual}}$ |
|---|---|---|
| Không gian số đọc | Góc phân cực, chu kỳ π, 2 trạng thái trực giao | Pha tích lũy Φ ∈ ℝ (liên tục, không tuần hoàn nếu đếm vòng) |
| Tiến hóa của hệ khi đọc đồng hồ | Unita, thuận nghịch, độ trùng khớp 1 | Có tiêu tán (Carnot, $\gamma_x$), không thuận nghịch |
| Ba pha Quá khứ → Hiện tại → Tương lai | Không có. Mọi thời điểm tương đương, có thể "chiếu ngược" | Có, và đồng hồ xếp đúng thứ tự ở $k_BT\le1$ MeV (100%) |
| Nguồn của mũi tên thời gian | Không có | Tiêu tán tích lũy $E_{\text{diss}}$ (không phải bản thân đồng hồ) |
| Giới hạn | Chỉ phân biệt 2 thời điểm trong nửa chu kỳ | Hỏng khi nhiễu mạnh (5 MeV: sai thứ tự 54.5%) |

## 4. Hai điểm trong bản thảo PDF nên chỉnh thêm

1. **Mục V.C.2 (Outlook) vẫn giữ câu cũ** "γₓ = 0.70ω₀ … corresponds physically to the strong surface-peaked absorption component … which dominates". Câu này mâu thuẫn với mục V.B.4 mới, vốn gọi 0.70 là "effective upper-bound baseline". Nên bỏ câu cũ ở V.C.2.
2. **Mục V.B.4** ghi "$r_{\text{rms}}\approx3.8$–$4.8$ fm cho $\gamma_{\text{eff}}\approx0.35$–$0.45$". Theo [wd-trung-binh-kiem-chung.md](wd-trung-binh-kiem-chung.md), các quỹ đạo 2s và 1d ($r_{\text{rms}}=4.78$ fm) chỉ cho 0.22–0.30. Khoảng đúng cho $r_{\text{rms}}=3.8$–4.8 fm là **0.22–0.45**. Riêng 0.35–0.45 chỉ đúng cho 0f và 0g.

Các số khác trong mục V.B.4 đều đúng: $\dot\theta_d/\dot\theta_{\text{br}}$ = 0.821 → 0.892, tốc độ đợt 0.156 → 0.231.

## Câu hỏi còn mở

- **Phiên bản lượng tử.** Có muốn làm phiên bản lượng tử của phần 2 không? Tức là trạng thái lịch sử Page–Wootters, trong đó môi trường (nơi nhận năng lượng Carnot) cũng nằm trong "vũ trụ tĩnh". Trạng thái điều kiện của hạt sẽ **mất độ tinh khiết** khi năng lượng thoát vào môi trường. Đây là cách chuẩn để đưa tính không thuận nghịch vào Page–Wootters, nhưng tốn tính toán hơn.
- **Đồng hồ phân cực có đếm vòng?** Nếu học thuyết muốn giữ hình ảnh "một photon phân cực", cần giải thích photon đó đếm vòng bằng cách nào. Ví dụ dùng bậc tự do mô-men quỹ đạo (OAM), hoặc một chuỗi photon liên tiếp.
