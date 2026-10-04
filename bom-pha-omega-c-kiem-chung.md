# Bơm pha lượng tử và tần số cắt $\hbar\omega_c\gtrsim50$ MeV: giải và kiểm chứng số liệu

> Tiếp nối [u-hermite-drude-kiem-chung.md](u-hermite-drude-kiem-chung.md).
> Số liệu tính bằng [giai_so_bom_pha.py](giai_so_bom_pha.py) (cần numpy). Chạy `python3 giai_so_bom_pha.py`, ≈ 1 phút trên 8 lõi.
> **Đơn vị:** $\hbar=m=\omega_0=1$ ($\hbar\omega_0=10$ MeV), $I=1$, $\kappa=1/3$, $\Delta E=1$, $g=0.1$, $k_BT=1$ MeV.
> **Cơ sở tính:** $\ell\in[-8,8]$, $n\le6$.

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\mathcal D_{\text{pump}}=\Gamma_{\text{pump}}\big(e^{i\theta}\rho e^{-i\theta}-\tfrac12\{e^{-i\theta}e^{i\theta},\rho\}\big)$ với $\Gamma_{\text{pump}}(\ell)=\Gamma_0\max(0,1-\ell^2/\ell_0^2)$ | **Viết như vậy thì không có bão hòa.** $\Gamma_{\text{pump}}(\ell)$ là hệ số vô hướng, không phải toán tử. Vì $e^{-i\theta}e^{i\theta}=1$, bộ sinh chỉ còn $\Gamma(e^{i\theta}\rho e^{-i\theta}-\rho)$, không phụ thuộc $\ell$. Muốn có bão hòa, phải đưa $\ell$ vào **toán tử nhảy**: $L_p=\sqrt{\Gamma_0}\,e^{i\theta}f(\hat\ell)$, với $f(\ell)=\sqrt{\max(0,1-\ell^2/\ell_0^2)}$. Đã mô phỏng cả hai dạng |
| Bơm + làm nguội Davies $\hat u_-$ tái lập chu trình tuần hoàn | **✓ Có chu trình dừng.** Bơm và làm nguội cân bằng nhau: ⟨ℓ⟩ dừng 0.65–2.38, có đợt làm nguội lặp lại. Số bước mỗi đợt **có dạng hình học** (CV đo khớp CV hình học cùng trung bình: 0.57/0.52, 0.37/0.39, 0.54/0.60) |
| …và "tương đương cơ chế Van der Pol cổ điển" | **Không khớp định lượng.** Số bước mỗi đợt TB chỉ **1.2–1.6** (cổ điển **5.1**). Năng lượng mỗi bước làm nguội TB **0.57–2.05 $\hbar\omega_0$** (cổ điển **0.40**). Ở cấu hình có rô-to quay gần tốc độ mô hình (⟨ℓ⟩ = 2.4): **2.05 $\hbar\omega_0$**, gấp 5 lần. Thời gian chờ giữa các đợt: CV 0.66–1.56, không nhất quán dạng mũ |
| Cần $\hbar\omega_c\gtrsim50$–$100$ MeV để giữ $\gamma_x\approx0.7$ trên dải $\omega_0\to\Omega$ | **Đúng một phần.** Tỉ lệ ma sát giữ lại tại $\Omega=2.63$: 50 MeV giữ **78%**, 100 MeV giữ **94%**. Giữ ≥ 90% cần **$\hbar\omega_c\ge79$ MeV**. Vậy 50 MeV chưa đủ, còn ~80–100 MeV thì đủ |
| "$\hbar\omega_{\text{cut}}\approx m_\sigma/4\approx50$–$100$ MeV" | **Sai số học và thiếu căn cứ.** Với $m_\sigma c^2=400$–500 MeV thì $m_\sigma/4=$ **100–125 MeV**, không phải 50–100. Hệ số /4 cũng không được suy ra |
| Trao đổi $2\pi$ tạo phổ liên tục với tần số cắt 50–100 MeV | **Không khớp ngưỡng.** Phổ liên tục $2\pi$ bắt đầu từ **$2m_\pi c^2=279$ MeV**. Nếu lấy thang đó làm tần số cắt thì ma sát giữ lại 99.1% |
| Vùng quasi-deuteron (40–100 MeV) | Thang năng lượng có thật, nhưng đó là **hấp thụ photon** bởi cặp pn. Việc dùng nó làm **bể ma sát** cho tọa độ $x$ của nucleon là một **giả thuyết mới**, chưa được suy ra |
| Hệ quả chưa nêu: tần số cắt cao làm hạt bị "mặc áo" | Cân bằng chính xác ở 1 MeV: $\langle p^2\rangle$ tăng **0.79 → 0.90 → 1.09** khi $\hbar\omega_c$ = 50 → 100 → 280 MeV (so với 0.5 khi không ghép). Tức là động năng điểm không của hạt tăng 60–120%. Phần này cần đưa vào ngân sách năng lượng |

## 1. Tần số cắt

| $\hbar\omega_c$ (MeV) | Giữ lại tại $\omega_0$ | Giữ lại tại $\Omega=2.63$ |
|---|---|---|
| 4.5 (độ rộng GDR) | 16.8% | 2.8% |
| 25 | 86.2% | 47.5% |
| 50 | 96.2% | **78.3%** |
| **79** | 98.4% | **90.0%** |
| 100 | 99.0% | 93.5% |
| 280 ($2m_\pi$) | 99.9% | 99.1% |

**Cân bằng chính xác** (chuyển động Brown lượng tử, $k_BT=1$ MeV, $\gamma=0.7$):

| $\hbar\omega_c$ | $\langle x^2\rangle$ | $\langle p^2\rangle$ | Tích |
|---|---|---|---|
| 50 MeV | 0.4349 | 0.7942 | 0.345 |
| 100 MeV | 0.4288 | 0.9037 | 0.388 |
| 280 MeV | 0.4237 | 1.0904 | 0.462 |

## 2. Bơm pha + làm nguội Davies (MCWF)

**Thiết lập:** 16 quỹ đạo × 1500 đơn vị thời gian (≈ 239 $T_0$), bỏ 20% đầu. Làm nguội: $L_c=\sqrt{0.5}\,\hat u_-$. Một "đợt" là các bước làm nguội cách nhau dưới $T_0$.

| Cấu hình | ⟨ℓ⟩ dừng | Nhảy bơm / làm nguội (mỗi $T_0$) | Năng lượng mỗi bước làm nguội: TB (10–90%) | Bước/đợt: TB, CV (CV hình học) | Phân bố bước/đợt | Chờ giữa đợt: TB, CV |
|---|---|---|---|---|---|---|
| Nguyên văn (không bão hòa), $\Gamma_0=0.05$ | 0.88 | 0.34 / 0.20 | 0.85 (−0.80–2.37) | 1.37, 0.57 (0.52) | 1:321, 2:76, 3:15, 4:12, 5:3, 6:1 | 6.33 $T_0$, 1.56 |
| Bão hòa, $\Gamma_0=0.05$ | 0.65 | 0.29 / 0.13 | **0.57** (−0.66–1.63) | 1.18, 0.37 (0.39) | 1:272, 2:48, 3:6 | 8.36 $T_0$, 1.23 |
| Bão hòa, $\Gamma_0=0.2$ | **2.38** | 0.68 / 0.52 | **2.05** (0.79–3.17) | 1.57, 0.54 (0.60) | 1:607, 2:255, 3:98, 4:25, 5:8, 7:1 | 2.69 $T_0$, 0.66 |
| **Cổ điển (Part I)** | ~2.6 | — | **0.40** | **5.12**, 0.90 | hình học $p=0.195$ | ≈ 4.6 $T_0$, ≈ 1 |

**Đọc bảng:**

1. **Có chu trình lượng tử**, và dạng hình học của số bước mỗi đợt xuất hiện tự nhiên. Đây là điểm khớp định tính.
2. **Có đánh đổi không thoát được:**
   - Muốn rô-to quay gần tốc độ mô hình (ℓ ≈ 2.4) thì mỗi bước làm nguội mất ~2 $\hbar\omega_0$. Lý do: chuyển dời giữa các mức riêng của $H$ khi ℓ ~ 2–3 cách nhau cỡ đó.
   - Muốn năng lượng mỗi bước nhỏ (0.57) thì rô-to gần như đứng yên (ℓ ≈ 0.65).
   - Không cấu hình nào đạt đồng thời ℓ ≈ 2.6 **và** 0.40 $\hbar\omega_0$ mỗi bước. Đây là hệ quả trực tiếp của việc rô-to chỉ mang ~3ħ.
3. **Số bước mỗi đợt nhỏ hơn cổ điển 3–4 lần**, vì mỗi bước lượng tử lấy đi nhiều năng lượng hơn.
4. **Một số bước "làm nguội" có $\Delta\langle H\rangle<0$** (10% dưới −0.66). Nguyên nhân: trạng thái là chồng chập nhiều mức riêng. Bước nhảy Davies hạ mức riêng nhưng có thể làm tăng kỳ vọng năng lượng của gói sóng.
5. **Dạng nguyên văn (không bão hòa) không chạy vọt ra biên cơ sở** trong mô phỏng, vì làm nguội đã cân bằng bơm. Tuy vậy, phần bão hòa là cần thiết nếu muốn kiểm soát ⟨ℓ⟩ độc lập với tốc độ làm nguội.

## 3. Đề xuất cập nhật đoạn Outlook cho Part I

Thêm vào đoạn Outlook ở [u-hermite-drude-kiem-chung.md](u-hermite-drude-kiem-chung.md):

> "(iv) Adding an incoherent, saturating phase pump $L_p=\sqrt{\Gamma_0}\,e^{i\theta}f(\hat\ell)$ to the Davies cooling channel produces a stationary quantum limit cycle with geometrically distributed relaxation bursts. However, the bursts contain only 1.2–1.6 jumps on average, and in the regime where the rotor turns at the classical rate ($\langle\ell\rangle\approx2.4$) each jump removes $\approx2\,\hbar\omega_0$, compared with $\approx0.4\,\hbar\omega_0$ per classical impact. (v) Retaining $\ge90\%$ of the classical friction at the rotor frequency requires a bath cutoff $\hbar\omega_c\gtrsim80$ MeV; identifying the corresponding physical continuum (e.g., quasi-deuteron or meson-exchange modes) remains an open question."

## Câu hỏi còn mở

- **Nguồn thật của bể có $\hbar\omega_c\sim80$–280 MeV.** Cần một phép tính độc lập, ví dụ phần ảo của thế quang học ở vùng 20–30 MeV, để ra hàm phổ $J(\omega)$ của nucleon trong hạt nhân. Chưa nên gán theo quasi-deuteron hay $m_\sigma$.
- **Có nên xem sai lệch lượng tử–cổ điển là kết quả chính của Part II không?** Các số liệu ở đây cho thấy với rô-to ~3ħ, thống kê đợt (1.2–1.6 bước, ~2 $\hbar\omega_0$ mỗi bước) khác hẳn cổ điển (5.1 bước, 0.4 $\hbar\omega_0$). Đây có thể là **dự đoán** của Part II, thay vì một khuyết điểm phải khắc phục.
