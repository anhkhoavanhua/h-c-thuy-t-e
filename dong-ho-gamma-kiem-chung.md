# Kiểm kiến trúc "đồng hồ hạt nhân gamma" Page–Wootters cho $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), A ~ 70, $\hbar\omega_0$ = 10 MeV, $T_0$ = 4.14·10⁻²² s. Nối tiếp [dong-ho-hat-nhan-page-wootters-kiem-chung.md](dong-ho-hat-nhan-page-wootters-kiem-chung.md).
> **Số liệu:** [giai_so_dong_ho_gamma.py](giai_so_dong_ho_gamma.py), chạy ≈ 3 phút. Phần G4 dùng mô phỏng Langevin Part I.

## Tóm tắt

| Đề xuất | Kết quả số | Kết luận |
|---|---|---|
| GDR làm đồng hồ ($\Delta E_C$ ≈ 3–11 MeV) | GDR sống **0.27–0.40 $T_0$**, tắt trước khi dao động hết **một** chu kỳ (0.47–0.70 chu kỳ) | ✗ Không làm đồng hồ được |
| Chuỗi xung gamma τ ≤ $T_0$/25 do hạt nhân tự phát | Poisson: cần 3200 photon trong 4 $T_0$ (≥ 3.2 GeV). Chuỗi đều lý tưởng: ≥ 183 MeV. Cả hai **vượt E\* = 20–30 MeV khoảng 6–160 lần** | ✗ Hạt nhân không đủ năng lượng tự làm đồng hồ |
| Dải quay hạt nhân | Tích tắc $T_0$/25 cần **J ≈ 600–1500**. Hạt nhân A ~ 70 vỡ trước J ~ 40–60 | ✗ |
| "Bắt trọn 100% chuỗi va chạm Carnot" với σ_t = 0.07 $T_0$ | Mức **đợt**: 99%. Mức **từng va chạm**: 88.5%. 16% cặp va chạm cách nhau < 0.01 $T_0$ | ✗ Chỉ đúng ở mức đợt |
| "Pha Quá khứ: liên kết cứng, tinh khiết cao, tiến hóa unita" | Bản lượng tử: liên kết nhòe ≥ 1.57κ; Tr ρ² ≈ 0.25 ở trạng thái dừng; có bơm và bể nên không unita | ✗ Chỉ là hình ảnh cổ điển |
| "ΔE ≈ 0.412 $\hbar\omega_0$ tiêu tán" ở pha Tương lai | 0.412 là năng lượng **mỗi va chạm** (mô hình chặn gãy p\* = 0.3), không phải mỗi đợt (cả đợt ≈ 1.3–2.5 $\hbar\omega_0$) | ✗ Nhầm đơn vị |
| Cầu nối ²²⁹Th: "E_c là hiệu hai đại lượng hàng chục MeV" | Phép triệt tiêu là giữa **hiệu** năng lượng Coulomb và hạt nhân của hai trạng thái, mỗi hiệu cỡ **0.1–0.5 MeV**, không phải hàng chục MeV. K ≈ 1.2·10⁴–6·10⁴ | ✓ Ý đúng, cần sửa số |

## G1. GDR làm đồng hồ

$E_G = 31.2A^{-1/3}+20.6A^{-1/6}$ = 17.7 MeV với A = 70; độ rộng thực nghiệm Γ ≈ 4–6 MeV.

| Γ (MeV) | Thời gian sống ħ/Γ | = $T_0$ | Số chu kỳ trước khi tắt | Q = E_G/Γ |
|---|---|---|---|---|
| 4 | 1.65·10⁻²² s | 0.40 | 0.70 | 4.4 |
| 5 | 1.32·10⁻²² s | 0.32 | 0.56 | 3.5 |
| 6 | 1.10·10⁻²² s | 0.27 | 0.47 | 3.0 |

Đồng hồ cần ≥ 25 tích tắc mỗi $T_0$ và giữ được qua ≥ 5 $T_0$, tức ≥ 125 tích tắc kết hợp. GDR tắt trước khi dao động hết một chu kỳ. **Trải năng lượng rộng (Γ lớn) là điều kiện cần cho độ phân giải, nhưng nó cũng là nguyên nhân khiến GDR tắt nhanh.** Không mode tập thể tắt nhanh nào làm đồng hồ được.

## G2. Chuỗi xung gamma phát từ chính hạt nhân

Mục tiêu: σ_t = 0.07 $T_0$ (99% xếp đúng 3 pha) trên một chu trình khoảng 4 $T_0$.

| Cách phát | Số photon | Năng lượng tối thiểu | So với E\* = 20–30 MeV |
|---|---|---|---|
| Poisson (phân rã gamma thật) | k ≥ (t/σ)² = **3192** | ≥ 3.2 GeV (1 MeV/photon) | thiếu **106–160 lần** |
| Chuỗi đều lý tưởng (τ ≤ √12 σ = 0.245 $T_0$) | 16 | 16 × ħ/2σ (11.2 MeV) = **183 MeV** | thiếu **6–9 lần** |

Thêm vào đó, hạt nhân phát gamma theo kiểu **Poisson**, không đều. Một chuỗi đều cần một nguồn kết hợp **bên ngoài** hạt nhân. **Một hạt nhân không thể tự mang một đồng hồ đủ nhanh cho chính nó.**

## G3. Dải quay

Tích tắc $T_0$/25 đòi $\hbar\omega_{\text{rot}}=dE/dJ\ge$ **250 MeV**.

| E(2⁺) (MeV) | ħ²/2𝓘 (MeV) | J cần | Chu kỳ quay ở J = 20 |
|---|---|---|---|
| 0.5 | 0.083 | 1499 | 2.9 $T_0$ |
| 1.0 | 0.167 | 749 | 1.5 $T_0$ |
| 1.2 | 0.200 | 624 | 1.2 $T_0$ |

Hạt nhân vỡ trước J ~ 40–60. Hơn nữa, một dải quay là **đồng hồ tuần hoàn**, nên nó vướng đúng vấn đề đếm vòng đã chứng minh ở Part III (sàn sai số do trộn vòng).

## G4. "Bắt trọn 100% va chạm Carnot"?

Từ Langevin Part I (851 cặp va chạm liên tiếp trong đợt): khoảng cách trung vị 0.25 $T_0$, phân vị 10% là 0.0005 $T_0$; **16%** số cặp cách nhau dưới 0.01 $T_0$.

- Với σ_t = 0.07 $T_0$: **88.5%** số cặp va chạm được xếp đúng thứ tự, và chỉ 69.6% cách nhau đủ xa (> 2σ).
- **99% chỉ đúng ở mức đợt** (Quá khứ / Hiện tại / Tương lai), không đúng ở mức từng va chạm.

## G5. Cầu nối với đồng hồ ²²⁹Th

Chuyển dời 8.36 eV nhỏ vì có sự triệt tiêu giữa **hiệu** năng lượng Coulomb ΔV_C và **hiệu** năng lượng hạt nhân của đồng phân so với trạng thái cơ bản. Mỗi hiệu chỉ cỡ 0.1–0.5 MeV, không phải hàng chục MeV. Hệ số khuếch đại K ≈ |ΔV|/E_c:

| \|ΔV\| (MeV) | K | δE_vt/E_vt nhỏ nhất đo được (đồng hồ chính xác 10⁻¹⁸) |
|---|---|---|
| 0.1 | 1.2·10⁴ | 8.4·10⁻²³ |
| 0.3 | 3.6·10⁴ | 2.8·10⁻²³ |
| 0.5 | 6.0·10⁴ | 1.7·10⁻²³ |

(|ΔV| ở đây là khoảng giá trị minh họa cho bậc độ lớn. Con số chính xác lấy từ các phân tích độ nhạy hằng số cấu trúc tinh tế của ²²⁹Th.)

- **Độ nhạy cực cao** (≈ 10⁻²³), nhưng chỉ với **đúng loại tín hiệu**:
  - dao động của $E_{\text{vt}}$ ở tần số động cơ (2.4·10²¹ Hz) bị lấy trung bình trong 1 s hỏi đồng hồ, chỉ còn **≤ 10⁻²²** lần;
  - một độ dịch **hằng** của $E_{\text{vt}}$ bị hấp thụ vào định nghĩa ν₀, nên không quan sát được tuyệt đối.
- **Đồng hồ Th chỉ kiểm được biến thiên chậm** (≲ kHz) của $E_{\text{vt}}$ theo thời gian, theo vị trí, hoặc theo môi trường. Học thuyết hiện **chưa** dự đoán biến thiên loại này.

## Kiến trúc khả thi rút ra từ số liệu

1. **Không có đồng hồ vật lý nào nằm trong hạt nhân** mà phân giải được ba pha ở 99%. Mọi ứng viên (GDR, chuỗi gamma, dải quay) đều thua vì **ngân sách năng lượng** (E\* ≈ 20–30 MeV) hoặc **thời gian sống**.
2. Trong Page–Wootters, đồng hồ τ = $T_0$/25 vẫn **hợp lệ như một thiết bị lý tưởng**: ràng buộc WdW thỏa tới 4·10⁻¹⁵ ở Part III. Nhưng cần nói rõ nó không phải một hạt nhân gamma cụ thể.
3. **Đồng hồ nội tại thật sự có sẵn là chính bộ đếm bản ghi Carnot** (Part III): mỗi bản ghi là một tích tắc không thuận nghịch. Tốc độ ≈ 0.5 bản ghi/$T_0$, nên độ phân giải ≈ 1–2 $T_0$, đủ xếp đúng khoảng 70% ở mức đợt (bảng N1 của bài trước). Đây là **đồng hồ nhiệt động**, thô nhưng có thật, và nó đồng thời là mũi tên thời gian.
4. **Đồng hồ ²²⁹Th** là công cụ **thực nghiệm bên ngoài** cho các đại lượng dừng hoặc biến thiên chậm, với độ nhạy δE_vt/E_vt ~ 10⁻²³. Cần học thuyết đưa ra dự đoán biến thiên chậm của $E_{\text{vt}}$.

## Câu hỏi còn mở
- Học thuyết có cơ chế nào làm $E_{\text{vt}}$ biến thiên **chậm** (theo thời gian vũ trụ, trường hấp dẫn, hay môi trường) không? Nếu có, đồng hồ Th kiểm được tới 10⁻²³. Nếu không, Th không kiểm được học thuyết.
