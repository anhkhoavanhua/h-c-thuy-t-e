# "Quá khứ = sóng, Hiện tại = hạt, ngưỡng = cả hai": giải số cho $\mathcal E_{\text{dual}}$

> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I, và đoạn diễn giải người dùng gửi (sóng ở $\mathcal E_{\text{past}}$, hạt ở $\mathcal E_{\text{present}}$, đồng thời ở ngưỡng $t_g$).
> Số liệu tính bằng [giai_so_song_hat.py](giai_so_song_hat.py). Lệnh `python3 giai_so_song_hat.py` chạy ≈ 75 s trên 8 lõi. Script chỉ import các module cũ, không sửa file nào.

## Vấn đề phương pháp

Trong PDF (7), **cả ba pha đều là cơ học cổ điển của chất điểm**. $x(t)$ và $\theta(t)$ xác định tại mọi lúc, kể cả trong pha khóa. "Sóng" ở pha Quá khứ là một **chuyển động tập thể tuần hoàn cổ điển**, không phải sóng vật chất. Vì vậy lưỡng tính sóng–hạt không thể đọc ra từ chính các phương trình (1)–(22).

Muốn kiểm, phải dùng các đại lượng lượng tử. Mình đo ba thứ:
- bước sóng de Broglie ở va chạm Carnot;
- trạng thái lượng tử của mô hình V.C.1 trong từng pha;
- mật độ mức so với Ericson.

## A. Hạt có định xứ tại điểm va chạm Carnot không?

Mô phỏng Langevin, tham số Bảng I, $k_BT=1$ MeV: 8 quỹ đạo × 60 $T_0$, 355 va chạm. Năng lượng mỗi va chạm $\Delta E=\tfrac12m_{\text{eff}}u^2$ trung bình **3.85 MeV**, trung vị 3.70 MeV.

Bước sóng de Broglie của chuyển động tương đối là $\lambda=h/(m_{\text{eff}}u)$:

| $m_{\text{eff}}$ | $\lambda$ trung vị | $\lambda$ (10% nhỏ nhất) | $\lambda/\kappa$ | $\lambda/R_{\text{hạt nhân}}$ | Va chạm có $\lambda<\kappa$ |
|---|---|---|---|---|---|
| 0.9 | **15.7 fm** | 9.9 fm | **23** | 3.2 | **0%** |
| 1.0 | 14.9 fm | 9.4 fm | 22 | 3.0 | 0% |

($\kappa=0.68$ fm là biên độ của cả vùng chuyển động, $R_{\text{hạt nhân}}=4.95$ fm.)

**Kết quả:**
- Bước sóng của "viên bi" lớn hơn **22 lần** vùng nó chuyển động và **3 lần** cả hạt nhân.
- Muốn $\lambda<\kappa$, mỗi va chạm phải mang **> 197 $\hbar\omega_0$ ≈ 2 GeV**, gấp ~500 lần giá trị thực.
- Vậy hình ảnh "định xứ tại điểm va chạm, giống viên bi va vào bia" **không đúng** ở thang này. Ở pha Hiện tại, hạt về mặt lượng tử là một sóng trải rộng hơn cả hạt nhân.

## B. Trạng thái lượng tử thực sự trong từng pha

Dùng mô hình MCWF của V.C.1: 16 quỹ đạo × 1500 đơn vị thời gian. Mình chia thời gian thành hai pha:
- **Quá khứ:** đã hơn $T_0$ kể từ bước làm nguội gần nhất.
- **Hiện tại:** trong vòng $T_0$ sau một bước làm nguội.

| Pha | Tỉ lệ thời gian | $\lvert\langle e^{i\theta}\rangle\rvert$ (độ xác định pha) | $\Delta\ell$ | $\Delta x$ của hạt | $\Delta x/\kappa$ | PR (số trạng thái dừng tham gia) |
|---|---|---|---|---|---|---|
| Quá khứ | 57% | **0.30** (90%: 0.49) | 1.08 | 0.718 $\ell_0$ | **2.15** | 3.8 |
| Hiện tại | 43% | **0.30** (90%: 0.49) | **0.72** | 0.713 $\ell_0$ | **2.14** | 3.1 |

**Cách đọc các cột:**
- $\lvert\langle e^{i\theta}\rangle\rvert$ = 1 nghĩa là pha rô-to hoàn toàn xác định, = 0 nghĩa là pha trải đều cả vòng.
- $\Delta x$ của trạng thái cơ bản dao động tử là 0.707 $\ell_0$.
- PR ≈ 1 nghĩa là trạng thái gần như là một sóng dừng (trạng thái riêng của $H$).

**Kết quả:**
1. **Quá khứ không phải "sóng có pha xác định".** Pha rô-to chỉ xác định 0.30: ba phần tư vòng tròn là bất định. Lý do là các trạng thái khóa gần với trạng thái riêng $\ell$, tức sóng dừng quanh vòng, không có góc pha cụ thể. "Pha đồng bộ xác định" là đặc tính của quỹ đạo cổ điển, không có ở mức ~3ħ (xem S3 trong [schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md): pha nhòe sau ≤ 0.4 $T_0$).
2. **Hiện tại không phải "hạt định xứ".** $\Delta x$ của hạt **giữ nguyên** 0.71 $\ell_0$ (2.1κ), bằng độ rộng điểm không, ở cả hai pha. Ngay sau bước nhảy, trạng thái còn **gần sóng dừng hơn** ($\Delta\ell$ 0.72 so với 1.08, PR 3.1 so với 3.8), vì bước làm nguội Davies chiếu về các trạng thái riêng thấp hơn.
3. **Hai pha gần như không phân biệt được về mặt sóng–hạt.** Khác biệt nhỏ duy nhất lại đi **ngược** hướng khẳng định: pha Hiện tại "sóng" hơn pha Quá khứ.
4. Năng lượng của bước nhảy nằm ở rô-to (99.7%, xem [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md)). Hạt gần như không bị "đập" chút nào.

## C. "Phổ liên tục" Ericson có đến từ mô hình không?

| $E$ (MeV) | Số mức của mô hình 2 bậc tự do ≤ E | Mật độ exciton Ericson $A=70$ (/MeV) | Khoảng cách mức của mô hình |
|---|---|---|---|
| 10 | 2 | 9.0·10⁷ | ~2.6 MeV |
| 20 | 5 | 1.8·10¹¹ | ~2.6 MeV |
| 40 | 15 | 1.3·10¹⁵ | ~1.9 MeV |

Mô hình chỉ có **15 mức** dưới 40 MeV, cách nhau khoảng 2 MeV. Hạt nhân A = 70 có ~10¹⁵ trạng thái mỗi MeV. Mật độ Ericson $\rho\propto E^{n-1}$ là tính chất của **70 nucleon**, không sinh ra từ rô-to và hạt. Bản thảo (V.A) **gán** Ericson từ bên ngoài, đúng như cách làm tiền cân bằng chuẩn. Câu "điều tiết bởi phổ liên tục của trường chân không dao động" không có cơ sở trong mô hình.

## D. Còn "ngưỡng = đồng thời cả hai"?

| Ý trong đoạn diễn giải | Trạng thái |
|---|---|
| Tọa độ và động lượng liên tục tại $t_g$, gia tốc nhảy bậc (Eq. 11) | ✓ Đúng, đã kiểm cổ điển ($\Delta p\to0$, $\Delta\ddot x=-\lambda/m$). Nhưng đây là tính chất của **ràng buộc một phía trong cơ học cổ điển**, không phải lưỡng tính sóng–hạt |
| Ở ngưỡng hệ "vừa sóng vừa hạt" | Phiên bản lượng tử đúng là một **chồng chập** giữa thành phần còn khóa và thành phần đã thoát. Xác suất thoát tăng **liên tục** theo ℓ: 1.4% ở ℓ = 3, 14% ở ℓ = 6 ([luong-tu-schrodinger-dirac-pw-wdw.md](luong-tu-schrodinger-dirac-pw-wdw.md), S3). Đó là chồng chập "khóa / gãy", không phải chồng chập "sóng / hạt" |
| Mỗi va chạm nhả gói 4.0–4.6 MeV | Cổ điển: TB 3.85 MeV ở lần chạy này. Lượng tử: gói thật là **20–35 MeV** theo bậc thang rô-to ([mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md)) |

## Kết luận

| Khẳng định | Kết quả số |
|---|---|
| Quá khứ = sóng có pha xác định | ✗ Trạng thái lượng tử là sóng dừng **không có** pha xác định ($\lvert\langle e^{i\theta}\rangle\rvert=0.30$). Pha xác định chỉ có trong quỹ đạo cổ điển |
| Hiện tại = hạt định xứ tại điểm va chạm | ✗ $\lambda_{\text{dB}}\approx15$ fm = 22κ = 3 $R_{\text{hạt nhân}}$; $\Delta x$ không đổi (2.1κ) |
| Ngưỡng = đồng thời cả hai | Phần đúng là sự liên tục của $(x,p)$ (cổ điển) và chồng chập khóa/gãy (lượng tử). Không phải lưỡng tính sóng–hạt |
| Exciton được điều tiết bởi phổ liên tục chân không | ✗ Mô hình có 15 mức dưới 40 MeV; mật độ Ericson đến từ 70 nucleon và được gán từ ngoài |

**Cách diễn đạt có cơ sở:** ba pha khác nhau ở **số bậc tự do và tính khả nghịch**, không ở tính sóng hay hạt:
- **Quá khứ:** 1 bậc tự do, bảo toàn.
- **Hiện tại:** 2 bậc tự do, có va chạm tiêu tán.
- **Tương lai:** trở lại 1 bậc tự do, với năng lượng thấp hơn.

Ở cả ba pha, hạt về mặt lượng tử đều là sóng có độ rộng ≈ 2κ. Nếu muốn đưa hình ảnh sóng–hạt vào bài báo, nên dùng đúng chỗ này: **ở mức ~3ħ, mô hình cổ điển là một bản ghi đo liên tục** ([pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md)), không phải ảnh chụp của một hạt định xứ.
