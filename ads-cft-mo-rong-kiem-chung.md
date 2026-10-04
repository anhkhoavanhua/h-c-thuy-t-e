# AdS/CFT cho $\mathcal E_{\text{dual}}$, phần mở rộng: ngưỡng gãy, độ sâu xuyên tâm, Hawking–Page

> **Cơ sở:** bản thảo [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I. Nối tiếp [ads-cft-toan-anh-kiem-chung.md](ads-cft-toan-anh-kiem-chung.md), với nền $L_{\text{AdS}}=\hbar c/\hbar\omega_0=19.73$ fm $=9.69\,\ell_0$.
> Số liệu tính bằng [giai_so_ads_cft_mo_rong.py](giai_so_ads_cft_mo_rong.py). Lệnh `python3 giai_so_ads_cft_mo_rong.py F G` chạy tức thì. `H` ≈ 4 phút, `K` ≈ 3 phút trên 8 lõi. Phần mô phỏng gọi lại `_traj` của [giai_so_dong_ho_photon.py](giai_so_dong_ho_photon.py) và `run` của module Langevin chung, **không sửa** file cũ nào.

## F. Ngưỡng gãy $R$ là lực thoát khỏi giếng AdS

Trong giới hạn phi tương đối tính, hấp dẫn AdS kéo hạt về tâm với lực $F(r)=m\omega_0^2r=kr$. Đây chính là số hạng $kx$ trong bản thảo.

| Đại lượng | Giá trị |
|---|---|
| Lực AdS tại $r=\ell_0$: $k\ell_0=\hbar\omega_0/\ell_0$ | 4.911 MeV/fm = 786.8 N |
| $R=E_{\text{vt}}/\ell_0$ | 2360 N $=\mathbf{3.0\times k\ell_0}$ |
| Bán kính tương ứng $r_{\text{br}}=R/k=3\ell_0=9\kappa$ | **6.11 fm** $=0.31\,L=1.24\,R_{\text{hạt nhân}}$ |

**Đọc theo AdS:** liên kết gãy khi ứng suất ly tâm vượt lực hấp dẫn AdS tại bán kính $3\ell_0$, tức ngay ngoài mặt hạt nhân. Điều kiện tự đồng dạng $E_{\text{vt}}/\Delta E=3$ vì vậy có nghĩa hình học: "tường chân không" nằm ở $r=3\ell_0$.

**Tính lại $\dot\theta_{\text{br}}$** theo cách của bản thảo: lấy $\dot\theta$ tại $\theta=0$, cho $\dot\theta(\theta)$ biến thiên theo bảo toàn năng lượng trên đa tạp ràng buộc, rồi dùng Eq. (7).

| Mô hình | $\dot\theta_{\text{br}}$ ($\omega_0$) | Lệch |
|---|---|---|
| Bản thảo (phi tương đối tính, thế điều hòa) | **3.2048** (bản thảo 3.205 ✓) | |
| + thế tĩnh AdS chính xác $V=kL^2(\sqrt{1+x^2/L^2}-1)$ | 3.2047 | −8·10⁻⁵ |
| + động học tương đối tính ($\gamma^3m$, $v_{\max}/c=0.11$) | 3.2021 | **−0.08%** |

Kết quả Part I của bản thảo đứng vững khi nhúng vào AdS. Hiệu chỉnh lớn nhất là ~0.1%, đến từ động học tương đối tính chứ không từ độ cong của AdS.

## G. Độ sâu xuyên tâm (UV/IR)

Mỗi mức năng lượng dao động $E$ ứng với bán kính quay đầu $r_{\max}=\ell_0\sqrt{2E/\hbar\omega_0}$ trong AdS.

| Trạng thái | $E$ (MeV) | $r_{\max}$ (fm) | $r_{\max}/L$ |
|---|---|---|---|
| Một va chạm Carnot $\langle\Delta E\rangle$ | 4.3 | 1.89 | 0.096 |
| Điểm không $\hbar\omega_0/2$ | 5.0 | 2.04 | 0.103 |
| Hạt bị kéo ở ngưỡng, $\tfrac12m(\kappa\dot\theta_{\text{br}})^2$ | 5.71 | 2.18 | 0.110 |
| $S_n$ (bắt đầu phát nucleon) | 8.0 | 2.58 | 0.131 |
| Trần năng lượng đợt (V.B.3) | 20.4 | 4.11 | 0.208 |
| $E_{\text{vt}}$ | 30.0 | **4.99** | 0.253 |

- Toàn bộ động lực học Past → Present → Future diễn ra ở $r\le0.25\,L$. Đây là vùng hồng ngoại sâu của AdS, nơi hình học gần như phẳng cộng một thế điều hòa. Điều này giải thích vì sao hiệu chỉnh AdS ở mục F chỉ cỡ 10⁻⁴.
- $E_{\text{vt}}=30$ MeV ứng với $r_{\max}=4.99$ fm, **trùng bán kính hạt nhân** $1.2A^{1/3}=4.95$ fm (lệch 0.8%). Ở biên độ này, hạt dao động chạm đúng mặt hạt nhân.

## H. Hawking–Page trong AdS₄

Lỗ đen AdS₄ có $T(r_h)=(1+3r_h^2/L^2)/(4\pi r_h)$ và năng lượng tự do $F\propto r_h(1-r_h^2/L^2)$. Giải số:

| | Không thứ nguyên ($\hbar\omega_0$) | MeV | Giải tích |
|---|---|---|---|
| $T_{\min}$ (dưới mức này không có lỗ đen) | 0.27566 | **2.76** | $\sqrt3/2\pi$ ✓ |
| $T_{\text{HP}}$ (lỗ đen thắng AdS nhiệt) | 0.31831 | **3.18** | $1/\pi$ ✓ |

$T_{\text{HP}}=\hbar\omega_0/\pi$ **không có tham số tự do** khi đã cố định $L=\hbar c/\hbar\omega_0$. Trong AdS/QCD, chuyển pha Hawking–Page ứng với chuyển pha giải giam.

### H′. Quét $k_BT$ bằng mô phỏng Langevin (Bảng I, 8 quỹ đạo × 50 $T_0$ mỗi mức)

| $k_BT$ (MeV) | $k_BT/T_{\text{HP}}$ | Đợt | Đồng hồ đọc "Tương lai trước Hiện tại" | Sai số $\hat t$ trung vị ($T_0$) |
|---|---|---|---|---|
| 1.0 | 0.31 | 43 | **0%** | 0.63 |
| 2.0 | 0.63 | 50 | 18.0% | 0.67 |
| 2.5 | 0.79 | 49 | 16.3% | 1.12 |
| 3.0 | 0.94 | 29 | **62.1%** | **10.3** |
| 3.5 | 1.10 | 40 | 55.0% | 15.2 |
| 4.0 | 1.26 | 40 | 37.5% | 16.4 |
| 5.0 | 1.57 | 43 | 58.1% | 7.0 |

Nhìn riêng bảng này thì đẹp. Sai số thời gian nhảy từ ~1 lên ~10 $T_0$ giữa 2.5 và 3.0 MeV, đúng khoảng $T_{\min}$–$T_{\text{HP}}$ = 2.76–3.18 MeV. Trên mức đó, tỉ lệ đợt đảo thứ tự dao động quanh ~50%, tức đồng hồ đoán ngẫu nhiên.

### K. Kiểm xem có phải trùng hợp không

Mình đổi các tham số **không liên quan đến** $\omega_0$. $T_{\text{HP}}$ giữ nguyên 3.18 MeV trong cả hai trường hợp:

| Cấu hình | $k_BT$ = 1 MeV | 2 MeV | 3 MeV | 4 MeV |
|---|---|---|---|---|
| Bảng I (mục H′) | 0% / 0.63 | 18% / 0.67 | 62% / 10.3 | 38% / 16.4 |
| $\gamma_x=0.35$ | 0% / 0.28 | 0% / 0.57 | 18% / 1.03 | 28% / **30.8** |
| $\Delta E=0.5$, $E_{\text{vt}}=1.5$ ($\dot\theta_{\text{br}}=2.34$) | **40% / 62.0** | 49% / 26.4 | 77% / 40.7 | 60% / 15.2 |

*(Mỗi ô: tỉ lệ đợt "Tương lai trước Hiện tại" / sai số $\hat t$ trung vị theo $T_0$.)*

**Kết luận: trùng hợp.**
- Điểm mất mũi tên thời gian **dịch theo $\Delta E$ và $\gamma_x$**. Ma sát giảm một nửa thì điểm vỡ lên 3–4 MeV. Rào thế giảm một nửa thì đã vỡ ngay ở 1 MeV.
- Trong khi đó $T_{\text{HP}}$ đứng yên. Do đó việc Bảng I vỡ gần $T_{\text{HP}}$ là do cặp $\Delta E=\hbar\omega_0$, $\gamma_x=0.70$ tình cờ đặt nó ở đó, không phải hệ quả của AdS.
- **Không nên** viết trong bản thảo rằng "giới hạn $k_BT\lesssim1$ MeV là pha giam Hawking–Page". Câu đúng: mũi tên thời gian của động cơ được bảo vệ bởi rào $\Delta E$ và ma sát $\gamma_x$. Nhiệt độ Hawking–Page chỉ tình cờ nằm cùng thang.

## Tổng hợp: cái gì đứng, cái gì không

| Khẳng định | Trạng thái | Nguồn |
|---|---|---|
| $L_{\text{AdS}}=\hbar c/\hbar\omega_0=19.73$ fm tái lập $T_0$, $\ell_0$ | **Đứng** (chính xác) | Bài trước, A |
| Chu kỳ không phụ thuộc biên độ, phổ $\Delta+2n+l$, bậc suy biến vỏ, $A\sim70$ | **Đứng** (sai số ≤ 10⁻⁸) | Bài trước, B–C |
| $R$ = lực AdS tại $r_{\text{br}}=3\ell_0=6.11$ fm | **Đứng** (đại số, đọc hình học mới) | F |
| $\dot\theta_{\text{br}}=3.205$ bền khi nhúng vào AdS | **Đứng** (lệch ≤ 0.08%) | F |
| $E_{\text{vt}}$ ↔ biên độ chạm mặt hạt nhân (4.99 vs 4.95 fm) | **Trùng hợp theo $A$** (chỉ khớp quanh $A\approx75$) | G |
| Giới hạn nhiệt của mũi tên thời gian = Hawking–Page | **Trùng hợp** (điểm vỡ dịch theo $\Delta E$, $\gamma_x$) | H′, K |
| Cận Bekenstein cho đợt Carnot | Chỉ thỏa khi tính khối lượng nghỉ | Bài trước, E |
| Đối ngẫu toàn ảnh thật sự cho động cơ | **Không** (không có CFT N lớn; rô-to phá SO(d,2)) | Bài trước, F |

## Đoạn đề xuất thêm vào bản thảo (tiếng Anh)

> In the AdS identification, the breakdown threshold acquires a geometric meaning: since the non-relativistic AdS restoring force is $F(r)=m\omega_0^2r$, the condition $|\lambda|\ge R=E_{\text{vt}}/\ell_0$ corresponds to exceeding the confining force at $r_{\text{br}}=3\ell_0\approx6.1$ fm, just outside the nuclear surface. Embedding the constrained dynamics in exact global AdS with relativistic kinematics shifts $\dot\theta_{\text{br}}$ by only $-0.08\%$ (3.2048 → 3.2021 $\omega_0$).

### Kiểm trùng hợp của mục G theo số khối $A$

Với $\hbar\omega_0=41A^{-1/3}$ MeV và $E_{\text{vt}}=3\hbar\omega_0$, ta có $r_{\max}(E_{\text{vt}})=\sqrt6\,\ell_0=2.46\,A^{1/6}$ fm, còn $R_{\text{hạt nhân}}=1.2\,A^{1/3}$ fm. Hai lũy thừa khác nhau, nên chúng chỉ bằng nhau tại **$A\approx75$**:

| $A$ | $\hbar\omega_0$ (MeV) | $r_{\max}(E_{\text{vt}})$ (fm) | $R_{\text{hạt nhân}}$ (fm) |
|---|---|---|---|
| 16 | 16.3 | 3.91 | 3.02 |
| 40 | 12.0 | 4.56 | 4.10 |
| **70** | 9.95 | **5.00** | **4.95** |
| 120 | 8.31 | 5.47 | 5.92 |
| 208 | 6.92 | 6.00 | 7.11 |

Như vậy việc khớp 0.8% là do bản thảo chọn $A\sim70$. Đây không phải một quy luật.

## Câu hỏi còn mở
- Chỉ số "bước Φ tăng" trong script (100% ở 1 MeV) khác định nghĩa "θ̇ > 0 trong 91.3% thời gian" của bản thảo. Script lấy mẫu Φ mỗi 0.05, nên dao động nhiệt ngắn hơn bị làm trơn. Hai con số không mâu thuẫn nhau.
