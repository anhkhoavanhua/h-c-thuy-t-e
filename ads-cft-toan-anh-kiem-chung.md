# Nguyên lý toàn ảnh AdS/CFT cho học thuyết $\mathcal E_{\text{dual}}$: giải số

> **Cơ sở:** bản thảo [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf), Bảng I: $\hbar\omega_0=10$ MeV, $A\sim70$, $E_{\text{vt}}=30$ MeV, $\ell_0=2.037$ fm, $\kappa=\ell_0/3$, $T_0=4.136\times10^{-22}$ s.
> Số liệu tính bằng [giai_so_ads_cft.py](giai_so_ads_cft.py). Lệnh `python3 giai_so_ads_cft.py` chạy ≈ 8 s. Script không dùng module cũ nào.

Bản thảo (7) **chưa có** mục AdS/CFT. Vì vậy mình dựng ánh xạ từ chính tham số của bản thảo, rồi kiểm từng mắt xích bằng số.

## 0. Cầu nối: AdS là một bẫy điều hòa

Với AdS toàn cục, $g_{00}=-(1+r^2/L^2)$. Giới hạn phi tương đối tính cho thế Newton $\Phi=c^2r^2/2L^2$, tức một dao động tử có

$$\omega_0=\frac{c}{L_{\text{AdS}}}\quad\Longrightarrow\quad L_{\text{AdS}}=\frac{\hbar c}{\hbar\omega_0}.$$

Như vậy hạt $x$ của $\mathcal E_{\text{dual}}$ (dao động tử $\hbar\omega_0$) chính là một nucleon khối lượng $m_N$ chuyển động trong AdS. Đây là mắt xích **duy nhất** được cố định hoàn toàn, không có tham số tự do. Mọi thứ bên dưới suy ra từ nó.

## A. Từ điển thang

| Đại lượng | Giá trị | Nhận xét |
|---|---|---|
| $L_{\text{AdS}}=\hbar c/\hbar\omega_0$ | **19.733 fm** | |
| $2\pi L/c$ | 4.1357·10⁻²² s | **Trùng** $T_0$ của Bảng I |
| $\ell_0=\sqrt{\hbar/m\omega_0}$ | 2.0364 fm | $=\sqrt{ƛ_C\,L}$ (trung bình nhân của bước sóng Compton 0.210 fm và $L$), khớp tới 4 chữ số |
| $\ell_0/L$ | 0.103 | Vận tốc $v/c$ của hạt ở biên độ $\ell_0$. Phi tương đối tính, sai số ~1% |
| $E_{\text{vt}}$ | $=3.000\,\hbar c/L$ | Điều kiện tự đồng dạng $E_{\text{vt}}/\Delta E=3$ trở thành $E_{\text{vt}}=3\hbar c/L_{\text{AdS}}$ |
| Bán kính hạt nhân $1.2A^{1/3}$ | 4.95 fm $=0.25\,L$ | Hạt nhân nằm sâu trong vùng "gần tâm" của AdS |

**Sửa giả thuyết cũ** (trong [thang-vu-tru-va-hoi-phuc-kiem-chung.md](thang-vu-tru-va-hoi-phuc-kiem-chung.md)): nếu đặt $E_{\text{vt}}=\hbar c/L_{\text{AdS}}$ thì $L_{\text{AdS}}=6.58$ fm. Giá trị này **nhỏ đúng 3 lần** so với $L$ mà bẫy điều hòa đòi hỏi. Không thể giữ cả hai. Chỉ $L=\hbar c/\hbar\omega_0$ tái lập được $T_0$ và $\ell_0$, nên cần chọn nó.

## B. Trắc địa trong AdS: đẳng thời chính xác

Tích phân RK4 các trắc địa kiểu thời gian xuyên tâm trong AdS toàn cục:

| $r_{\max}/L$ | $r_{\max}$ (fm) | $v_{\max}/c$ | Chu kỳ thời gian toàn cục / $2\pi L/c$ | Chu kỳ thời gian riêng / $2\pi L/c$ |
|---|---|---|---|---|
| 0.034 ($=\kappa$) | 0.68 | 0.034 | **1.00000000** | 1.000000 |
| 0.103 ($=\ell_0$) | 2.04 | 0.103 | **1.00000000** | 1.000000 |
| 0.25 | 4.93 | 0.243 | 1.00000000 | 1.000000 |
| 1 | 19.7 | 0.707 | 1.00000000 | 0.999999 |
| 10 | 197 | **0.995** | 1.00000000 | 0.999998 |

Chu kỳ **không phụ thuộc biên độ**, kể cả khi $v\to c$. Trong AdS, mọi trắc địa kiểu thời gian xuất phát từ một điểm đều hội tụ lại sau $\pi L/c$. Tính đẳng thời của dao động tử hạt nhân ($T_0$ không đổi) vì vậy là một tính chất **hình học chính xác** của AdS, không chỉ là xấp xỉ bậc hai.

## C. Từ điển AdS/CFT: $\Delta(\Delta-d)=m^2L^2$

Giải số phương trình Klein–Gordon toàn cục cho trường khối lượng $m_N$ ($m_NLc/\hbar=93.892$) bằng sai phân hữu hạn, có ngoại suy Richardson.

**AdS₂/CFT₁** (hạt $x$ một chiều, đúng như bản thảo): $\Delta=\tfrac12+\sqrt{\tfrac14+m^2L^2}=94.393231$.

| $n$ | $\omega$ số | $\Delta+n$ | Sai | $E-m_Nc^2$ (MeV) | Dao động tử $\hbar\omega_0(n+\tfrac12)$ |
|---|---|---|---|---|---|
| 0 | 94.393231 | 94.393231 | 2·10⁻¹⁰ | 5.0133 | 5.0000 |
| 1 | 95.393231 | 95.393231 | 9·10⁻¹⁰ | 15.0133 | 15.0000 |
| 2 | 96.393231 | 96.393231 | 2·10⁻⁹ | 25.0133 | 25.0000 |
| 4 | 98.393231 | 98.393231 | 8·10⁻⁹ | 45.0133 | 45.0000 |

**AdS₄/CFT₃** (hạt nhân ba chiều): $\Delta=\tfrac32+\sqrt{\tfrac94+m^2L^2}=95.403881$, $\omega=\Delta+2n+l$ với sai số ≤ 4·10⁻⁹.

| Vỏ $N=2n+l$ | $(n,l)$ | Bậc suy biến AdS₄ | Dao động tử 3D | $E-m_Nc^2$ (MeV) |
|---|---|---|---|---|
| 0 | (0,0) | 1 | 1 | 15.12 |
| 1 | (0,1) | 3 | 3 | 25.12 |
| 2 | (1,0) (0,2) | 6 | 6 | 35.12 |
| 3 | (1,1) (0,3) | 10 | 10 | 45.12 |
| 4 | (2,0) (1,2) (0,4) | 15 | 15 | 55.12 |

**Kết quả:**
1. Phổ AdS **trùng** phổ dao động tử của bản thảo, chỉ lệch một dịch tương đối tính **13 keV** (AdS₂) hoặc **120 keV** (AdS₄), tức ≈ $d^2(\hbar\omega_0)^2/8m_Nc^2$.
2. Bậc suy biến vỏ khớp hoàn toàn. Số nucleon lấp đầy vỏ (×2 spin) là 2, 8, 20, 40, **70**, nên $A\sim70$ của Bảng I đúng là vỏ $N=4$ đóng.
3. Về phía biên (CFT), nucleon ứng với một toán tử có chiều dị thường $\Delta\approx94.4$ (CFT₁) hoặc 95.4 (CFT₃). Mỗi mức kích thích $\hbar\omega_0$ của hạt là một "hậu duệ" $\Delta+n$ của toán tử đó.

## D. Thời gian trong AdS tuần hoàn: nối với đồng hồ Page–Wootters

| Kiểm tra | Kết quả |
|---|---|
| Pha $e^{-i\Delta\tau}$ sau $\tau=2\pi$ | $\Delta \bmod 1=0.393$, **không nguyên** |
| Bó sóng ~7 mức: $\lvert\langle\psi(0)\vert\psi(\pi)\rangle\rvert^2$ | 0.000000 (hạt ở điểm đối cực) |
| $\lvert\langle\psi(0)\vert\psi(2\pi)\rangle\rvert^2$, $\lvert\langle\psi(0)\vert\psi(4\pi)\rangle\rvert^2$ | **1.000000**, 1.000000 |

1. Trên hyperboloid AdS gốc, thời gian đóng thành vòng ($\tau\sim\tau+2\pi$). Khi $\Delta$ không nguyên, trường của nucleon **không xác định** trên đó. Phải dùng **phủ phổ quát**, nơi $\tau\in\mathbb R$.
2. Dù vậy, trạng thái vật lý của hạt vẫn lặp lại chính xác sau mỗi $2\pi L/c$. Một hạt trong AdS là đồng hồ tuần hoàn, **giống hệt kết luận OAM** ở [pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md): không tự đếm được vòng.
3. Hình học AdS đòi thời gian trên $\mathbb R$ (phủ phổ quát). Mục V.C.4 của bản thảo đòi pha tích lũy $\Phi\in\mathbb R$. Hai đòi hỏi **trùng nhau**. Bộ đếm vòng (chuỗi xung) chính là thứ phân biệt các tấm của phủ phổ quát.

## E. Cận Bekenstein cho đợt Carnot

$S\le2\pi ER/\hbar c$ với $R=4.95$ fm, so với số trạng thái exciton Ericson ($g=A/13$).

| $E$ (MeV) | Cận, chỉ năng lượng kích thích | $\ln$(số exciton trong 1 MeV) | Cận, kể cả khối lượng nghỉ |
|---|---|---|---|
| 8.0 | 1.26 nat | 13.2 nat | 10 351 nat |
| 13.0 | 2.05 nat | 17.7 nat | 10 352 nat |
| 20.4 (trần đợt) | 3.21 nat | 22.2 nat | 10 353 nat |

- Cận Bekenstein được thỏa với dư ~500 lần, **nhưng chỉ khi tính cả khối lượng nghỉ** (đúng định nghĩa của cận).
- Nếu chỉ tính năng lượng kích thích Carnot thì bị vi phạm 7–10 lần. Vì vậy **không** được viết "entropy của đợt bị chặn bởi diện tích/năng lượng đợt".
- Hạt nhân còn cách xa việc bão hòa cận toàn ảnh. Nó không phải một hệ "kiểu lỗ đen".

## F. Những gì toàn ảnh **không** cho

1. **AdS/CFT đầy đủ cần một CFT biên có N lớn, liên kết mạnh.** Kết quả C chỉ là phần **động học** của từ điển: nhóm đối xứng SO(d,2) cùng liên hệ $\Delta\leftrightarrow m$. Điều này đúng cho bất kỳ trường tự do nào trong AdS. Một hạt cộng một rô-to không phải là CFT có đối ngẫu hấp dẫn.
2. **Rô-to và ràng buộc phá đối xứng bảo giác.** Các số hạng $\Delta E\cos\theta$, $gx\sin\theta$, ràng buộc $x=\kappa\sin\theta$ và va chạm Carnot đều không thuộc SO(d,2). Pha khóa ép hạt chạy ở $\dot\theta_d\approx2.6\,\omega_0$, **khác** tần số riêng $c/L$ của AdS. Trong bức tranh toàn ảnh, $\mathcal E_{\text{dual}}$ là một **nhiễu loạn cưỡng bức** lên nền AdS, không phải một phần của hình học.
3. **Ryu–Takayanagi** (entropy vướng víu bằng diện tích mặt cực tiểu) cần một lý thuyết trường ở biên. Mô hình 0+1 chiều không có vùng con để chia, nên không kiểm được. Entropy rô-to–hạt $S_{vN}\approx0.08$ ở bài trước là vướng víu giữa hai bậc tự do, không phải entropy RT.

## G. Đoạn văn có thể thêm vào bản thảo (tiếng Anh)

> *Holographic scale identification.* In the non-relativistic limit, global AdS acts as an isotropic harmonic trap with $\omega_0=c/L_{\text{AdS}}$. Matching $\hbar\omega_0=10$ MeV fixes $L_{\text{AdS}}=19.73$ fm, reproducing $T_0=2\pi L_{\text{AdS}}/c$ and $\ell_0=\sqrt{ƛ_C L_{\text{AdS}}}$, with $E_{\text{vt}}=3\hbar c/L_{\text{AdS}}$. The global Klein–Gordon spectrum of a nucleon-mass field, $\omega=\Delta+2n+l$ with $\Delta(\Delta-3)=m_N^2L_{\text{AdS}}^2$, reproduces the shell-model oscillator levels and degeneracies up to a relativistic shift of 0.12 MeV. The rotor coupling and Carnot impacts explicitly break the SO(3,2) symmetry, so this identification is kinematic and does not imply a holographic dual of the engine.

## H. Lỗi còn trong bản thảo (7)

1. Có **5 chỗ "[cite: 14]"** sót lại, ở V.B.4, V.C.2 (2 chỗ) và V.C.4 (2 chỗ). Danh mục tài liệu chỉ có [1]–[10], nên phải xóa các chỗ này hoặc thay bằng trích dẫn đúng.
2. **V.C.4** ghi "100% temporal fidelity". Con số này chỉ đúng khi xếp thứ tự theo **đợt**. Khoảng 36% cặp gãy→khóa kề nhau cách nhau dưới 0.01 $T_0$ và không chuỗi xung nào tách được. Nên thêm "at the burst level".

## Câu hỏi còn mở
- Muốn đưa rô-to vào hình học (ví dụ θ là góc của AdS₃ với biên $S^1$), cần một phổ tuyến tính theo $\ell$. Rô-to của bản thảo có phổ $\ell^2/2I$ (bậc hai) nên chưa khớp. Có muốn thử hướng này không?
