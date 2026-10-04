# Part III (LaTeX): số liệu, cách xử lý ba giới hạn, nguồn từng con số

> **Bản thảo:** [paper_part3.tex](paper_part3.tex) → [paper_part3.pdf](paper_part3.pdf). Định dạng revtex4-2, cùng kiểu Part I. Biên dịch bằng `tectonic paper_part3.tex`: sạch, không còn overfull box.
> **Số liệu mới:** [giai_so_part3.py](giai_so_part3.py), gồm `N1` ≈ 4 phút, `N2` ≈ 1 phút, `N3` ≈ 2 phút.
> **Cơ sở:** [HocThuyetE (7).pdf](../../HocThuyetE%20(7).pdf) (Part I).

## 1. Ba giới hạn được nêu và cách Part III xử lý

| Giới hạn | Cách xử lý trong Part III | Số liệu |
|---|---|---|
| **WdW chỉ chạy được 7 tích tắc** (môi trường $3^k$) | **Gỡ chính xác, không xấp xỉ.** Mọi đại lượng người quan sát trong đo được chỉ phụ thuộc $\rho_S(k)=\Lambda^k(\rho_0)$ (ánh xạ Kraus rút gọn) và trọng số Born của nhánh bản ghi. Cả hai tính được trong chiều của S | Khớp trạng thái lịch sử **đầy đủ 260 253 chiều** tới **1.4·10⁻¹⁵**. Kéo dài tới **60 $T_0$ = 1500 tích tắc**. 200 nhánh Born: ⟨ℓ⟩ = 2.19 ± 0.07 so với 2.262 |
| **Chưa giải WdW liên tục** | Đồng hồ lý tưởng $\hat p_c$ cho nghiệm giải tích $\Psi(c)=e^{-i\hat H_{SE}c}\Psi(0)$. Bản rời rạc hội tụ về Lindblad khi $\tau\to0$ | Sai số giảm **đúng một nửa** mỗi lần chia đôi τ (tỉ số 1.95–2.04), tức đúng $O(\tau)$ |
| **AdS/CFT không có CFT N lớn** | **Không tuyên bố toàn ảnh.** Trình bày là phép nhúng **động học**, ghi rõ phạm vi, kèm phép thử bác bỏ trùng hợp Hawking–Page | L = 19.733 fm; đẳng thời 10⁻⁸; phổ KG ≤ 8·10⁻⁹; θ̇_br lệch −0.08%. Trùng hợp $T_{\text{HP}}$ bị bác bỏ khi đổi $\gamma_x$, $\Delta E$ |

**Vì sao cách xử lý này chịu được phản biện:**
- Không có điều gì cần "chứng minh đối ngẫu toàn ảnh trọn vẹn", vì bài **không tuyên bố** đối ngẫu.
- WdW được giải **chính xác** cho mọi đại lượng quan sát được. Việc giới hạn 7 tích tắc chỉ dùng để kiểm chéo.

## 2. Số liệu mới chạy cho Part III

### N1. Trạng thái lịch sử, ánh xạ rút gọn, giới hạn liên tục
- $\lVert\sum_jK_j^\dagger K_j-1\rVert=3\cdot10^{-14}$.
- $k\le7$: $\max\lVert\operatorname{Tr}_E\langle k\vert\Psi\rangle\rangle-\Lambda^k(\rho_0)\rVert_1/2=$ **1.4·10⁻¹⁵**, trên S 119 chiều với $E=3^7$.

**Hội tụ về Lindblad (khoảng cách vết):**

| τ | t = 1 $T_0$ | 2 $T_0$ | 4 $T_0$ |
|---|---|---|---|
| $T_0/25$ | 4.06·10⁻² | 3.33·10⁻² | 2.03·10⁻² |
| $T_0/50$ | 1.99·10⁻² | 1.64·10⁻² | 1.04·10⁻² |
| $T_0/100$ | 9.90·10⁻³ | 8.16·10⁻³ | 5.26·10⁻³ |
| $T_0/200$ | 4.94·10⁻³ | 4.08·10⁻³ | 2.65·10⁻³ |

**Dài hạn** ($\tau=T_0/25$):

| $t/T_0$ | ⟨ℓ⟩ | ⟨H⟩ | Tr ρ² | S(E) |
|---|---|---|---|---|
| 0 | 3.000 | 4.970 | 1.000 | 0 |
| 1 | 2.464 | 4.683 | 0.377 | 1.250 |
| 20 | 2.263 | 4.122 | 0.241 | 1.615 |
| 60 | 2.262 | 4.120 | 0.241 | 1.615 |

Entropy bản ghi tăng đơn điệu rồi bão hòa khi động cơ vào chu trình dừng.

### N2. Bản ghi môi trường là nhánh của |Ψ⟩⟩
200 nhánh lấy mẫu theo Born, mỗi nhánh 60 $T_0$:
- Mỗi nhánh **thuần** (lệch < 10⁻¹⁵).
- Mỗi bản ghi Carnot: trung bình **1.93 $\hbar\omega_0$**, trung vị 2.17 (V.C.1: 2.05).
- Đợt bản ghi: **1.51** bản ghi/đợt (V.C.1: 1.57). Thời gian chờ giữa các đợt trung bình 2.84 $T_0$, CV 0.70.

### N3. Đồng hồ rô-to có tác dụng J

| J (ħ) | 1 − F (một vòng) | 1 − F (WdW) |
|---|---|---|
| 3 | 0.054 | 0.307 |
| 6 | 0.026 | 0.101 |
| 12 | 0.020 | 0.046 |
| 30 | 0.015 | 0.048 |
| 60 | 0.014 | 0.050 |

**Phát hiện: có sàn sai số.** Phản lực lên đồng hồ chỉ trội khi J ≲ 12ħ. Trên mức đó 1 − F **bão hòa**, và sàn này **không** do độ phân giải đồng hồ:
- Ở J = 30, tăng σ_ℓ từ 3.2 lên 10 làm 1 − F **xấu đi**, từ 0.015 lên 0.026.
- Sàn ở cửa sổ một vòng khớp với rò rỉ từ các vòng lân cận, $e^{-1/(2\cdot0.35^2)}=0.017$.
- Sàn ở giới hạn WdW (0.05) là do trộn kết hợp $t$ với $t+nT_{\text{rot}}$.

Kết luận: đồng hồ tuần hoàn **không bao giờ** tái lập Schrödinger chính xác. Đây là lý lẽ định lượng cho mục V.C.4 (cần bộ đếm). Ban đầu mình thử khớp luật lũy thừa $J^p$, nhưng các điểm bão hòa làm phép khớp vô nghĩa, nên mình **không** đưa nó vào bài.

## 3. Nguồn các số liệu khác trong bài

| Mục trong bài | Tài liệu nguồn |
|---|---|
| II.A Phổ trên đa tạp ràng buộc, nhòe pha, hồi sinh | [schrodinger-rang-buoc-kiem-chung.md](schrodinger-rang-buoc-kiem-chung.md) |
| II.B Phân rã bước nhảy (92% rô-to, 99.7% năng lượng) | [mcwf-tach-buoc-kiem-chung.md](mcwf-tach-buoc-kiem-chung.md) |
| III.A Ràng buộc WdW 4·10⁻¹⁵ (S 36 chiều) | [schrodinger-pw-wdw-e7-tong-hop.md](schrodinger-pw-wdw-e7-tong-hop.md) |
| III.E OAM, chuỗi xung, 36% rung gãy/khóa | [pw-luong-tu-moi-truong-dong-ho-dem-vong.md](pw-luong-tu-moi-truong-dong-ho-dem-vong.md) |
| IV.A Sieve với kênh Davies thật, tốc độ nhảy 0.023 vs 0.182 | [nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md) mục D |
| IV.B Darwinism R = 3.3, Born | [nguoi-quan-sat-kiem-chung.md](nguoi-quan-sat-kiem-chung.md) mục B |
| IV.C Englert V, D | [rung-dong-m-song-hat-kiem-chung.md](rung-dong-m-song-hat-kiem-chung.md) |
| V AdS, Hawking–Page | [ads-cft-toan-anh-kiem-chung.md](ads-cft-toan-anh-kiem-chung.md), [ads-cft-mo-rong-kiem-chung.md](ads-cft-mo-rong-kiem-chung.md) |

## 4. Những gì bài **cố ý không** tuyên bố
- **Không** có đối ngẫu toàn ảnh của động cơ. Mục V tự nêu điều này trong đoạn "Scope".
- **Không** giải thích vì sao ra đúng một kết quả cụ thể. Mục IV.B ghi rõ đây là vấn đề diễn giải.
- **Không** gắn giới hạn $k_BT\lesssim1$ MeV với Hawking–Page, vì phép thử cho thấy đó là trùng hợp.
- Mô hình va chạm là Markov, sai số $O(\tau)$. HEOM phi Markov để cho các phần sau (mục Limitations).

## Câu hỏi còn mở
- Bài trích Part I là "manuscript (2026)". Nếu Part I đã nộp hoặc có arXiv, cần thay bằng số tham chiếu thật.
- Bài nhắc "Part III" nhưng không trích Part II. Nếu có Part II, cần bổ sung vào phần giới thiệu.
- Email lấy theo PDF (7) (khoa181101@gmail.com). Bản `paper_prc.tex` cũ dùng một email khác. Bạn cần chọn một email thống nhất.
