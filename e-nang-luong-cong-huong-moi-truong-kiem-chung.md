# E năng lượng bên ngoài có cấu trúc: cộng hưởng quyết định thời gian có chạy hay không

> **Nối tiếp:** [e-nang-luong-ben-ngoai-pho-kiem-chung.md](e-nang-luong-ben-ngoai-pho-kiem-chung.md). Ở bản đó môi trường là bể phẳng, và động cơ nhả năng lượng ra theo 23 vạch rô-to.
> **Câu hỏi:** nếu E năng lượng bên ngoài **có cấu trúc riêng**, tức chỉ nhận năng lượng quanh một vài tần số, thì điều gì xảy ra với thời gian (h), dòng bản ghi và năng lượng của động cơ?
> **Số liệu:** [giai_so_e_nang_luong_cong_huong.py](giai_so_e_nang_luong_cong_huong.py), kết quả thô [e_nang_luong_cong_huong_ket_qua.txt](e_nang_luong_cong_huong_ket_qua.txt), chạy khoảng 9 giây (phần `extra` khoảng 13 giây).

## Mô hình

Dùng Part III (S 36 chiều, bơm a = 1). Kênh Carnot được nhân thêm một mật độ phổ có đỉnh:

$$L_C(\omega_p)=\sqrt{0.5}\;V\big[\mathrm{triu}(V^\dagger uV)\odot\sqrt{f(\omega_{ij})}\big]V^\dagger,\qquad f(\omega)=\varepsilon+(1-\varepsilon)\frac{\gamma^2}{(\omega-\omega_p)^2+\gamma^2}$$

với γ = 0.15 ħω0 (độ rộng đỉnh) và ε = 0.05 (sàn 5%).
- Đỉnh cao bằng bể phẳng, nên ở đúng ω_p môi trường nhận năng lượng như bể chuẩn; xa đỉnh nó chỉ nhận 5%.
- Kiểm: khi không có đỉnh, kênh dựng lại khớp kênh gốc tới 5.6·10⁻¹⁶.

## Kết quả: quét vị trí đỉnh ω_p

| Môi trường | J_C (/T0) | P nhả (ħω0/T0) | **h (nat/T0)** | ΔE vượt chân không (ħω0) | Vạch nhả mạnh nhất |
|---|---|---|---|---|---|
| Bể phẳng (chuẩn) | 0.494 | 1.370 | 5.23 | 4.26 | 3.54, 2.47, 2.45 |
| ω_p = 0.50 | 0.058 | 0.186 | 0.73 | 8.17 | 3.54 (qua sàn) |
| ω_p = 1.00 | 0.057 | 0.190 | 0.78 | 8.20 | 3.54 |
| ω_p = 1.21 | 0.059 | 0.194 | 0.77 | 8.19 | 3.54 |
| ω_p = 1.75 | 0.062 | 0.206 | 0.76 | 8.16 | 3.54 |
| ω_p = 2.00 | 0.067 | 0.220 | 0.83 | 8.14 | 3.54 |
| **ω_p = 2.46** (vạch rô-to) | **0.126** | 0.367 | **1.46** | 7.72 | 3.54, **2.47, 2.45** |
| ω_p = 3.00 | 0.115 | 0.387 | 1.27 | 7.83 | 3.54 |
| **ω_p = 3.54** (vạch rô-to chính) | **0.375** | **1.295** | **3.46** | **6.06** | **3.54 (1.21)** |
| ω_p = 4.20 | 0.094 | 0.323 | 1.19 | 7.94 | 3.54 |
| ω_p = 5.00 | 0.063 | 0.213 | 0.80 | 8.17 | 3.54 |
| ω_p = 6.00 | 0.064 | 0.237 | 0.89 | 8.11 | 3.54, 5.99, 6.01 |
| **Hai đỉnh 2.46 + 3.54** | 0.369 | 1.154 | **3.64** | 6.08 | 3.54, 2.47, 2.45 |

## Đọc kết quả

1. **Cộng hưởng là có thật và rất sắc.** Khi đỉnh môi trường trùng vạch rô-to chính 3.54 ħω0:
   - dòng bản ghi tăng **6.5 lần** (0.375 so với khoảng 0.06 khi lệch);
   - độ mở tương lai tăng **4.5 lần** (3.46 so với khoảng 0.75).

   Vạch thứ hai 2.46 ħω0 cho cộng hưởng yếu hơn (h = 1.46). Ở các vị trí khác chỉ còn mức nền do sàn 5%.
2. **Môi trường lệch cộng hưởng làm thời gian gần đóng băng.** h giảm còn 0.73–0.89, khoảng 15% so với bể phẳng. Động cơ **không xả được năng lượng** nên năng lượng dồn lại: ΔE ≈ 8.2 ħω0, gần **gấp đôi** bể phẳng (4.26). Hệ bị kẹt ở trên, nghĩa là Hiện tại thắng: đúng chế độ "một pha thắng ⟹ đóng băng" của [xung-dot-qua-khu-hien-tai](xung-dot-qua-khu-hien-tai-mo-bang-thoi-gian-kiem-chung.md), nhưng lần này do **cấu trúc của E năng lượng bên ngoài** gây ra.
3. **Môi trường hai đỉnh** (đặt đúng cả hai vạch rô-to) cho h = 3.64, nhỉnh hơn một đỉnh ở 3.54. Khi môi trường có "nhiều năng lượng khác nhau" **trùng với bậc thang của động cơ**, thời gian chạy tốt hơn.
4. **Bể phẳng vẫn cho h lớn nhất (5.23)**, vì nó nhận mọi tần số với cùng cường độ, trong khi đỉnh hẹp γ = 0.15 cắt bớt phần lớn. Thứ quyết định không phải "có đỉnh hay không" mà là **bao nhiêu phổ của môi trường phủ lên các vạch của động cơ**.

## Kết luận

E năng lượng bên ngoài **chọn lọc** thời gian:
- Động cơ chỉ thoát khỏi 0 và giữ được thời gian chạy khi môi trường bên ngoài có năng lượng **tại đúng các bậc thang rô-to** của nó (2.46 và 3.54 ħω0, tức 25 và 35 MeV).
- Nếu môi trường "nói một ngôn ngữ năng lượng khác" (đỉnh lệch), động cơ bị nhốt với năng lượng dồn ứ, bản ghi gần như ngừng, và thời gian gần đóng băng.
- Cộng hưởng giữa vũ trụ bên trong và E năng lượng bên ngoài là điều kiện để thời gian chảy.

**Câu tiếng Anh cho bản thảo:**

> The outer environment selects whether time flows. Giving the Carnot bath a peaked spectral density (width 0.15 ħω₀) shows a sharp resonance: with the peak on the main rotor line (3.54 ħω₀) the record current is 6.5 times and the future-opening rate h 4.5 times larger than off resonance, with a secondary resonance at 2.46 ħω₀. Off resonance the engine cannot discharge its energy, which piles up to ≈ 8.2 ħω₀ above the vacuum (twice the flat-bath value) while h drops to ≈ 0.75 nat per T₀: the Present wins and time nearly freezes. Time keeps running only when the energy outside matches the rotor ladder inside.

## Giải tiếp 1: độ rộng đỉnh γ

Chạy `python3 giai_so_e_nang_luong_cong_huong.py extra`, khoảng 13 giây. Bể Carnot có đỉnh, bơm phẳng.

| γ (ħω0) | Một đỉnh 3.54: J_C | h | Hai đỉnh 2.46 + 3.54: J_C | h |
|---|---|---|---|---|
| 0.05 | 0.372 | 3.38 | 0.364 | 3.54 |
| 0.15 | 0.375 | 3.49 | 0.369 | 3.58 |
| 0.30 | 0.384 | 3.63 | 0.380 | 3.78 |
| 0.60 | 0.413 | 3.95 | 0.410 | 4.08 |
| 1.00 | 0.445 | 4.36 | 0.443 | 4.46 |
| 2.00 | 0.479 | 4.84 | 0.478 | 4.92 |
| 4.00 | 0.491 | 5.10 | 0.490 | 5.15 |
| *Bể phẳng* | *0.494* | *5.23* | | |

- h tăng **đơn điệu** theo γ và tiến dần tới bể phẳng. Với chiều cao đỉnh cố định bằng bể phẳng, môi trường có đỉnh **không bao giờ vượt** bể phẳng: phủ rộng hơn thì tốt hơn.
- Ngay cả đỉnh rất hẹp (γ = 0.05) vẫn giữ được 65–68% h của bể phẳng. Chỉ cần đúng vạch rô-to là đủ cho phần lớn dòng thời gian.
- Hai đỉnh luôn nhỉnh hơn một đỉnh (+2–5%).
- **Lưu ý:** phép so này giữ **chiều cao** đỉnh bằng nhau. Nếu giữ **tổng cường độ ghép** bằng nhau (đỉnh hẹp thì cao hơn), kết luận có thể khác. Chưa kiểm.

## Giải tiếp 2: nguồn bơm cũng có cấu trúc (γ = 0.15)

| Bơm | Bể Carnot | J_C (/T0) | P vào (ħω0/T0) | **h (nat/T0)** | ΔE (ħω0) |
|---|---|---|---|---|---|
| phẳng | phẳng | 0.494 | 1.370 | **5.22** | 4.26 |
| đỉnh 2.46 + 3.54 | phẳng | 0.080 | 0.190 | 1.05 | 1.30 |
| phẳng | đỉnh 2.46 + 3.54 | 0.369 | 1.154 | 3.56 | 6.08 |
| **đỉnh 2.46 + 3.54** | **đỉnh 2.46 + 3.54** | 0.171 | 0.501 | **2.21** | 3.46 |
| đỉnh 1.0 (lệch) | đỉnh 2.46 + 3.54 | 0.055 | 0.104 | 1.01 | 2.41 |
| đỉnh 2.46 + 3.54 | đỉnh 1.0 (lệch) | 0.052 | 0.167 | **0.66** | **7.35** |

- **Giả thuyết "cả hai cùng cộng hưởng thì nhanh hơn bể phẳng" bị bác bỏ.** Hai phía cùng cộng hưởng cho h = 2.21, chỉ bằng 42% so với cả hai phẳng (5.22).
- **Bơm có cấu trúc gây hại nhiều hơn bể có cấu trúc:** h = 1.05 so với 3.56. Bơm phải đẩy rô-to lên qua **mọi** bậc ℓ → ℓ+1 trên đường tới ngưỡng gãy. Các bậc đó có năng lượng khác 2.46 và 3.54, nên một nguồn chỉ cấp ở hai tần số sẽ chặn đường lên.
- **Cấp vào đúng nhưng xả ra lệch** cho tình trạng tệ nhất: h = 0.66 và năng lượng dồn **7.35 ħω0**. Động cơ nạp được mà không xả được, nên Hiện tại thắng và thời gian gần đóng băng.
- **Cấp vào lệch, xả ra đúng:** h = 1.01, năng lượng thấp (2.41). Động cơ không nạp đủ, nên Quá khứ thắng.
- **Cảnh báo:** kênh bơm tạo kết hợp lượng tử giữa các trạng thái riêng. Lọc nó trong cơ sở riêng (như ở đây) làm đổi bản chất của kênh, không chỉ đổi phổ. Kết quả phần này là định tính.

## Kết luận cập nhật

1. E năng lượng bên ngoài phải có năng lượng **tại các vạch rô-to** ở phía **nhận** (Carnot). Lệch cộng hưởng thì thời gian gần đóng băng.
2. Ở phía **cấp** (bơm), môi trường cần **phổ rộng** để đẩy rô-to qua mọi bậc lên ngưỡng. Một nguồn hẹp thì chặn đường lên.
3. Trong mọi cấu hình đã thử, **môi trường phổ rộng ở cả hai phía cho thời gian chạy nhanh nhất** (h = 5.22). Mỗi lần thu hẹp hay làm lệch phổ đều kéo hệ về phía một pha thắng:
   - nạp không xả ⟹ Hiện tại thắng;
   - xả không nạp ⟹ Quá khứ thắng.

   Đây là lời giải frozen time một lần nữa, nhưng lần này nhìn từ phía E năng lượng bên ngoài.

## Giới hạn

- Mật độ phổ được đưa vào theo kiểu Davies (nhân biên độ chuyển dời với √f(ω)), chưa phải một mô hình bể vi mô có bộ nhớ. Với đỉnh hẹp, các hiệu ứng phi Markov có thể quan trọng.
- Kết quả phụ thuộc γ và sàn ε. Với ε = 0, có thể xuất hiện trạng thái tối và h giảm thêm. Chưa quét hai tham số này.
- Kênh bơm giữ nguyên phẳng. Chỉ môi trường **nhận** năng lượng (Carnot) có cấu trúc.

## Câu hỏi còn mở

- Đã giải: quét γ (đỉnh không vượt bể phẳng khi giữ cùng chiều cao) và bơm có cấu trúc (cả hai cùng cộng hưởng vẫn thua bể phẳng: 2.21 so với 5.22).
- Còn mở: so sánh khi giữ **tổng cường độ ghép** bằng nhau thay vì chiều cao đỉnh; và dùng bể có bộ nhớ (phi Markov, HEOM) cho các đỉnh hẹp.
