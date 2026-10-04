# Nguồn lái tương tác mạnh, cơ chế Kramers và $\gamma_x=\text{const}$: giải và kiểm chứng số liệu

> Tiếp nối [scaling-A-kiem-chung.md](scaling-A-kiem-chung.md).
> Số liệu tính bằng [giai_so_kramers_A.py](giai_so_kramers_A.py) (chạy `python3 giai_so_kramers_A.py`, ≈ 25 phút trên 8 lõi).
> **Bộ thông số chuẩn:** $I=1,\ m=k=1,\ \chi=1,\ \kappa=1/3,\ R=3,\ g=0.1,\ e=0,\ \tau=10,\ \gamma_x=0.70,\ \Omega_0=1.05\,\dot\theta_{\text{gãy}}$ ($\dot\theta_{\text{gãy}}=3.2048$).

## Tóm tắt kết luận

| Khẳng định | Kết quả kiểm chứng |
|---|---|
| $\alpha_{\text{eff}}=\alpha_{\pi NN}\approx14$ ⟹ $\Gamma/\omega_0\approx0.10$, $\tau\approx10$ | **Số học đúng** ($\Gamma/\omega_0=0.0995$, $\tau=10.05$) |
| Cản/lái do "bức xạ" meson ở 10 MeV | **Không thể về động học.** Muốn phát một pion thật cần $m_\pi c^2=139.6$ MeV, gấp 14 lần $\hbar\omega_0=10$ MeV. Công thức kiểu Larmor chỉ áp cho lượng tử không khối lượng. Trao đổi pion ảo không tạo ra tiêu tán bức xạ |
| Giá trị 14 | **Phụ thuộc cách viết hằng số ghép.** Dạng giả-vô hướng: $g^2/4\pi\approx13.6$–$14$. Dạng giả-vectơ tương đương: $f^2/4\pi\approx0.075$, cho **$\tau=1877$**. Hai dạng liên hệ qua $g^2/4\pi=(2m_p/m_\pi)^2f^2/4\pi=13.56$. Chọn dạng nào để ra $\tau=10$ là một lựa chọn, không phải suy ra |
| Gãy là vượt rào Kramers, $\Delta K=1.24$ | **Bị bác bỏ bằng số.** Tốc độ gãy dừng (đã bỏ giai đoạn khởi động) **gần như không đổi** khi $k_BT$ tăng 8 lần: 0.158 → 0.158 → 0.143 → 0.140 → 0.146 mỗi $T_0$. Khớp Arrhenius $\ln r=a+b/k_BT$ cho **$b=+0.013$**, trong khi Kramers dự đoán **$b=-\Delta K=-1.675$**. Kramers dự đoán $r(1\text{ MeV})/r(5\text{ MeV})=1.5\times10^{-6}$; **đo được 1.12** |
| "$N\approx15$–$20$ lần gãy khớp tuyệt đối Monte Carlo (10.6–22.9)" | **Không phải phép tính.** Kramers cần một prefactor $\omega_a\omega_b/(2\pi\gamma_\theta)$ nhưng $\omega_a,\omega_b$ không được xác định. Con số 15–20 được **gán**, không tính ra. Ngoài ra 10.6–22.9 lần là số đếm **có cả giai đoạn khởi động** (xuất phát trên ngưỡng) |
| $\lambda_{\text{mfp}}\propto R_{\text{hạt nhân}}\propto A^{1/3}$ ⟹ $\gamma_x=\text{const}$ | **Không có căn cứ.** Quãng đường tự do trung bình của nucleon trong chất hạt nhân được quyết định bởi mật độ và năng lượng ($\sim1/\rho\sigma$), gần như không đổi theo $A$. Đây là một giả thiết thêm vào |
| Với $\gamma_x=\text{const}$: "$\eta\approx0.40$ cho mọi hạt nhân" | **Sai.** Đo với $\gamma_x=0.70$, tỉ số 3: **$\eta=0.478$ (O), 0.418 (Ca), 0.422 (Pb)**. $\eta$ vẫn thay đổi theo $A$, qua nhiễu $k_BT/\hbar\omega_0(A)$ |
| Bảng: O 6.5 (+6%), Ca 4.8 (**+20%**), Pb 2.8 (+7%) | Với $\eta$ giả định 0.40: Ca là $16.4\cdot40^{-1/3}=4.795$, nên lệch **+28.3%**, không phải +20%. **Với $\eta$ đo được: O 7.78 (+27%), Ca 5.01 (+34%), Pb 2.92 (+12%)** |
| 8 tham số "cố định tiên nghiệm", gồm $E_{\text{vt}}=30$ MeV cho mọi hạt nhân | **Không tự nhất quán với luật $A^{-1/3}$.** Giữ $E_{\text{vt}}=30$ MeV cố định thì tỉ số $E_{\text{vt}}/\hbar\omega_0$ thay đổi theo $A$ (1.84 → 4.34). Mô phỏng khi đó cho **O 4.24 (−31%), Ca 4.13 (+10%), Pb 3.45 (+32%)**, tức gần như phẳng ($\approx A^{-0.08}$). Ngoài ra, độ sâu Woods–Saxon điển hình là ~50 MeV, không phải 30 MeV |
| Danh sách tham số | Còn thiếu $g=0.1$, $\Omega_0=1.05\,\dot\theta_{\text{gãy}}$ và $k_BT$, cũng là các tham số chọn tay |

## 1. τ từ tương tác mạnh

| Hằng số ghép | $\Gamma/\omega_0=\tfrac23\alpha\,\hbar\omega_0/m_pc^2$ | $\tau$ |
|---|---|---|
| $g^2/4\pi\approx14$ (giả-vô hướng) | 0.0995 | **10.05** |
| $f^2/4\pi\approx0.075$ (giả-vectơ) | 0.00053 | **1 877** |

Ngưỡng phát pion thật: $m_\pi c^2/\hbar\omega_0=14$. Dưới ngưỡng, không có kênh phát pion thật.

## 2. Kiểm Kramers: tốc độ gãy ở trạng thái dừng

**Thiết lập:**
- Xuất phát ở $\dot\theta=0.82\,\dot\theta_{\text{gãy}}$, tức dưới ngưỡng. Bỏ $10T_0$ đầu làm giai đoạn khởi động.
- Mỗi $k_BT$: đo 16 quỹ đạo × $90T_0$, tổng $1440\,T_0$.
- "Đợt thoát" là lần gãy xảy ra sau khi đã khóa ổn định ít nhất $T_0$; điều này loại các lần gãy lại ngay sau khi vừa khóa.
- Không nhiễu (tất định): **0 lần gãy** trong $50T_0$, và $\dot\theta$ dừng ở 2.32–2.62.

| $k_BT$ | Tổng số gãy | Đợt thoát | $r$ (mỗi $T_0$) |
|---|---|---|---|
| 1 MeV (0.1) | 1 146 | 227 | 0.1576 |
| 2 MeV (0.2) | 1 803 | 228 | 0.1583 |
| 3 MeV (0.3) | 1 904 | 206 | 0.1431 |
| 5 MeV (0.5) | 2 591 | 202 | 0.1403 |
| 8 MeV (0.8) | 3 082 | 210 | 0.1458 |

**Khớp Arrhenius:** độ dốc $b=+0.013$. Kramers dự đoán $b=-\Delta K=-1.675$, với $\Delta K=\tfrac12I(\dot\theta_{\text{gãy}}^2-\dot\theta_{\text{dừng}}^2)$ và $\dot\theta_{\text{dừng}}=2.631$ khi $\gamma_x=0.70$.

**Đọc kết quả:**
- Chỉ cần một chút nhiễu (1 MeV) là đã gãy, nhưng tốc độ **bão hòa** ở ≈ 0.15 lần mỗi $T_0$, tức khoảng một đợt mỗi 6.7 $T_0$, bất kể nhiệt độ.
- Vậy tốc độ gãy **không** do độ cao rào năng lượng quyết định. Nhiều khả năng nó bị giới hạn bởi thời gian nạp lại năng lượng sau mỗi lần tái khóa.
- Nhiễu mạnh hơn chỉ làm tăng số lần **gãy lại ngay sau khi khóa**: tổng số gãy tăng 1 146 → 3 082, trong khi số đợt thoát giữ nguyên.
- Hình ảnh "trạng thái giả bền, vượt rào 12.4 MeV theo $e^{-\Delta K/k_BT}$" **không đúng** trong mô hình này.

## 3. $\eta$ với $\gamma_x=0.70$ (nhiễu 1 MeV, $N=32$, $10T_0$)

| Hạt nhân | Giả thiết | $E_{\text{vt}}/\hbar\omega_0$ | $\eta$ đo | $\Delta E$ mô hình | $16.4A^{-1/3}$ (đề xuất) | 3⁻ thực nghiệm | Lệch |
|---|---|---|---|---|---|---|---|
| ¹⁶O | tỉ số 3 cố định | 3.000 | 0.4779 | 7.78 MeV | 6.51 | 6.130 | **+26.8%** |
| ¹⁶O | $E_{\text{vt}}=30$ MeV cố định | 1.844 | 0.2606 | 4.24 MeV | 6.51 | 6.130 | **−30.8%** |
| ⁴⁰Ca | tỉ số 3 cố định | 3.000 | 0.4175 | 5.01 MeV | 4.80 | 3.737 | **+33.9%** |
| ⁴⁰Ca | $E_{\text{vt}}=30$ MeV cố định | 2.502 | 0.3441 | 4.13 MeV | 4.80 | 3.737 | **+10.4%** |
| ²⁰⁸Pb | tỉ số 3 cố định | 3.000 | 0.4221 | 2.92 MeV | 2.77 | 2.614 | **+11.7%** |
| ²⁰⁸Pb | $E_{\text{vt}}=30$ MeV cố định | 4.335 | 0.4986 | 3.45 MeV | 2.77 | 2.614 | **+32.0%** |

**Đọc kết quả:**

1. **$\gamma_x=\text{const}$ không làm $\eta$ hằng số.** Nhiễu 1 MeV tính theo đơn vị $\hbar\omega_0$ thay đổi theo $A$ (0.061 ở O, 0.145 ở Pb), nên $\eta$ vẫn thay đổi 0.42–0.48.
2. **Hai cách đọc "8 tham số cố định" cho hai kết quả trái ngược:**
   - Giữ **tỉ số** $E_{\text{vt}}/\hbar\omega_0=3$: $\Delta E$ giảm theo $A$, nhưng lệch +12% tới +34%.
   - Giữ **$E_{\text{vt}}=30$ MeV** (đúng như bảng tham số đề xuất ghi): $\Delta E$ gần phẳng theo $A$, lệch −31% tới +32%.

   Không cách nào cho luật $16.4A^{-1/3}$ khớp 6–7% như bảng đề xuất.

## 4. Khung đã sửa

| Mục | Đề xuất | Sau kiểm chứng |
|---|---|---|
| Nguồn lái | Tương tác mạnh, $\tau\approx10$ | Số học đúng với $g^2/4\pi$. Nhưng không có kênh bức xạ pion ở 10 MeV, và dạng giả-vectơ cho $\tau\approx1900$ |
| Cơ chế gãy | Kramers, $\Delta K=12.4$ MeV | ✗. Tốc độ gãy không phụ thuộc $k_BT$ ($b=+0.013$ so với −1.675) |
| $\gamma_x$ | $=\text{const}$ vì $\lambda_{\text{mfp}}\propto A^{1/3}$ | Giả thiết không có căn cứ. Kể cả khi chấp nhận, $\eta$ vẫn không hằng số |
| Luật tỉ lệ | $16.4A^{-1/3}$, lệch 6–20% | ✗. Mô phỏng: lệch +12% tới +34% (tỉ số 3) hoặc −31% tới +32% ($E_{\text{vt}}=30$ MeV) |

## Câu hỏi còn mở

- **Nếu gãy không phải Kramers, cái gì đặt ra tốc độ ≈ 0.15 lần mỗi $T_0$?** Có thể là thời gian nạp lại của động cơ sau mỗi tái khóa. Cần đo tốc độ này theo $\tau$ và $\gamma_x$ (không theo $k_BT$) để xác định.
- **Chọn một bộ tham số duy nhất.** Tỉ số $E_{\text{vt}}/\hbar\omega_0$ cố định, hay $E_{\text{vt}}$ cố định? Hai lựa chọn cho luật tỉ lệ khác hẳn nhau. Lựa chọn nên được quyết định **trước khi** so với phổ 3⁻.
