# 火星热氧冕输运模型计划

## 1. 项目目标

本项目拟建立一个由 MGITM 背景大气驱动的火星热氧输运模型，模拟由分子氧离子解离复合

\[
\mathrm{O_2^+ + e^- \rightarrow O + O}
\]

产生的超热氧原子。主要科学输出为：

1. 热氧数密度随高度的剖面 \(n_{\mathrm{hot\,O}}(r)\)
2. 可选的二维或三维热氧冕分布
3. 热氧逃逸概率 \(P_{\mathrm{esc}}(z,E)\)
4. 光化学氧逃逸通量和全球逃逸率
5. 不同季节 \(L_s\) 和太阳活动条件下的 O 冕变化

计算内核计划使用 Julia，数据预处理、诊断和作图使用 Python。

## 2. 结论：这是 Monte Carlo 模型吗？

是，但文献中存在三个相关而不完全相同的实现。

### 2.1 Lillis et al. (2017)

Lillis et al. 使用三维 Monte Carlo 热原子输运模型。模型在给定产生高度生成各向同性的热 O，随机抽取初始能量和方向，随后计算其在球对称中性大气中的碰撞、散射、能量损失和重力运动。该模型的直接目标是得到随产生高度变化的逃逸概率。

他们在每个高度持续产生粒子，直到有 2500 个粒子逃逸。逃逸概率定义为

\[
P_{\mathrm{esc}}(z)=
\frac{N_{\mathrm{escaped}}}{N_{\mathrm{spawned}}}.
\]

这一实现适合计算逃逸概率和逃逸通量，但仅记录逃逸与否不足以直接给出完整 O 冕密度剖面。

### 2.2 Rahmati (2016)

Ali Rahmati 的博士论文同时使用：

1. 一维 two-stream/Liouville 模型，用于快速计算热 O 通量、逃逸率和分布函数
2. 三维 Monte Carlo 模型，用于完整追踪热 O 的碰撞区和无碰撞区轨迹，并验证 two-stream 结果

其 Monte Carlo 模型包含弹道轨道、逃逸轨道和极少量卫星轨道，并在球壳上记录粒子的位置和速度以重建分布函数。论文中的模型将外边界延伸到约 110,000 km。

### 2.3 Lee et al. (2015)

Lee et al. 将 MGITM 与 M-AMPS 粒子模型单向耦合，直接模拟热氧冕及其逃逸。当前文件夹中的 12 组 MGITM 数据与该论文的四个季节和三个太阳活动水平设计一致。因此，对于本项目的“O 冕剖面”目标，Lee et al. 的 MGITM 与 Monte Carlo 耦合框架是最直接的总体参考，Lillis et al. 和 Rahmati 的工作适合确定产生率、碰撞物理和数值验证方法。

## 3. 当前本地材料

### 3.1 文献

1. `Lillis-Photochemical escape of oxygen from Mar.pdf`
2. `Rahmati_ku_0099D_14448_DATA_1.pdf`

### 3.2 MGITM 数据

`MGITM` 文件夹包含 12 个三维背景大气文件：

1. 四个季节：\(L_s=0^\circ,90^\circ,180^\circ,270^\circ\)
2. 三个太阳活动水平：F070、F130、F200
3. 网格：72 个经度点，36 个纬度点，62 个高度点
4. 水平分辨率：\(5^\circ \times 5^\circ\)
5. 垂直分辨率：2.5 km
6. 高度范围：98.75 至 251.25 km
7. 数密度单位：\(\mathrm{m^{-3}}\)

每个文件包含：

\[
T_n,\ T_i,\ T_e,\ n_{\mathrm{CO_2}},\ n_{\mathrm O},\ n_{\mathrm{N_2}},
n_{\mathrm{CO}},\ n_{\mathrm{O_2}},\ n_{\mathrm{O_2^+}},
n_{\mathrm{O^+}},\ n_{\mathrm{CO_2^+}},\ n_e
\]

以及中性风 \(U_N,V_N,W_N\)。这些变量足以计算 O2+ 解离复合产生率并构造主要背景碰撞种类。

## 4. 热 O 产生模型

### 4.1 解离复合反应率和热 O 产生率

本项目固定采用用户提供的分段表达式。O2+ 解离复合的体积反应率定义为

\[
\alpha(z)=
\begin{cases}
1.95\times10^{-7}\,
n_e(z)n_{\mathrm{O_2^+}}(z)
\left(\dfrac{300}{T_e(z)}\right)^{0.70},
& T_e<1200\ \mathrm{K},\\[8pt]
7.39\times10^{-8}\,
n_e(z)n_{\mathrm{O_2^+}}(z)
\left(\dfrac{1200}{T_e(z)}\right)^{0.56},
& T_e>1200\ \mathrm{K}.
\end{cases}
\]

当 \(n_e\) 和 \(n_{\mathrm{O_2^+}}\) 使用 \(\mathrm{cm^{-3}}\) 时，
\(\alpha\) 的单位为 \(\mathrm{cm^{-3}\,s^{-1}}\)。在数值实现中将
\(T_e=1200\ \mathrm{K}\) 归入低温或高温分支均不会产生跳变，建议明确采用
\(T_e\leq1200\ \mathrm{K}\) 的低温分支。

每次解离复合反应产生两个 O 原子，因此热 O 原子的总产生率为

\[
Q_{\mathrm{hot\,O}}(z)=2\alpha(z).
\]

也可先定义不含密度的反应率系数

\[
k(T_e)=
\begin{cases}
1.95\times10^{-7}
\left(\dfrac{300}{T_e}\right)^{0.70},
& T_e\leq1200\ \mathrm{K},\\[8pt]
7.39\times10^{-8}
\left(\dfrac{1200}{T_e}\right)^{0.56},
& T_e>1200\ \mathrm{K},
\end{cases}
\]

其中 \(k\) 的单位为 \(\mathrm{cm^3\,s^{-1}}\)，并使用

\[
\alpha=n_en_{\mathrm{O_2^+}}k(T_e).
\]

实现时必须统一单位。MGITM 数密度为 \(\mathrm{m^{-3}}\)，建议 Julia 内部全部采用 SI 单位。转换关系为

\[
1\ \mathrm{cm^3\,s^{-1}}=10^{-6}\ \mathrm{m^3\,s^{-1}}.
\]

### 4.2 解离复合支路

本项目固定采用以下四个非忽略反应支路和分支比：

| 产物 | 释放能量 | 分支比 |
|---|---:|---:|
| \(\mathrm{O(^3P)+O(^3P)}\) | 6.99 eV | 26.5% |
| \(\mathrm{O(^1D)+O(^3P)}\) | 5.02 eV | 47.3% |
| \(\mathrm{O(^1D)+O(^1D)}\) | 3.06 eV | 20.4% |
| \(\mathrm{O(^1D)+O(^1S)}\) | 0.83 eV | 5.8% |

四个分支比之和为 100%。Monte Carlo 抽样时使用累计概率

\[
[0,0.265),\quad
[0.265,0.738),\quad
[0.738,0.942),\quad
[0.942,1].
\]

第一版模型可先在质心系中将释放能量平均分配给两个相同质量的 O 原子。对应的单个 O 初始动能分别约为

\[
3.495,\quad 2.510,\quad 1.530,\quad 0.415\ \mathrm{eV}.
\]

之后再加入 O2+ 热运动和电子热运动，使初始能谱随 \(T_i\) 和 \(T_e\) 展宽。两个产物 O 在质心系中的速度方向相反，因此实现时应将一次反应的两个 O 作为相关粒子生成，或者使用具有等价统计权重的单粒子抽样方法。

## 5. 碰撞与输运物理

### 5.1 背景碰撞种类

第一版至少包含：

1. O 与 CO2
2. O 与 O
3. O 与 N2
4. O 与 CO

Lillis et al. 使用的代表性总截面为：

| 碰撞 | 总截面 |
|---|---:|
| O 与 CO2 | \(2.0\times10^{-14}\ \mathrm{cm^2}\) |
| O 与 O | \(0.6\times10^{-14}\ \mathrm{cm^2}\) |
| O 与 N2 | \(1.8\times10^{-14}\ \mathrm{cm^2}\) |
| O 与 CO | \(1.8\times10^{-14}\ \mathrm{cm^2}\) |

Rahmati 使用能量依赖形式

\[
\sigma_s(E)=\alpha_s E^{-0.2},
\]

并采用 Kharchenko et al. (2000) 的强前向峰化微分截面。建议先实现常数总截面版本作为基线，再实现能量和散射角依赖版本。

### 5.2 Monte Carlo 自由程

局地总碰撞频率因子为

\[
\Lambda^{-1}(\mathbf r,E)
=\sum_s n_s(\mathbf r)\sigma_s(E).
\]

可使用事件驱动抽样：

\[
\tau=-\ln \xi,\qquad
\int_0^\ell \sum_s n_s(\mathbf r(l))\sigma_s(E)\,dl=\tau,
\]

其中 \(\xi\) 为 \(0\) 至 \(1\) 的均匀随机数。该方法比固定小步长碰撞概率法更精确，也更适合 Julia 并行计算。

碰撞种类按

\[
P_s=
\frac{n_s\sigma_s}{\sum_j n_j\sigma_j}
\]

随机选择。碰撞后的能量和方向由质心系散射角以及动量、能量守恒确定。

### 5.3 重力与边界

粒子在碰撞之间按照火星中心引力运动：

\[
\frac{d^2\mathbf r}{dt^2}
=-\frac{GM_{\mathrm M}}{r^3}\mathbf r.
\]

需要定义：

1. 下边界，建议第一版使用 100 或 120 km，到达后视为热化并吸收
2. 外边界，调试时可先取数个火星半径，正式 O 冕计算再扩大
3. 热化阈值，例如 \(E<0.01\) eV，需通过敏感性测试确认
4. 逃逸判据，应同时检查外边界和总机械能，而不只检查局地速度

## 6. 如何得到 O 冕密度剖面

这是本项目与只计算逃逸概率的模型之间最重要的区别。

对于稳态源，每个模拟粒子代表一定的物理产生率权重 \(w_i\)。粒子穿过径向球壳 \(j\) 时，在该球壳内的驻留时间为 \(\Delta t_{ij}\)。热 O 数密度可由驻留时间估计：

\[
n_{\mathrm{hot\,O}}(r_j)
=\frac{1}{V_j}
\sum_i w_i\Delta t_{ij},
\]

其中

\[
V_j=\frac{4\pi}{3}
\left(r_{j+1}^3-r_j^3\right).
\]

若只模拟局地柱或有限立体角，则 \(V_j\) 必须替换为对应网格体积。权重 \(w_i\) 来自 MGITM 网格单元内的产生率：

\[
w_i =
\frac{Q_{\mathrm{hot\,O},k}V_k}{N_{\mathrm{test},k}}.
\]

同时累计速度或能量直方图，可以得到

\[
f(r,E),\quad n(r,E),\quad
n_{\mathrm{ballistic}}(r),\quad
n_{\mathrm{escaping}}(r).
\]

因此第一版 Monte Carlo 输出必须包含粒子驻留时间，而不能只输出最终状态。

## 7. 建议的建模层级

### 阶段 A：一维球对称基线模型

从每个 MGITM 文件中提取接近日下点的垂直剖面，建立球对称背景。该模型用于：

1. 快速验证产生率
2. 重现 Lillis 的逃逸概率随高度变化
3. 与 Rahmati 的热 O 分布函数和逃逸率数量级比较
4. 得到第一组 O 冕径向剖面

这是最适合立即开始的版本。

### 阶段 B：二维轴对称模型

按太阳天顶角对 MGITM 场进行平均，保留昼夜差异，得到

\[
n_{\mathrm{hot\,O}}(r,\mathrm{SZA}).
\]

### 阶段 C：完整三维 MGITM 耦合

直接使用经度、纬度和高度网格，模拟季节、地方时、中性风和全球非对称性。这一阶段最接近 Lee et al. 的 M-GITM 与 M-AMPS 框架。

## 8. Julia 与 Python 分工

### Julia

建议模块：

```text
HotOTransport.jl
├── Constants.jl
├── MGITMReader.jl
├── BackgroundAtmosphere.jl
├── DissociativeRecombination.jl
├── CrossSections.jl
├── ParticleSource.jl
├── CollisionSampler.jl
├── Trajectory.jl
├── Tallies.jl
└── Simulation.jl
```

Julia 负责：

1. MGITM 数据读取或读取 Python 预处理后的标准文件
2. 热 O 产生与测试粒子加权
3. 碰撞和轨迹推进
4. 多线程或分布式 Monte Carlo
5. 驻留时间、能谱和逃逸统计
6. 将结果写为 HDF5 或 NetCDF

随机数种子、物理参数和输入文件名必须写入输出元数据，保证可重复性。

### Python

Python 负责：

1. 检查和可视化 MGITM 输入
2. 绘制 \(Q_{\mathrm{hot\,O}}(z)\)
3. 绘制 \(P_{\mathrm{esc}}(z,E)\)
4. 绘制 O 冕密度剖面和能谱
5. 比较 12 个季节与太阳活动案例
6. 绘制 Monte Carlo 收敛性和不确定度

所有图中文字使用 Arial，数学公式保留数学字体。

## 9. 第一版数值验证

第一版不能直接追求完整三维结果，应先通过以下测试：

1. 无碰撞两体轨道测试，机械能和角动量守恒
2. 均匀背景中的自由程分布满足指数分布
3. 碰撞前后总动量和总能量守恒
4. O 与不同质量靶粒子的平均能量损失合理
5. 无碰撞、各向同性源的最大直接逃逸比例符合几何预期
6. 逃逸概率随产生高度和初始能量总体单调增加
7. 增加粒子数时，密度和逃逸率按 \(N^{-1/2}\) 收敛
8. 常数截面结果与能量依赖截面结果的差异可解释
9. 计算结果对下边界、外边界和热化阈值不敏感
10. 太阳高活动条件下的 O 冕总体高于太阳低活动条件

## 10. 主要风险与需要谨慎处理的问题

1. O 与 CO2 的弹性碰撞截面是不确定度最大的输入之一。Lillis et al. 指出其误差可能接近因子 2。
2. MGITM 顶边界只有 251.25 km。热 O 源大部分位于约 180 至 280 km，因此需要谨慎外推 O2+、电子和中性背景，或证明 251.25 km 以上源项可忽略。
3. MGITM 中的热 O 产生率是三维且随地方时变化的。日下点球对称模型只是基线，不应被解释为完整全球 O 冕。
4. O 与 O 碰撞会产生次级热 O。若忽略反冲 O，逃逸概率和高空密度可能偏低。
5. 稳态密度必须使用具有绝对物理权重的驻留时间累计。仅统计穿越次数不能直接得到绝对数密度。
6. 高空弹道粒子的驻留时间很长，需要避免少数高权重轨迹造成过大 Monte Carlo 方差。

## 11. 推荐的第一步实现范围

第一版建议选择 `MGITM_LS000_F070_150901.dat` 的日下点附近剖面，并完成：

1. 读取和单位转换
2. 计算 O2+ 解离复合产生率
3. 四个反应支路的初始能量抽样
4. 球对称背景中的 O 与 CO2、O、N2、CO 碰撞
5. 火星引力轨迹
6. 径向驻留时间累计
7. 热 O 密度剖面、逃逸概率和逃逸率
8. Python 诊断图

基线通过后再批量运行 12 个 MGITM 条件。

## 12. 参考文献

1. Lillis, R. J., et al. (2017), Photochemical escape of oxygen from Mars: First results from MAVEN in situ data, *Journal of Geophysical Research: Space Physics*, 122, 3815-3836, https://doi.org/10.1002/2016JA023525.
2. Rahmati, A. (2016), *Oxygen Exosphere of Mars: Evidence from Pickup Ions Measured by MAVEN*, PhD dissertation, University of Kansas.
3. Lee, Y., Combi, M. R., Tenishev, V., Bougher, S. W., and Lillis, R. J. (2015), Hot oxygen corona at Mars and the photochemical escape of oxygen: Improved description of the thermosphere, ionosphere, and exosphere, *Journal of Geophysical Research: Planets*, 120, 1880-1892, https://doi.org/10.1002/2015JE004890.
4. Fox, J. L., and Hać, A. B. (2009), Photochemical escape of oxygen from Mars: A comparison of the exobase approximation to a Monte Carlo method, *Icarus*, 204, 527-544, https://doi.org/10.1016/j.icarus.2009.07.005.
5. Fox, J. L., and Hać, A. B. (2014), The escape of O from Mars: Sensitivity to the elastic cross sections, *Icarus*, 228, 375-385, https://doi.org/10.1016/j.icarus.2013.10.014.
6. Kharchenko, V., et al. (2000), Angular and energy-dependent O-O collision data used by the referenced hot oxygen transport models.
