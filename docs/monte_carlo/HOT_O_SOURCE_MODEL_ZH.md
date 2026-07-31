# 热 O 高度和初生能量分布

## 1. 这一部分要计算什么

源模型需要得到三个量：

1. 每个高度、每单位体积、每秒产生多少个热 O，即总产生率
   $Q_{\mathrm{hotO}}(z)$。
2. 在该高度产生的热 O 落入各能量格的概率，即 $P_k(z)$。
3. 每个能量格内的热 O 产生率，即 $Q_k(z)$。

对于能量范围 $[E_k,E_k+\Delta E)$，分箱概率定义为

```math
P_k(z)=
\frac{N_k(z)}{N_{\mathrm{tot}}(z)}.
```

$P_k(z)$ 没有单位，并满足

```math
\sum_k P_k(z)=1.
```

当前图的左侧面板画的是 $P_k(z)$，因此左图不除以能量格宽度
$\Delta E$。

右图使用

```math
Q_k(z)=
Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z)P_k(z).
```

因此，右图也不除以 $\Delta E$。所有能量格的产生率之和等于该高度的
总产生率：

```math
\sum_k Q_k(z)=Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z).
```

这里的 $Q_k(z)$ 是刚产生时的源分布，不是经过碰撞传输后的 O 冕密度。

## 2. 高度产生率

主要光化学源是解离复合反应

```math
\mathrm{O_2^+}+e\rightarrow\mathrm{O}+\mathrm{O}.
```

反应事件率为

```math
R_{\mathrm{DR}}(z)=
n_e(z)n_{\mathrm{O_2^+}}(z)\alpha[T_e(z)].
```

每次反应产生两个 O，所以热 O 原子产生率为

```math
Q_{\mathrm{hotO}}(z)=
2n_e(z)n_{\mathrm{O_2^+}}(z)\alpha[T_e(z)].
```

当前使用的反应系数为

```math
\alpha(T_e)=
\begin{cases}
1.95\times10^{-7}
\left(\dfrac{300}{T_e}\right)^{0.70},
& T_e<1200\ \mathrm{K},\\[6pt]
7.39\times10^{-8}
\left(\dfrac{1200}{T_e}\right)^{0.56},
& T_e>1200\ \mathrm{K},
\end{cases}
```

单位为 $\mathrm{cm^3\,s^{-1}}$。

因此，$n_e$ 和 $n_{\mathrm{O_2^+}}$ 主要控制产生率的高度形状，$T_e$ 通过反应系数进一步调节产生率。

下图给出默认 MGITM 基准大气的 $\mathrm{O_2^+}$ 数密度、温度以及热 O 总体积产生率。右侧面板 c 是本节使用的 $Q_{\mathrm{hotO}}(z)$，单位为 $\mathrm{cm^{-3}\,s^{-1}}$。

![MGITM 热 O 高度产生率](../../examples/figures/mgitm_ls000_f070_profiles.png)

## 3. 四个反应分支

模型使用四个不可忽略的解离复合分支。

| 产物 | 总释放能量 | 分支概率 | 两个 O 平分时单个 O 的基本能量 |
|---|---:|---:|---:|
| O($^3P$) + O($^3P$) | 6.99 eV | 26.5% | 3.495 eV |
| O($^1D$) + O($^3P$) | 5.02 eV | 47.3% | 2.510 eV |
| O($^1D$) + O($^1D$) | 3.06 eV | 20.4% | 1.530 eV |
| O($^1D$) + O($^1S$) | 0.83 eV | 5.8% | 0.415 eV |

如果电子和 $\mathrm{O_2^+}$ 都完全静止，而且不考虑振动能，那么每个分支会在能量和高度图上形成一条很细的竖直线。

真实分布不是四条无限细的线，因为反应物在反应前具有热速度，而且 $\mathrm{O_2^+}$ 可能处于不同振动态。

## 4. Monte Carlo 在这里做什么

对于 MGITM 的每一个高度 $z$，重复生成许多次解离复合事件。每次事件依次进行：

1. 按分支概率抽取一个反应通道。
2. 根据 $T_e(z)$ 抽取一个电子速度。
3. 根据 $T_i(z)$ 抽取一个 $\mathrm{O_2^+}$ 速度。
4. 按振动态布居抽取振动量子数 $v$。
5. 在反应物质心系中计算两个 O 的相对速度。
6. 抽取一个各向同性的产物方向。
7. 把两个 O 的速度转换回 LAB 系。
8. 记录 LAB 系中每个 O 的动能。

重复次数足够多以后，将动能放入能量直方图。每个能量格的计数除以总计数，
得到左图使用的分箱概率 $P_k(z)$。将 $P_k(z)$ 乘以该高度的总产生率，
得到右图使用的每格产生率 $Q_k(z)$。

Monte Carlo 不是为了计算总产生率。总产生率由密度和反应系数直接计算。Monte Carlo 的作用是把反应分支、热速度、振动态和随机方向共同转换成 LAB 系能量概率分布。

## 5. 如何由 $T_e$ 和 $T_i$ 抽取速度

当前模型把电子和 $\mathrm{O_2^+}$ 的背景整体速度都设为零：

```math
\mathbf u_{e,\mathrm{bulk}}
=
\mathbf u_{i,\mathrm{bulk}}
=0.
```

温度只决定相对于零整体速度的热运动。电子和 $\mathrm{O_2^+}$ 都使用归一化
的三维 Maxwellian 速度分布：

```math
f_s(\mathbf v\mid T_s)
=
\left(\frac{m_s}{2\pi k_{\mathrm B}T_s}\right)^{3/2}
\exp\left[
-\frac{m_s|\mathbf v-\mathbf u_s|^2}
{2k_{\mathrm B}T_s}
\right].
```

这里 $\mathbf u_s$ 是整体速度。当前电子和 $\mathrm{O_2^+}$ 均设置为

```math
\mathbf u_s=(0,0,0).
```

该概率密度在完整三维速度空间中的积分为 1：

```math
\int_{\mathbb R^3}
f_s(\mathbf v\mid T_s)\,d^3v=1.
```

如果使用物理数密度形式 $n_sf_s$，其速度空间积分为 $n_s$。源粒子的 Monte
Carlo 只需要归一化概率分布，所以抽样本身不乘数密度。数密度已经包含在
$Q_{\mathrm{hotO}}(z)$ 和宏粒子权重中。

MarsHotO 使用与 TestParticle.jl 相同的热速度约定：

```math
v_{\mathrm{th},s}
=
\sqrt{\frac{2k_{\mathrm B}T_s}{m_s}}.
```

这对应 TestParticle.jl 的构造方式：

```julia
u_bulk = [0.0, 0.0, 0.0]
p = n * kB * T
vdf = TP.Maxwellian(u_bulk, p, n; m=mass)
```

MarsHotO 不需要依赖 TestParticle.jl，包内的
`sample_maxwellian_velocity` 使用相同的数学定义直接完成抽样。

每个笛卡尔速度分量的标准差为

```math
\sigma_{v,s}
=
\frac{v_{\mathrm{th},s}}{\sqrt{2}}
=
\sqrt{\frac{k_{\mathrm B}T_s}{m_s}}.
```

因此，对每一种反应物分别抽取三个独立标准正态随机数：

```math
\xi_x,\xi_y,\xi_z\sim\mathcal N(0,1),
```

并计算

```math
\mathbf v_s
=
\mathbf u_s
+
\sqrt{\frac{k_{\mathrm B}T_s}{m_s}}
(\xi_x,\xi_y,\xi_z).
```

当 $\mathbf u_s=0$ 时，三个分量具有相同方差，所以所得三维分布是各向同性的。
速度方向会自然覆盖整个球面，不需要另外抽取极角和方位角。

对应的速度大小概率密度为

```math
P_s(v)
=
4\pi v^2
\left(\frac{m_s}{2\pi k_{\mathrm B}T_s}\right)^{3/2}
\exp\left(-\frac{m_sv^2}{2k_{\mathrm B}T_s}\right),
\qquad v\ge0.
```

总动能概率密度为

```math
p_s(E\mid T_s)
=
\frac{2}{\sqrt{\pi}}
\frac{\sqrt{E}}{(k_{\mathrm B}T_s)^{3/2}}
\exp\left(-\frac{E}{k_{\mathrm B}T_s}\right),
\qquad E\ge0.
```

它满足

```math
\int_0^\infty p_s(E\mid T_s)\,dE=1,
\qquad
\langle E_s\rangle=\frac{3}{2}k_{\mathrm B}T_s.
```

300 K 条件下的速度分量、总动能和方向余弦抽样见下图：

![300 K Maxwellian 速度抽样](../../examples/figures/thermal_energy_sampling_300K.png)

电子质量远小于 $\mathrm{O_2^+}$ 质量，所以相同温度和相同动能下电子速度更大。分别抽到 $\mathbf v_e$ 和 $\mathbf v_i$ 后，反应物质心速度为

```math
\mathbf V_{\mathrm{COM}}=
\frac{m_e\mathbf v_e+m_i\mathbf v_i}{m_e+m_i}.
```

反应物相对动能为

```math
E_{\mathrm{rel}}=
\frac{1}{2}\mu
\left|\mathbf v_e-\mathbf v_i\right|^2,
\qquad
\mu=\frac{m_em_i}{m_e+m_i}.
```

$T_e$ 和 $T_i$ 随高度变化，所以 $\mathbf V_{\mathrm{COM}}$ 和 $E_{\mathrm{rel}}$ 的统计分布也随高度变化。这就是能量峰宽度随高度变化的来源之一。

## 6. 振动态如何加入

当前配置采用振动量子间隔

```math
\Delta E_{\mathrm{vib}}=0.23\ \mathrm{eV}
```

以及 $v=0$ 到 $8$ 的布居比例：

```text
0.800, 0.074, 0.043, 0.035, 0.025,
0.015, 0.0047, 0.00027, 0.00021
```

代码先将这些比例归一化，再随机抽取 $v$。该事件携带的附加振动能为

```math
E_{\mathrm{vib}}=v\Delta E_{\mathrm{vib}}.
```

对释放能为 $E_b$ 的反应分支，总可用平动能为

```math
E_{\mathrm{avail}}=
E_b+E_{\mathrm{rel}}+E_{\mathrm{vib}}.
```

振动能使每个基本反应峰向高能侧产生附加结构和展宽。

这是当前模型采用的能量预算近似。更精细的模型还可以让振动态影响不同解离通道的分支比例。

## 7. 质心系中的两个 O

两个产物都是 O，质量相同。在产物质心系中，它们的速度大小相同、方向相反，因此每个 O 得到一半的总可用能量：

```math
E_{\mathrm O,COM}=\frac{E_{\mathrm{avail}}}{2}.
```

每个 O 的速度大小为

```math
u=
\sqrt{\frac{E_{\mathrm{avail}}}{m_{\mathrm O}}},
```

其中 $E_{\mathrm{avail}}$ 在计算时需要从 eV 转换为 J。

抽取一个各向同性单位向量 $\hat{\mathbf n}$ 后，两个产物的 LAB 速度为

```math
\mathbf v_{\mathrm O,1}
=\mathbf V_{\mathrm{COM}}+u\hat{\mathbf n},
```

```math
\mathbf v_{\mathrm O,2}
=\mathbf V_{\mathrm{COM}}-u\hat{\mathbf n}.
```

最后计算 LAB 动能：

```math
E_{\mathrm O,LAB}
=\frac{1}{2}m_{\mathrm O}
\left|\mathbf v_{\mathrm O}\right|^2.
```

同一个 COM 能量在 LAB 系中不再对应唯一能量，因为 $\mathbf V_{\mathrm{COM}}$ 与 $\hat{\mathbf n}$ 的夹角每次不同。这会进一步展宽能量峰。

## 8. 如何得到高度和能量二维图

对每个高度分别生成大量事件。设 $N_k(z)$ 是落入第 $k$ 个能量格的热 O
样本数，$N_{\mathrm{tot}}(z)$ 是该高度的总样本数，当前能量格宽度为

```math
\Delta E=0.025\ \mathrm{eV}.
```

### 8.1 左图：每个能量格的概率

左图直接使用

```math
P_k(z)=
\frac{N_k(z)}
{N_{\mathrm{tot}}(z)}.
```

它表示一个新产生的热 O 落入第 $k$ 个 0.025 eV 能量格的概率。该量没有
单位，并且在每个高度满足

```math
\sum_k P_k(z)=1.
```

因此左图不除以 $\Delta E$。其 colorbar 标注为
`Probability per 0.025 eV bin`。

### 8.2 右图：每个能量格的产生率

右图使用

```math
Q_k(z)=
Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z)P_k(z).
```

MGITM 源剖面中的总产生率最初以 $\mathrm{cm^{-3}\,s^{-1}}$ 计算。绘制
右图前将其转换为 $\mathrm{m^{-3}\,s^{-1}}$：

```math
Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z)
=
10^6Q_{\mathrm{hotO}}^{(\mathrm{cm^{-3}})}(z).
```

右图中的每个数值表示第 $k$ 个 0.025 eV 能量格内，每立方米、每秒产生的
热 O 数量。它的单位是 $\mathrm{m^{-3}\,s^{-1}}$，并满足

```math
\sum_k Q_k(z)
=
Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z).
```

左图和右图都不除以 $\Delta E$。

图中：

* 横坐标是初生 O 的 LAB 能量。
* 纵坐标是产生高度。
* 左图颜色是分箱概率 $P_k(z)$。
* 右图颜色是 $\log_{10}Q_k(z)$，$Q_k$ 的单位为
  $\mathrm{m^{-3}\,s^{-1}}$。
* 图像插值只用于让色块显示平滑，不会增加新的物理信息。

![热 O 初生能量图](../../examples/figures/mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_energy_maps.png)

## 9. 与 Rahmati 和 Lillis 模型的关系

Rahmati 和 Lillis 所用源模型的核心思想相同：

1. 由解离复合率确定高度源强。
2. 由分支能量确定主要能量峰。
3. 由反应物热速度和两体反应运动学计算 LAB 系初生能量。
4. 在每个高度直接生成带权重的解离复合事件，每个事件产生一对 O，随后分别调用单粒子输运模型。

当前 MarsHotO 显式加入给定的 $\mathrm{O_2^+}$ 振动态布居。$P_k(z)$ 和 $Q_k(z)$ 是事件采样得到的诊断量，用于绘图和源模型验证，不作为输运程序再次独立抽取单个 O 能量的输入。这样可以保留同一反应中两个 O 的共同分支、振动态、COM 速度和严格反向关系。

## 10. 最重要的区别

```text
Q_hotO(z)      每单位体积、每秒产生多少个热 O
P_k(z)         新产生的热 O 落入第 k 个能量格的概率，无量纲
Q_k(z)         第 k 个能量格内的产生率，单位为 m^-3 s^-1
O 冕分布       经过运动和碰撞后，实际存在的热 O 高度、能量和速度分布
```

$Q_k(z)$ 由化学事件的统计直方图得到，O 冕分布由同一批带权重事件产生的 O 直接进行输运得到。

## 11. 输运 Monte Carlo 的宏粒子权重

默认源高度为 100 至 250 km，间隔为 1 km。定义解离复合事件率

```math
R_{\mathrm{DR}}(z_i)=n_e(z_i)n_{\mathrm{O_2^+}}(z_i)k[T_e(z_i)].
```

每个高度生成

```math
N_{\mathrm{event},i}=10000
```

个反应事件。每个事件产生两个初级热 O，因此总初级粒子数仍为

```math
151\times20000=3.02\times10^6.
```

第 $i$ 个源高度对应的球壳体积为

```math
V_i=
\frac{4\pi}{3}
\left[
(R_{\mathrm M}+z_{i,+})^3
-
(R_{\mathrm M}+z_{i,-})^3
\right].
```

该球壳每秒发生的真实反应事件数为

```math
S_{\mathrm{event},i}=R_{\mathrm{DR}}(z_i)V_i.
```

每个模拟事件的权重为

```math
w_{\mathrm{event},i}=
\frac{R_{\mathrm{DR}}(z_i)V_i}{N_{\mathrm{event},i}}.
```

$w_{\mathrm{event},i}$ 的单位是 $\mathrm{s^{-1}}$。每个事件独立抽样电子和
$\mathrm{O_2^+}$ 的 Maxwell 热速度。bulk velocity 设为零仅表示两个 Maxwell
分布的均值为零，每个事件的反应物 COM 速度通常不为零。两个产物在该 COM
系中严格反向，并且都继承同一个事件权重。随后两个 O 分别调用单粒子输运。
如果碰撞靶粒子是 O，反冲 O 作为次级热 O 加入队列，并继承母粒子的权重。

## 12. 如何由轨迹得到 O 冕高度和能量分布

当权重为 $w_p$ 的粒子在高度格 $i$ 和能量格 $k$ 中停留
$\Delta t_p$ 时，驻留时间估计器累加

```math
C_{ik}\mathrel{+}=w_p\Delta t_p.
```

$C_{ik}$ 表示稳态条件下该高度和能量格中存在的真实热 O 数量。完成所有
轨迹后，除以诊断球壳体积：

```math
n_{ik}=
\frac{\sum_p w_p\Delta t_p}
{V_i}.
```

$n_{ik}$ 的单位是 $\mathrm{m^{-3}}$ per energy bin。当前输出不除以
能量格宽度。对能量格求和得到总热 O 数密度：

```math
n_i=\sum_k n_{ik}.
```

如果需要每个高度归一化后的 O 冕能量概率，则计算

```math
F_{ik}=
\frac{n_{ik}}
{\sum_j n_{ij}}.
```

$F_{ik}$ 是经过重力和碰撞输运后的能量概率，不等于初生源概率
$P_k(z_i)$。
