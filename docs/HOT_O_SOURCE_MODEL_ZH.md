# 热 O 高度和初生能量分布

## 1. 这一部分要计算什么

源模型要得到两个量：

1. 每个高度每秒产生多少个热 O，即 $Q_{\mathrm{hotO}}(z)$。
2. 在该高度产生的热 O 具有怎样的能量概率分布，即 $p(E\mid z)$。

二者相乘得到高度和能量的二维产生率：

```math
Q(E,z)=Q_{\mathrm{hotO}}(z)\,p(E\mid z).
```

其单位为 $\mathrm{cm^{-3}\,s^{-1}\,eV^{-1}}$。对能量积分后应回到总产生率：

```math
\int Q(E,z)\,dE=Q_{\mathrm{hotO}}(z).
```

这里的 $Q(E,z)$ 是刚产生时的源分布，不是经过碰撞传输后的 O 冕密度 $n(E,z)$。

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

重复次数足够多以后，将动能放入能量直方图并归一化，就得到 $p(E\mid z)$。

Monte Carlo 不是为了计算总产生率。总产生率由密度和反应系数直接计算。Monte Carlo 的作用是把反应分支、热速度、振动态和随机方向共同转换成 LAB 系能量概率分布。

## 5. 如何由 $T_e$ 和 $T_i$ 抽取速度

温度为 $T$、质量为 $m$ 的 Maxwell 分布，可以通过独立抽取三个笛卡尔速度分量实现：

```math
v_x,v_y,v_z\sim
\mathcal{N}\left(0,\frac{k_{\mathrm B}T}{m}\right).
```

也就是说，每个分量的标准差是

```math
\sigma_v=\sqrt{\frac{k_{\mathrm B}T}{m}}.
```

电子质量很小，所以相同温度下电子速度远大于 $\mathrm{O_2^+}$ 速度。但是反应运动学还要乘以质量，因此不能只凭速度大小判断其对最终 LAB 能量展宽的贡献。

设抽到的速度为 $\mathbf v_e$ 和 $\mathbf v_i$。反应物质心速度为

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

对每个高度分别生成大量事件，然后执行以下归一化：

```math
p(E_k\mid z)
\approx
\frac{N_k(z)}
{N_{\mathrm{tot}}(z)\Delta E},
```

其中 $N_k$ 是落入第 $k$ 个能量区间的样本数。

再计算

```math
Q(E_k,z)=Q_{\mathrm{hotO}}(z)p(E_k\mid z).
```

绘图时：

* 横坐标是初生 O 的 LAB 能量。
* 纵坐标是产生高度。
* 颜色是 $\log_{10}Q(E,z)$。
* 图像插值只用于让色块显示平滑，不会增加新的物理信息。

![热 O 初生能量图](../examples/figures/mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_energy_maps.png)

## 9. 与 Ali 和 Lillis 模型的关系

Ali 和 Lillis 所用源模型的核心思想相同：

1. 由解离复合率确定高度源强。
2. 由分支能量确定主要能量峰。
3. 由反应物热速度和两体反应运动学计算 LAB 系初生能量。
4. 将得到的 $Q(E,z)$ 输入热原子传输模型。

当前 MarsHotO 在此基础上显式加入了给定的 $\mathrm{O_2^+}$ 振动态布居。源模型产生的是初始粒子，后续碰撞不会反过来改变源项，而是把 $Q(E,z)$ 转换为 O 冕的 $n(E,z)$、$n(z)$ 和速度分布。

## 10. 最重要的区别

```text
Q_hotO(z)      每单位体积、每秒产生多少个热 O
p(E|z)         新产生的热 O 能量概率密度
Q(E,z)         新产生热 O 的高度和能量谱
n(E,z)         经过运动和碰撞后，实际存在的热 O 高度和能量分布
```

$Q(E,z)$ 由化学源模型得到，$n(E,z)$ 必须通过热 O 传输 Monte Carlo 模型计算。
