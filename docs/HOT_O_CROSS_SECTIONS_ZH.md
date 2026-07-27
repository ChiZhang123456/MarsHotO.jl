# 热 O 与中性大气的碰撞截面

## 1. 碰撞截面是什么

碰撞截面 $\sigma$ 可以直观理解为一个粒子被另一个粒子“碰到”的有效面积。它不是原子的真实几何面积，而是对某种碰撞过程发生可能性的定量描述。

对于能量为 $E$ 的热 O，目标成分 $s$ 的局地碰撞系数为

```math
\kappa_s(E,z)=n_s(z)\sigma_s(E),
```

其中：

* $n_s(z)$ 是目标成分的数密度，单位为 $\mathrm{m^{-3}}$。
* $\sigma_s(E)$ 是热 O 与该成分的总碰撞截面，单位为 $\mathrm{m^2}$。
* $\kappa_s$ 的单位为 $\mathrm{m^{-1}}$。

所有中性成分的贡献相加：

```math
\kappa_{\mathrm{tot}}(E,z)=\sum_s n_s(z)\sigma_s(E).
```

平均自由程为

```math
\lambda(E,z)=\frac{1}{\kappa_{\mathrm{tot}}(E,z)}.
```

所以，中性密度越高或碰撞截面越大，平均自由程越短，热 O 越容易碰撞。

## 2. 当前采用的总碰撞截面

MarsHotO 当前在 `data/cross_sections/ali_total_cross_sections.toml` 中设置了以下 3 eV 参考截面。

| 目标成分 | $\sigma(3\ \mathrm{eV})$，$\mathrm{cm^2}$ | 当前 MGITM 中是否存在 |
|---|---:|---|
| O | $6.4\times10^{-15}$ | 是 |
| CO | $1.8\times10^{-14}$ | 是 |
| N₂ | $1.8\times10^{-14}$ | 是 |
| O₂ | $1.8\times10^{-14}$ | 是 |
| Ar | $1.2\times10^{-14}$ | 否 |
| CO₂ | $2.0\times10^{-14}$ | 是 |

单位换算为

```math
1\ \mathrm{cm^2}=10^{-4}\ \mathrm{m^2}.
```

当前采用 Ali 模型中的近似能量依赖：

```math
\sigma_s(E)=\sigma_s(3\ \mathrm{eV})
\left(\frac{E}{3\ \mathrm{eV}}\right)^{-0.2}.
```

这一关系表示能量升高时总截面缓慢减小。它是模型参数化，不应被理解为对所有能量范围都精确成立。

Lillis 等人列出的参考值为 O 撞击 CO₂、O、N₂ 和 CO 时分别采用 $2.0\times10^{-14}$、$0.6\times10^{-14}$、$1.8\times10^{-14}$ 和 $1.8\times10^{-14}\ \mathrm{cm^2}$。当前 O 截面的 $6.4\times10^{-15}\ \mathrm{cm^2}$ 与其量级一致。

## 3. MGITM 中实际有哪些中性成分

当前放入项目的 MGITM 文件包含 CO₂、O、N₂、CO 和 O₂。它不包含 Ar。

因此，即使截面配置中保留了 Ar，只要 $n_{\mathrm{Ar}}(z)=0$，Ar 对总碰撞系数的贡献就是零：

```math
n_{\mathrm{Ar}}(z)\sigma_{\mathrm{Ar}}(E)=0.
```

Ar 不会被 Monte Carlo 模型抽中。保留 Ar 配置的意义是以后可以使用包含 Ar 密度的其他大气输入，而不是假设当前 MGITM 已经提供 Ar。

## 4. 如何决定这一步是否碰撞

设粒子在当前位置走一步 $ds$。当这一步足够短时，发生碰撞的概率近似为

```math
P_{\mathrm{coll}}\approx ds\,\kappa_{\mathrm{tot}}
=ds\sum_s n_s\sigma_s.
```

代码抽取一个 $0$ 到 $1$ 的均匀随机数 $r$。如果

```math
r<P_{\mathrm{coll}},
```

就认为这一步发生碰撞。

Ali 给出的步长规则为：

* 如果 $\lambda<10\ \mathrm{km}$，使用 $ds=0.1\lambda$。
* 如果 $\lambda>10\ \mathrm{km}$，使用 $ds=1\ \mathrm{km}$。

在第一种情况中，$P_{\mathrm{coll}}\approx0.1$。这能避免一步太长而漏掉多个可能的碰撞。

如果需要对任意步长使用更严格的表达式，可以写成

```math
P_{\mathrm{coll}}=1-\exp(-ds\,\kappa_{\mathrm{tot}}).
```

当 $ds\,\kappa_{\mathrm{tot}}\ll1$ 时，它与 $ds\,\kappa_{\mathrm{tot}}$ 的线性近似相同。

## 5. 碰撞发生后，如何选择目标成分

已知这一步发生了碰撞，与成分 $s$ 碰撞的条件概率为

```math
P(s\mid\mathrm{coll})=
\frac{n_s(z)\sigma_s(E)}
{\sum_j n_j(z)\sigma_j(E)}.
```

例如，某高度的 CO₂ 密度远高于其他成分，即使各成分截面相近，大多数碰撞仍然会发生在 CO₂ 上。较高高度的 O 相对丰度上升后，O 碰撞会变得更重要。

这里仅决定“与谁碰撞”。碰撞以后向什么方向散射，以及损失多少能量，由微分碰撞截面和两体碰撞运动学决定，见 [LAB、COM、散射角和碰撞能量损失](HOT_O_SCATTERING_TWO_BODY_ZH.md)。

## 6. 角度截断对有效碰撞频率的影响

Ali 模型只显式追踪 COM 散射角 $\theta\ge10^\circ$ 的碰撞。小于 $10^\circ$ 的事件非常多，但每次只产生很小的方向和能量变化。

如果总截面包含所有角度，而模拟只保留 $\theta\ge\theta_{\min}$，则有效碰撞系数需要乘以保留比例

```math
f_{\mathrm{cut}}=
\frac{\displaystyle\int_{\theta_{\min}}^\pi
\frac{d\sigma}{d\Omega}\sin\theta\,d\theta}
{\displaystyle\int_0^\pi
\frac{d\sigma}{d\Omega}\sin\theta\,d\theta}.
```

对于当前使用的角分布

```math
\frac{d\sigma}{d\Omega}\propto
\sin^\beta\left(\frac{\theta}{2}\right),
\qquad \beta=-1.85,
```

有

```math
f_{\mathrm{cut}}
=1-\sin^{\beta+2}\left(\frac{\theta_{\min}}{2}\right).
```

当 $\theta_{\min}=10^\circ$ 时，保留比例约为 $0.307$。也就是说，约 30.7% 的总截面对应模型显式处理的散射角范围。

## 7. 当前输入的限制

1. 截面的能量依赖目前是统一的幂律参数化。
2. 不同目标成分目前使用共同形式的散射角分布。
3. MGITM 顶部约为 251 km。更高处的中性密度由代码进行对数线性外推，这是传输至外逸层时的重要不确定性。
4. 当前模型主要处理弹性碰撞，没有完整加入激发、电离、电荷交换等非弹性通道。

因此，现有配置适合建立可重复的基线模型。若要得到高精度逃逸率，需要进一步引入每种碰撞对的能量相关总截面和微分截面数据。
