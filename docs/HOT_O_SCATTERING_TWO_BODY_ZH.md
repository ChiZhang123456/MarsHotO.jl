# COM 散射角与两体碰撞

## 1. 当前角度分布约定

MarsHotO 读取
`data/cross_sections/scattering_angle_distribution.txt`。这个查找表来自
MarsASPEN，并由 Kallio 与 Barabash (2001) 的图数字化得到。

原始查找表把第二列标记为 H 或 H ENA 的 LAB 散射角。当前 MarsHotO 按照模型
约定，将这些角度数值直接作为热 O 碰撞的经验 COM 散射角使用。这是一项明确
的经验近似，不表示 Kallio 原文把该角度定义为 COM。

表格是一个逆累积分布：

```math
U\sim\mathcal U(0,1),
```

```math
\theta_{\mathrm{COM}}=F^{-1}(U).
```

程序在表格相邻点之间进行线性插值。表中角度约为 $0.12^\circ$ 到
$180^\circ$。不设置人为的最小角度截断。方位角独立抽取：

```math
\phi=2\pi U_\phi,\qquad U_\phi\sim\mathcal U(0,1).
```

## 2. LAB 和 COM 分别是什么

LAB 是相对火星静止的坐标系。热 O 的位置、速度、重力轨迹和最终输出都在
LAB 中描述。

COM 是两个碰撞粒子的质心参考系。碰撞前的质心速度为

```math
\mathbf V_{\mathrm{COM}}
=
\frac{m_1\mathbf v_1+m_2\mathbf v_2}
{m_1+m_2}.
```

当前模型采用静止靶近似：

```math
\mathbf v_2=0.
```

相对速度为

```math
\mathbf g=\mathbf v_1-\mathbf v_2.
```

模型抽取的 $\theta_{\mathrm{COM}}$ 是碰撞前后相对速度方向的夹角：

```math
\cos\theta_{\mathrm{COM}}
=
\frac{\mathbf g\cdot\mathbf g'}
{|\mathbf g||\mathbf g'|}.
```

## 3. Ali Rahmati 的能量损失公式

Rahmati 博士论文第 2.2.4 节公式 2.19 给出入射热 O 的相对能量损失：

```math
\boxed{
\frac{\Delta E}{E}
=
\frac{2m_1m_2}{(m_1+m_2)^2}
\left(1-\cos\theta_{\mathrm{COM}}\right)
}
```

其中：

* $m_1$ 是入射热 O 的质量。
* $m_2$ 是靶中性粒子的质量。
* $\theta_{\mathrm{COM}}$ 是 COM 散射角。

当前代码把查找表抽到的角度直接代入这个公式，不再把它当作 LAB 角，也不进行
LAB 到 COM 的角度转换。

对于 O 撞 O，$m_1=m_2$，所以

```math
\frac{\Delta E}{E}
=
\frac12\left(1-\cos\theta_{\mathrm{COM}}\right)
=
\sin^2\left(\frac{\theta_{\mathrm{COM}}}{2}\right).
```

因此：

| COM 角 | $\Delta E/E$ |
|---:|---:|
| $0^\circ$ | 0 |
| $60^\circ$ | 0.25 |
| $90^\circ$ | 0.50 |
| $120^\circ$ | 0.75 |
| $180^\circ$ | 1.00 |

这与把同一数值解释成 LAB 角完全不同。现在 O 撞 O 的角度可以覆盖整个
$0^\circ$ 到 $180^\circ$ COM 范围。

## 4. 如何计算碰撞后的两个速度

首先围绕碰撞前相对速度方向旋转 $\theta_{\mathrm{COM}}$ 和 $\phi$，得到
$\mathbf g'$。弹性碰撞满足

```math
|\mathbf g'|=|\mathbf g|.
```

然后转换回 LAB：

```math
\mathbf v_1'
=
\mathbf V_{\mathrm{COM}}
+
\frac{m_2}{m_1+m_2}\mathbf g',
```

```math
\mathbf v_2'
=
\mathbf V_{\mathrm{COM}}
-
\frac{m_1}{m_1+m_2}\mathbf g'.
```

这些公式同时保证动量守恒

```math
m_1\mathbf v_1+m_2\mathbf v_2
=
m_1\mathbf v_1'+m_2\mathbf v_2',
```

以及总动能守恒

```math
\frac12m_1v_1^2+\frac12m_2v_2^2
=
\frac12m_1v_1'^2+\frac12m_2v_2'^2.
```

热 O 损失的动能全部转移给靶粒子。

## 5. 次级热 O

如果靶粒子是 O，碰撞后的 $\mathbf v_2'$ 属于另一个 O 原子。如果它的能量
高于最低追踪能量，程序把它作为次级热 O 加入待追踪队列。次级 O 继承相同的
宏粒子权重。

## 6. 一次碰撞的计算顺序

1. 用中性密度和总截面判断是否发生碰撞。
2. 按 $n_s\sigma_s$ 选择靶成分。
3. 从 Kallio 查找表抽取一个数值，并将其解释为
   $\theta_{\mathrm{COM}}$。
4. 独立抽取方位角 $\phi$。
5. 在 COM 中旋转相对速度。
6. 转换回 LAB，得到入射粒子和靶粒子的碰撞后速度。
7. 在 LAB 中继续追踪入射热 O。
8. 如果靶粒子是 O，则按需要继续追踪次级热 O。

## 7. 模型近似

Kallio 与 Barabash (2001) 的原始分布描述高能 H 或 H ENA，并且原表标记为
LAB 角。当前模型把相同数值作为热 O 的 COM 角使用，目的是与 Rahmati 的
COM 能量损失公式和 COM 两体运动学直接结合。

未来获得 O 与 CO2、O、N2、CO 和 O2 各自的微分截面后，应使用专属的
COM 角分布替换当前经验分布。

## 参考文献

Kallio, E., and Barabash, S. (2001), Atmospheric effects of precipitating
energetic hydrogen atoms on the Martian atmosphere, Journal of Geophysical
Research: Space Physics, 106(A1), 165 to 177,
https://doi.org/10.1029/2000JA002003.

Rahmati, A. (2016), Oxygen Exosphere of Mars: Evidence from Pickup Ions
Measured by MAVEN, PhD dissertation, University of Kansas, Section 2.2.4.
