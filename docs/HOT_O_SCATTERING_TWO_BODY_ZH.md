# LAB 散射角与两体碰撞

## 1. 当前模型采用什么角分布

当前模型直接读取
`data/cross_sections/scattering_angle_distribution.txt`。这个文件来自
MarsASPEN，原始数据由 Kallio 与 Barabash (2001) 的 Figure 2 数字化得到。
表中两列分别是：

1. 0 到 1 之间的累积概率，也可以看成均匀随机数 $U$。
2. 入射粒子在 LAB 中的散射角 $\Theta_{\mathrm{LAB}}$，单位为度。

因此它本身就是一个逆累积分布表。抽样过程为

```math
U\sim\mathcal U(0,1),
```

```math
\Theta_{\mathrm{LAB}}=F^{-1}(U).
```

程序在表格相邻点之间进行线性插值。方位角独立抽取：

```math
\phi=2\pi U_\phi,\qquad U_\phi\sim\mathcal U(0,1).
```

表中最小角约为 $0.12^\circ$，最大角为 $180^\circ$。代码使用表中的完整
范围，不另外截断角度，也不对总碰撞截面乘角度保留比例。

## 2. LAB 是什么

LAB 就是相对火星静止的坐标系。MGITM 的高度、热 O 的速度、火星重力和
最终输出的能量都在这个坐标系中描述。

碰撞前，入射热 O 的速度为 $\mathbf v_1$。Rahmati 传输流程采用静止靶近似，
所以背景中性粒子的碰撞前速度设为零：

```math
\mathbf v_2=0.
```

当前查找表给出的 $\Theta_{\mathrm{LAB}}$ 是 $\mathbf v_1$ 与碰撞后入射粒子
速度 $\mathbf v_1'$ 之间的夹角：

```math
\cos\Theta_{\mathrm{LAB}}
=
\frac{\mathbf v_1\cdot\mathbf v_1'}
{|\mathbf v_1||\mathbf v_1'|}.
```

这里的角度不是 COM 散射角。程序不会先把它当作 COM 角，再转换回 LAB。

## 3. 随机角度如何变成新的速度方向

首先定义入射方向：

```math
\hat{\mathbf e}_0=\frac{\mathbf v_1}{|\mathbf v_1|}.
```

再构造两个与 $\hat{\mathbf e}_0$ 垂直的单位向量
$\hat{\mathbf e}_1$ 和 $\hat{\mathbf e}_2$。碰撞后的入射粒子方向为

```math
\hat{\mathbf e}'
=
\cos\Theta_{\mathrm{LAB}}\,\hat{\mathbf e}_0
+
\sin\Theta_{\mathrm{LAB}}
\left(
\cos\phi\,\hat{\mathbf e}_1
+
\sin\phi\,\hat{\mathbf e}_2
\right).
```

这个步骤只确定方向。碰撞后的速率还要由两体弹性碰撞关系计算。

## 4. LAB 角对应的碰撞后速率

设入射热 O 的质量为 $m_1$，静止靶粒子的质量为 $m_2$，并定义

```math
r=\frac{m_1}{m_2}.
```

对于当前模型中的 O 撞击 O、CO、N2、O2 和 CO2，均有 $m_1\le m_2$。
入射粒子碰撞后的速率比为

```math
\frac{v_1'}{v_1}
=
\frac{
r\cos\Theta_{\mathrm{LAB}}
+
\sqrt{1-r^2\sin^2\Theta_{\mathrm{LAB}}}
}
{1+r}.
```

因此

```math
\mathbf v_1'
=
v_1
\frac{
r\cos\Theta_{\mathrm{LAB}}
+
\sqrt{1-r^2\sin^2\Theta_{\mathrm{LAB}}}
}
{1+r}
\hat{\mathbf e}'.
```

靶粒子的反冲速度由动量守恒得到：

```math
\mathbf v_2'
=
\frac{m_1}{m_2}
\left(\mathbf v_1-\mathbf v_1'\right).
```

代码测试以下两个守恒关系：

```math
m_1\mathbf v_1
=
m_1\mathbf v_1'+m_2\mathbf v_2',
```

```math
\frac12m_1v_1^2
=
\frac12m_1v_1'^2+\frac12m_2v_2'^2.
```

所以入射热 O 损失的能量正好成为靶粒子的反冲能量。

## 5. 能量损失如何由散射角决定

入射粒子的剩余能量比例为

```math
\frac{E_1'}{E_1}
=
\left(\frac{v_1'}{v_1}\right)^2.
```

因此能量损失比例为

```math
\frac{\Delta E_1}{E_1}
=
1-
\left[
\frac{
r\cos\Theta_{\mathrm{LAB}}
+
\sqrt{1-r^2\sin^2\Theta_{\mathrm{LAB}}}
}
{1+r}
\right]^2.
```

小角度散射通常只引起很小的方向改变和能量损失。大角度散射通常会把更多
能量传给靶粒子。不同靶质量对应不同的能量损失曲线。

特别地，当热 O 撞击静止 O 时，$m_1=m_2$。如果靶 O 获得的反冲能量高于
最低追踪能量，它会作为次级热 O 加入粒子队列，并继承相同的宏粒子权重。

## 6. 一次碰撞的完整 Monte Carlo 顺序

1. 用中性密度和总截面判断本步是否发生碰撞。
2. 按 $n_s\sigma_s$ 选择靶成分。
3. 抽取 $U$，从查找表得到 LAB 散射角。
4. 独立抽取方位角 $\phi$。
5. 用 LAB 两体公式计算入射粒子的碰撞后速率。
6. 旋转入射方向，得到 $\mathbf v_1'$。
7. 用动量守恒计算靶粒子的反冲速度 $\mathbf v_2'$。
8. 在 LAB 中继续追踪入射热 O。
9. 如果靶粒子是 O，则按需要继续追踪次级热 O。

## 7. 当前近似与适用范围

Kallio 与 Barabash (2001) 的角分布针对火星大气中的高能 H 或 H ENA。当前
模型按照用户指定，将同一经验逆累积分布用于热 O 与所有中性成分的碰撞。
这是一个明确的模型近似，并不是 O 与各中性成分的专属微分截面。

后续如果获得 O 与 CO2、O、N2、CO 和 O2 各自的能量相关微分截面，应分别
建立角度和能量相关的散射表，替换当前共同分布。

## 参考文献

Kallio, E., and Barabash, S. (2001), Atmospheric effects of precipitating
energetic hydrogen atoms on the Martian atmosphere, Journal of Geophysical
Research: Space Physics, 106(A1), 165 to 177,
https://doi.org/10.1029/2000JA002003.
