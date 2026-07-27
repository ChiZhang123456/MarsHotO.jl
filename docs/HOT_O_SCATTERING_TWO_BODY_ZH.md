# LAB、COM、散射角和碰撞能量损失

## 1. 先说最重要的结论

散射角描述的是一次碰撞把两粒子的相对运动方向转了多少。MarsHotO 先在 COM 系中抽取散射角，再利用动量守恒和能量守恒，计算两个粒子碰撞后的速度，最后转换回 LAB 系继续追踪。

散射角越小，热 O 通常只轻微改变方向，能量损失很小。散射角越大，传给大气粒子的能量通常越多。

## 2. LAB 是什么

LAB 是相对于火星静止的坐标系。MGITM 的高度、热 O 的上升和下降、火星重力、逃逸速度以及最后画出的能量，都在 LAB 系中描述。

在一次碰撞前：

* 热 O 的 LAB 速度是 $\mathbf v_1$。
* 目标中性粒子的 LAB 速度是 $\mathbf v_2$。

Rahmati 模型采用静止靶近似，即 $\mathbf v_2=0$。

## 3. COM 到底是什么

COM 是 center of mass frame 的缩写，中文称为质心参考系。

两个粒子的质心速度为

```math
\mathbf V_{\mathrm{COM}}=
\frac{m_1\mathbf v_1+m_2\mathbf v_2}
{m_1+m_2}.
```

COM 系就是以 $\mathbf V_{\mathrm{COM}}$ 跟随这对粒子一起运动的参考系。在 COM 系中，两个粒子的总动量为零。它们总是沿相反方向运动：

```math
m_1\mathbf u_1+m_2\mathbf u_2=0,
```

其中

```math
\mathbf u_1=\mathbf v_1-\mathbf V_{\mathrm{COM}},
\qquad
\mathbf u_2=\mathbf v_2-\mathbf V_{\mathrm{COM}}.
```

COM 不是另一个物理地点，只是换了一个运动的观察者。这样做的好处是，弹性碰撞只需要改变相对速度的方向，不需要改变它的大小。

## 4. 一个直观例子

假设一颗 O 原子向右撞上一颗静止的 CO₂ 分子。

* 站在火星上看，这是 LAB 系。O 向右运动，CO₂ 起初静止。
* 跟着二者质心一起向右移动，这是 COM 系。在这个观察者看来，O 和 CO₂ 在碰撞前从相反方向靠近。
* 碰撞后，它们仍然从质心向相反方向离开，但离开的轴线可能转了一个角度。

这条轴线转过的角度，就是模型使用的 COM 散射角 $\theta$。

## 5. 散射角究竟是哪两个方向之间的夹角

碰撞前的相对速度为

```math
\mathbf g=\mathbf v_1-\mathbf v_2.
```

碰撞后的相对速度为

```math
\mathbf g'=\mathbf v_1'-\mathbf v_2'.
```

COM 散射角定义为

```math
\cos\theta=
\frac{\mathbf g\cdot\mathbf g'}
{|\mathbf g||\mathbf g'|}.
```

所以，$\theta$ 是碰撞前后相对速度方向之间的夹角。

它不一定等于热 O 在 LAB 系中的偏转角。热 O 的 LAB 偏转角应另外计算：

```math
\cos\Theta_{\mathrm{LAB}}=
\frac{\mathbf v_1\cdot\mathbf v_1'}
{|\mathbf v_1||\mathbf v_1'|}.
```

两者不同的原因是，从 COM 转回 LAB 时还要加上 $\mathbf V_{\mathrm{COM}}$。

## 6. 微分碰撞截面如何给出散射角概率

总碰撞截面回答“碰撞是否发生”。微分碰撞截面回答“发生碰撞后往哪个方向散射”。

当前角分布采用

```math
\frac{d\sigma}{d\Omega}
=\alpha
\sin^\beta\left(\frac{\theta}{2}\right),
\qquad \beta=-1.85.
```

因为 $\beta$ 为负数，当 $\theta$ 很小时，微分截面很大。这表示碰撞强烈偏向前向散射。

但是，抽取 $\theta$ 时不能直接把 $d\sigma/d\Omega$ 当成一维概率。球面上一圈角度对应的立体角为

```math
d\Omega=\sin\theta\,d\theta\,d\phi.
```

因此散射角的一维概率密度满足

```math
p(\theta)\propto
\frac{d\sigma}{d\Omega}\sin\theta.
```

$\sin\theta$ 是从球面立体角转换成角度概率时必须包含的 Jacobian。

## 7. 当前分布如何抽样

对 $\theta\ge\theta_{\min}$ 归一化后，当前分布可写为

```math
p(\theta)=
\frac{(\beta+2)
\sin^{\beta+1}(\theta/2)\cos(\theta/2)}
{2\left[1-\sin^{\beta+2}(\theta_{\min}/2)\right]}.
```

代码使用逆累积分布抽样。取均匀随机数 $U\in[0,1]$：

```math
\theta=
2\arcsin
\left[
s_{\min}^{\beta+2}
+U\left(1-s_{\min}^{\beta+2}\right)
\right]^{1/(\beta+2)},
```

其中

```math
s_{\min}=\sin\left(\frac{\theta_{\min}}{2}\right).
```

方位角没有偏好，因此独立抽取

```math
\phi=2\pi U_\phi.
```

当前采用 $\theta_{\min}=10^\circ$。被省略的更小角度事件很多，但单次能量变化很小。这个截断是计算效率和小角累积效应之间的近似。

## 8. 如何把抽到的角度变成新方向

先令

```math
\hat{\mathbf e}_0=\frac{\mathbf g}{|\mathbf g|}.
```

再构造两个与 $\hat{\mathbf e}_0$ 垂直、彼此也垂直的单位向量 $\hat{\mathbf e}_1$ 和 $\hat{\mathbf e}_2$。

抽到 $\theta$ 和 $\phi$ 后，新的相对速度为

```math
\mathbf g'
=|\mathbf g|
\left[
\cos\theta\,\hat{\mathbf e}_0
+\sin\theta
\left(
\cos\phi\,\hat{\mathbf e}_1
+\sin\phi\,\hat{\mathbf e}_2
\right)
\right].
```

弹性碰撞中 $|\mathbf g'|=|\mathbf g|$，只有方向发生改变。

## 9. 如何转换回 LAB 速度

两个粒子碰撞后的 LAB 速度为

```math
\mathbf v_1'
=\mathbf V_{\mathrm{COM}}
+\frac{m_2}{m_1+m_2}\mathbf g',
```

```math
\mathbf v_2'
=\mathbf V_{\mathrm{COM}}
-\frac{m_1}{m_1+m_2}\mathbf g'.
```

这两个式子同时保证线动量守恒：

```math
m_1\mathbf v_1+m_2\mathbf v_2
=m_1\mathbf v_1'+m_2\mathbf v_2',
```

以及弹性碰撞的总动能守恒：

```math
\frac12m_1v_1^2+\frac12m_2v_2^2
=
\frac12m_1v_1'^2+\frac12m_2v_2'^2.
```

注意，总动能守恒不表示热 O 自己的动能不变。热 O 损失的能量转移给了目标粒子。

## 10. 散射角为什么决定能量损失

当目标粒子在 LAB 系中静止时，入射热 O 的相对能量损失为

```math
\frac{\Delta E}{E}
=
\frac{2m_1m_2}{(m_1+m_2)^2}
(1-\cos\theta).
```

这个式子包含两个因素。

第一，角度因素：

```math
1-\cos\theta.
```

* $\theta=0^\circ$ 时，它等于 0，热 O 不损失能量。
* $\theta=90^\circ$ 时，它等于 1。
* $\theta=180^\circ$ 时，它等于 2，能量传递最大。

第二，质量因素：

```math
\frac{2m_1m_2}{(m_1+m_2)^2}.
```

两个粒子质量越接近，最大能量传递越有效。质量差异很大时，入射粒子更难把大部分能量交给目标粒子。

## 11. 数值例子

### 11.1 O 撞击 CO₂

取 $m_1=16$、$m_2=44$，若 $\theta=60^\circ$：

```math
\frac{\Delta E}{E}
=
\frac{2\times16\times44}{(16+44)^2}
(1-\cos60^\circ)
\approx0.1956.
```

一个 3 eV 的热 O 碰撞后剩余能量约为

```math
E'=3(1-0.1956)\approx2.41\ \mathrm{eV}.
```

### 11.2 O 撞击 O

两者质量相同，取 $\theta=60^\circ$：

```math
\frac{\Delta E}{E}
=\frac12(1-\cos60^\circ)=0.25.
```

热 O 把 25% 的初始能量传给目标 O。

如果发生 $180^\circ$ 的 COM 后向散射：

```math
\frac{\Delta E}{E}=1.
```

在静止靶和理想弹性碰撞条件下，入射 O 可以把全部动能传给目标 O。原来的目标 O 随后成为新的高速 O。

## 12. 为什么要追踪次级 O

当目标成分是 O 时，碰撞后的 $\mathbf v_2'$ 属于另一个 O 原子。如果它的能量仍高于最低追踪能量，它也属于非热 O，必须作为次级粒子继续追踪。

因此一次 O 和 O 碰撞可能把一条粒子轨迹变成两条：

```text
入射热 O
    ↓
与背景 O 碰撞
    ↓
原热 O 继续运动 + 被加速的次级 O 继续运动
```

这不会创造额外能量。两条轨迹的总动量和总能量仍然等于碰撞前的总量。

## 13. 从一次碰撞到完整传输模型

一次 Monte Carlo 碰撞的实际顺序为：

1. 用 LAB 位置和速度沿重力轨迹移动粒子。
2. 用中性密度和总截面决定是否碰撞。
3. 按 $n_s\sigma_s$ 选择目标成分。
4. 抽取 COM 散射角 $\theta$ 和方位角 $\phi$。
5. 旋转相对速度 $\mathbf g$ 得到 $\mathbf g'$。
6. 用两体公式计算 $\mathbf v_1'$ 和 $\mathbf v_2'$。
7. 在 LAB 系中继续追踪主粒子。
8. 如果目标是 O，并且次级 O 能量足够高，则继续追踪次级粒子。

反复执行这一过程，就可以累计热 O 在不同高度、能量和速度区间内的驻留时间，并得到 O 冕的 $n(E,z)$ 和速度分布。

![碰撞截面和散射物理](../examples/figures/hot_o_collision_cross_sections_and_scattering.png)

## 14. 当前模型近似

1. 当前主要处理两体弹性碰撞。
2. Rahmati 模型把目标中性粒子设为静止。
3. 不同成分暂时采用共同的前向峰化散射角分布。
4. 小于 $10^\circ$ 的散射没有逐次显式追踪。
5. 若以后获得每一种碰撞对的能量相关微分截面，应分别替换当前共同角分布。

这些近似不改变 LAB 和 COM 的定义，也不改变两体弹性碰撞的守恒关系。它们主要影响碰撞频率、角度统计和能量沉降速度。
