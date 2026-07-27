# MarsHotO 中热 O 两体碰撞的具体过程

## 1. 碰撞模型要解决什么问题

O₂⁺ 解离复合反应在火星热层中产生热 O。一个新产生的热 O 原子在运动过程中会与背景 O、CO、N₂、O₂、Ar 和 CO₂ 碰撞。每次碰撞都会改变热 O 的方向和动能，进而影响它最终向上逃逸、向下沉降，还是形成弹道或卫星轨道。

MarsHotO 使用 Monte Carlo 方法逐个追踪热 O：

```text
产生热 O
    ↓
沿轨迹运动并受到火星重力作用
    ↓
累计碰撞光学深度
    ↓
确定碰撞位置和目标成分
    ↓
抽取 COM 散射角和方位角
    ↓
完成两体弹性碰撞
    ↓
将速度转换回 LAB 坐标系
    ↓
继续追踪
```

本文重点解释从“判断碰撞发生”到“得到碰撞后速度”的完整过程。

## 2. 一次碰撞中的两个粒子

入射粒子是热 O，记为粒子 1：

$$
m_1=m_{\mathrm O},
\qquad
\mathbf v_1=\text{碰撞前热 O 的 LAB 速度}.
$$

目标粒子是背景中性成分，记为粒子 2：

$$
m_2=m_s,
\qquad
\mathbf v_2=\text{碰撞前目标粒子的 LAB 速度}.
$$

Ali Rahmati 的解析能量损失公式假设目标粒子静止：

$$
\mathbf v_2=0.
$$

MarsHotO 的完整碰撞运动学允许 $\mathbf v_2\ne0$。代码根据 MGITM 中性温度 $T_n$ 抽取目标粒子的 Maxwell 速度：

$$
v_{2x},v_{2y},v_{2z}
\sim
\mathcal N\left(0,\frac{k_{\mathrm B}T_n}{m_2}\right).
$$

## 3. LAB 和 COM 坐标系

### 3.1 LAB 坐标系

LAB 是火星大气所在的参考系。MGITM 密度和温度、火星重力、粒子位置、向上运动和向下运动都在 LAB 中定义。

### 3.2 COM 坐标系

COM 是质心坐标系。质心速度为：

$$
\boxed{
\mathbf V_{\mathrm{COM}}
=
\frac{m_1\mathbf v_1+m_2\mathbf v_2}{m_1+m_2}
}
$$

定义碰撞前相对速度：

$$
\boxed{
\mathbf g=\mathbf v_1-\mathbf v_2
}
$$

两个粒子在 COM 中的速度为：

$$
\mathbf v_{1,\mathrm{COM}}
=
\frac{m_2}{m_1+m_2}\mathbf g,
$$

$$
\mathbf v_{2,\mathrm{COM}}
=
-\frac{m_1}{m_1+m_2}\mathbf g.
$$

因此：

$$
m_1\mathbf v_{1,\mathrm{COM}}
+
m_2\mathbf v_{2,\mathrm{COM}}
=0.
$$

也就是说，在 COM 中两个粒子的动量大小相等、方向相反。这使弹性碰撞可以简化为对相对速度方向的一次旋转。

## 4. 如何确定下一次碰撞的位置

### 4.1 能量相关总截面

对背景成分 $s$，Ali 使用：

$$
\boxed{
\sigma_s(E)
=
\sigma_s(3\,\mathrm{eV})
\left(\frac{E}{3\,\mathrm{eV}}\right)^{-0.2}
}
$$

当前代码使用热 O 的 LAB 动能：

$$
E=\frac{1}{2}m_{\mathrm O}|\mathbf v_1|^2.
$$

如果目标静止，这就是相对碰撞能量。当前代码会为目标抽取热速度，但是计算总截面时仍采用热 O 的 LAB 动能。这与 Ali 的目标静止近似一致，是当前模型的一项明确近似。

### 4.2 局地碰撞系数

MGITM 给出位置 $\mathbf r$ 处各背景成分的密度 $n_s(\mathbf r)$。考虑最小散射角后，局地总碰撞系数为：

$$
\boxed{
\kappa(\mathbf r,E)
=
\sum_s n_s(\mathbf r)\sigma_s(E)f_\sigma
}
$$

$\kappa$ 的单位是 $\mathrm{m^{-1}}$。如果局地状态不变，平均自由程为：

$$
\lambda=\frac{1}{\kappa}.
$$

### 4.3 光学深度抽样

模型先产生均匀随机数 $u$，然后抽取下一次碰撞所需的光学深度：

$$
\boxed{
\tau_*=-\ln u,
\qquad
u\sim U(0,1)
}
$$

粒子每前进一段 $\Delta s$，累计：

$$
\tau\leftarrow\tau+\kappa\Delta s.
$$

当：

$$
\tau\ge\tau_*
$$

时，在当前位置执行碰撞。这个方法可以自然处理大气密度和碰撞截面沿轨迹变化的情况。

## 5. 如何选择目标成分

如果当前位置已经确定发生碰撞，目标成分 $s$ 的权重为：

$$
w_s=n_s\sigma_s(E)f_\sigma.
$$

目标成分的条件概率为：

$$
\boxed{
P(s\mid\mathrm{collision})
=
\frac{w_s}{\sum_j w_j}
}
$$

因此，目标选择同时取决于局地密度和碰撞截面。CO₂ 的截面较大并不意味着每次都与 CO₂ 碰撞。如果当前位置 O 密度远高于 CO₂，O 仍可能是最常被选中的目标。

## 6. 从微分截面得到散射角概率

Ali 使用 Kharchenko et al. (2000) 的 O 与 O 微分截面拟合：

$$
\boxed{
\frac{d\sigma}{d\Omega}
=
\alpha
\sin^\beta\left(\frac{\theta}{2}\right)
}
$$

其中：

$$
\alpha=0.36\times10^{-16}
\ \mathrm{cm^2\,sr^{-1}},
\qquad
\beta=-1.85.
$$

$\theta$ 是 COM 散射角。它不是热 O 最终在 LAB 中转过的角度。

### 6.1 为什么不能直接把微分截面当作角度概率

球坐标的立体角元为：

$$
d\Omega=\sin\theta\,d\theta\,d\phi.
$$

对方位角积分后：

$$
d\Omega=2\pi\sin\theta\,d\theta.
$$

所以散射角的概率必须满足：

$$
p(\theta)
\propto
\frac{d\sigma}{d\Omega}\sin\theta.
$$

$\sin\theta$ 是立体角 Jacobian。如果漏掉这一项，得到的角度分布会错误。

### 6.2 归一化概率密度

在：

$$
\theta_{\min}\le\theta\le\pi
$$

范围内，归一化概率密度为：

$$
\boxed{
p(\theta)
=
\frac{
(\beta+2)
\sin^{\beta+1}(\theta/2)
\cos(\theta/2)
}{
2\left[
1-\sin^{\beta+2}(\theta_{\min}/2)
\right]
}
}
$$

它满足：

$$
\int_{\theta_{\min}}^\pi p(\theta)\,d\theta=1.
$$

### 6.3 反 CDF 抽样

令：

$$
q=\beta+2=0.15.
$$

首先抽取：

$$
u_\theta\sim U(0,1).
$$

然后计算：

$$
A=\sin^q\left(\frac{\theta_{\min}}{2}\right),
$$

$$
x=A+u_\theta(1-A),
$$

$$
\boxed{
\theta
=
2\arcsin\left(x^{1/q}\right)
}
$$

方位角没有优先方向，因此：

$$
\boxed{
\phi=2\pi u_\phi,
\qquad
u_\phi\sim U(0,1)
}
$$

## 7. 10°截断及其物理含义

微分截面在小角度处很大，因此会产生大量几乎不改变方向和能量的碰撞。Ali 指出，小于约 $10^\circ$ 的碰撞对总体热 O 通量和逃逸率影响较小。MarsHotO 当前默认：

$$
\theta_{\min}=10^\circ.
$$

截断散射角时，总碰撞截面也必须按相同比例减小。保留下来的角积分比例为：

$$
\boxed{
f_\sigma(\theta_{\min})
=
1-\sin^{\beta+2}
\left(\frac{\theta_{\min}}{2}\right)
}
$$

当 $\theta_{\min}=10^\circ$ 时：

$$
f_\sigma\approx0.307.
$$

也就是说，采用 Ali 拟合函数时，$10^\circ$ 以上的碰撞约占拟合总截面的 30.7%。代码对碰撞频率和角度抽样使用同一个 $\theta_{\min}$，保证两部分一致。

## 8. 两体弹性碰撞的核心计算

### 8.1 构造相对速度坐标基

沿碰撞前相对速度定义：

$$
\mathbf e_0=\frac{\mathbf g}{|\mathbf g|}.
$$

再构造两个与 $\mathbf e_0$ 垂直并且彼此垂直的单位向量：

$$
\mathbf e_1,\qquad\mathbf e_2.
$$

$(\mathbf e_0,\mathbf e_1,\mathbf e_2)$ 构成局地正交基。

### 8.2 在 COM 中旋转相对速度

弹性碰撞不改变 COM 相对速度的大小：

$$
|\mathbf g'|=|\mathbf g|.
$$

抽到 $\theta$ 和 $\phi$ 后：

$$
\boxed{
\mathbf g'
=
|\mathbf g|
\left[
\cos\theta\,\mathbf e_0
+
\sin\theta
\left(
\cos\phi\,\mathbf e_1
+
\sin\phi\,\mathbf e_2
\right)
\right]
}
$$

可以把它理解为，把原本沿 $\mathbf e_0$ 的相对速度箭头转过 $\theta$，再用 $\phi$ 决定绕原方向转到哪个方位。

### 8.3 计算碰撞后的 COM 速度

$$
\mathbf v'_{1,\mathrm{COM}}
=
\frac{m_2}{m_1+m_2}\mathbf g',
$$

$$
\mathbf v'_{2,\mathrm{COM}}
=
-\frac{m_1}{m_1+m_2}\mathbf g'.
$$

### 8.4 转换回 LAB

给两个 COM 速度重新加上质心速度：

$$
\boxed{
\mathbf v'_1
=
\mathbf V_{\mathrm{COM}}
+
\frac{m_2}{m_1+m_2}\mathbf g'
}
$$

$$
\boxed{
\mathbf v'_2
=
\mathbf V_{\mathrm{COM}}
-
\frac{m_1}{m_1+m_2}\mathbf g'
}
$$

$\mathbf v'_1$ 是热 O 碰撞后的 LAB 速度。$\mathbf v'_2$ 是目标粒子碰撞后的 LAB 速度。

当前主输运程序继续追踪粒子 1。对于 O 与 O 碰撞，反冲目标 O 也可能成为次级热 O。碰撞函数已经返回 $\mathbf v'_2$，但是主输运程序目前还没有将它加入粒子队列。

## 9. 为什么动量和能量自动守恒

### 9.1 动量守恒

把碰撞后的 LAB 速度代入总动量：

$$
m_1\mathbf v'_1+m_2\mathbf v'_2.
$$

与 $\mathbf g'$ 有关的两项会互相抵消，因此：

$$
m_1\mathbf v'_1+m_2\mathbf v'_2
=
(m_1+m_2)\mathbf V_{\mathrm{COM}}
=
m_1\mathbf v_1+m_2\mathbf v_2.
$$

所以线性动量严格守恒。

### 9.2 动能守恒

COM 中的相对动能为：

$$
E_{\mathrm{COM}}
=
\frac{1}{2}\mu|\mathbf g|^2,
$$

其中约化质量为：

$$
\mu=\frac{m_1m_2}{m_1+m_2}.
$$

碰撞只旋转 $\mathbf g$，不改变它的大小，所以 COM 动能不变。质心速度也不变，因此 LAB 中两个粒子的总动能严格守恒。

热 O 的“能量损失”不是总能量消失，而是热 O 的一部分动能转移给了目标粒子。

## 10. 热 O 能量损失公式

当目标粒子在 LAB 中初始静止时：

$$
\boxed{
\frac{\Delta E}{E}
=
\frac{2m_1m_2}{(m_1+m_2)^2}
(1-\cos\theta)
}
$$

其中：

$$
\Delta E=E_1-E'_1.
$$

这个公式说明，能量损失由质量比和 COM 散射角共同决定。

### 10.1 小角度碰撞

当 $\theta$ 很小时：

$$
1-\cos\theta\approx\frac{\theta^2}{2}.
$$

所以：

$$
\frac{\Delta E}{E}\propto\theta^2.
$$

这解释了为什么小角度前向散射的数量很多，但是单次能量传递很小。

### 10.2 完全后向散射

当 $\theta=180^\circ$ 时：

$$
\left(\frac{\Delta E}{E}\right)_{\max}
=
\frac{4m_1m_2}{(m_1+m_2)^2}.
$$

如果两个粒子质量相等：

$$
m_1=m_2,
$$

则：

$$
\left(\frac{\Delta E}{E}\right)_{\max}=1.
$$

这表示热 O 与静止 O 发生完全后向的等质量碰撞时，入射 O 可以把全部动能传递给目标 O。

## 11. 一个具体数值例子

假设热 O 与静止 CO₂ 碰撞：

$$
m_1=16\ \mathrm{amu},
\qquad
m_2=44\ \mathrm{amu}.
$$

假设抽到的 COM 散射角为：

$$
\theta=60^\circ.
$$

质量因子为：

$$
\frac{2m_1m_2}{(m_1+m_2)^2}
=
\frac{2\times16\times44}{(16+44)^2}
\approx0.3911.
$$

又因为：

$$
1-\cos60^\circ=0.5,
$$

所以：

$$
\frac{\Delta E}{E}
\approx0.3911\times0.5
\approx0.1956.
$$

如果碰撞前热 O 的能量为 3 eV：

$$
\Delta E\approx0.587\ \mathrm{eV},
$$

$$
E'\approx2.413\ \mathrm{eV}.
$$

作为对比，如果目标是 O，则 $m_1=m_2=16\ \mathrm{amu}$。在相同的 $60^\circ$ COM 散射角下：

$$
\frac{\Delta E}{E}
=
\frac{1}{2}(1-\cos60^\circ)
=0.25.
$$

等质量碰撞的能量传递更有效。代码为目标抽取热速度后，最终结果会在静止目标解析结果附近产生小幅变化。

## 12. 一次碰撞在程序中的完整顺序

当前 `transport_particle!` 按以下顺序执行。

1. 根据热 O 位置计算高度。
2. 从 MGITM 剖面插值得到 $T_n$ 和各中性成分密度。
3. 根据热 O 的 LAB 速度计算动能。
4. 计算各成分的能量相关总截面。
5. 应用与 $\theta_{\min}$ 对应的有效截面比例。
6. 计算局地总碰撞系数 $\kappa$。
7. 粒子在火星重力作用下前进一步。
8. 累计光学深度 $\kappa\Delta s$。
9. 达到随机阈值后，按 $n_s\sigma_s$ 选择目标成分。
10. 根据 $T_n$ 和目标质量抽取目标 Maxwell 速度。
11. 从解析反 CDF 抽取 $\theta_{\mathrm{COM}}$。
12. 均匀抽取方位角 $\phi$。
13. 计算质心速度和碰撞前相对速度。
14. 在 COM 中旋转相对速度。
15. 将两个粒子的速度转换回 LAB。
16. 使用 $\mathbf v'_1$ 更新热 O。
17. 重置碰撞光学深度，继续追踪下一次碰撞。

## 13. 公式和 Julia 文件的对应关系

| 物理内容 | Julia 文件 | 主要函数 |
|---|---|---|
| 总截面随能量变化 | `src/cross_sections.jl` | `total_cross_section` |
| 局地碰撞系数 | `src/cross_sections.jl` | `collision_coefficient` |
| 目标成分选择 | `src/cross_sections.jl` | `choose_collision_target` |
| 微分截面 | `src/scattering.jl` | `differential_cross_section` |
| 有效角积分比例 | `src/scattering.jl` | `angular_cross_section_fraction` |
| 角度 PDF 和 CDF | `src/scattering.jl` | `scattering_angle_pdf`, `scattering_angle_cdf` |
| COM 散射角抽样 | `src/scattering.jl` | `sample_scattering_angle` |
| 两体碰撞速度变换 | `src/collision_kinematics.jl` | `elastic_collision` |
| 解析能量损失 | `src/collision_kinematics.jl` | `fractional_energy_loss` |
| 粒子输运循环 | `src/transport.jl` | `transport_particle!` |

## 14. 当前模型的主要假设和限制

1. 所有背景成分的总截面都采用 $E^{-0.2}$ 能量依赖。
2. 所有背景成分使用相同的 COM 角分布，该角分布来自 O 与 O 碰撞拟合。
3. 不同成分通过局地密度、总截面和目标质量产生不同结果。
4. 默认最小散射角为 $10^\circ$。
5. 截面能量目前使用热 O 的 LAB 动能。
6. 碰撞运动学允许目标粒子具有 Maxwell 热速度。
7. 主输运目前只继续追踪入射热 O。
8. O 与 O 碰撞产生的次级热 O 尚未加入粒子队列。
9. 当前碰撞为弹性碰撞，不包括激发、电离或化学反应导致的非弹性能量损失。

## 15. 最直观的理解

可以把 COM 两体碰撞理解为：

1. 先和两个粒子的质心一起运动。
2. 在这个参考系中，两个粒子的动量大小相等、方向相反。
3. 碰撞只把相对速度方向转过一个随机角度，速度大小不变。
4. 再把质心的整体运动加回来。
5. 由于目标质量不同，变回 LAB 后热 O 分到的速度和能量也不同。

Monte Carlo 并不是直接随机指定能量损失。模型随机抽取碰撞位置、碰撞对象、COM 散射角和方位角。能量损失随后由粒子质量、初始速度和两体守恒关系自动确定。

## 16. 参考资料

1. Rahmati, A. (2016), *Oxygen Exosphere of Mars: Evidence from Pickup Ions Measured by MAVEN*, PhD dissertation, University of Kansas, Sections 2.2.3 to 2.2.5.
2. Kharchenko, V., Dalgarno, A., Zygelman, B., and Yee, J. H. (2000), Energy transfer in collisions of oxygen atoms in the terrestrial atmosphere.
3. Fox, J. L., and Hać, A. (2014), collision cross sections used for hot oxygen transport calculations.
