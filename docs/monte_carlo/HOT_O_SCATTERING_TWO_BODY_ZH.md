# Rahmati COM 散射角与两体碰撞

## 1. 当前模型采用的散射角分布

MarsHotO 使用 Rahmati 博士论文第 2.2.3 和 2.2.4 节给出的模型。该模型采用
Kharchenko 等人计算的 O 与 O 微分碰撞截面，并使用解析函数拟合：

```math
\frac{d\sigma}{d\Omega}
=
\alpha\sin^\beta\left(\frac{\theta_{\mathrm{COM}}}{2}\right),
```

其中

```math
\alpha=0.36\times10^{-16}\ \mathrm{cm^2\,sr^{-1}},
\qquad
\beta=-1.85.
```

$\theta_{\mathrm{COM}}$ 是质心坐标系中的散射角。当前模型假设所有 O 与中性
成分碰撞都采用相同形式的 COM 角分布。

## 2. 为什么概率中要包含 $\sin\theta$

微分截面是单位立体角内的截面。球坐标中的立体角元为

```math
d\Omega
=
\sin\theta\,d\theta\,d\phi.
```

因此 COM 极角的一维概率密度不是单独的 $d\sigma/d\Omega$，而是

```math
p(\theta)
\propto
\frac{d\sigma}{d\Omega}\sin\theta.
```

代入 Rahmati 拟合并归一化后得到

```math
p(\theta)
=
\frac{\beta+2}{2}
\sin^{\beta+1}\left(\frac{\theta}{2}\right)
\cos\left(\frac{\theta}{2}\right),
\qquad
0\le\theta\le\pi.
```

虽然 $\theta$ 接近零时概率密度很大，但因为 $\beta+2=0.15>0$，从
$0^\circ$ 开始的积分仍然有限。

## 3. 不使用 10° 截断

Rahmati 的历史 two-stream 计算曾采用 $\theta_{\min}=10^\circ$。当前
MarsHotO 按照本项目约定令

```math
\theta_{\min}=0^\circ.
```

因此：

1. 保留所有小于 $10^\circ$ 的前向散射。
2. 散射角范围为 $0^\circ$ 到 $180^\circ$。
3. 碰撞频率使用完整总截面。
4. 不使用角度保留比例缩放总截面。

## 4. CDF 和逆变换抽样

从零度积分到 $\theta$，可以得到累积分布：

```math
F(\theta)
=
\sin^{\beta+2}\left(\frac{\theta}{2}\right).
```

计算机生成均匀随机数

```math
R\sim\mathcal U(0,1).
```

令 $R=F(\theta)$，反解得到

```math
\boxed{
\theta_{\mathrm{COM}}(R)
=
2\arcsin\left(
R^{1/(\beta+2)}
\right)
}
```

因为

```math
\frac{1}{\beta+2}
=
\frac{1}{0.15}
\approx6.667,
```

大部分随机数会产生非常小的 COM 散射角。只有 $R$ 非常接近 1 时，角度才会
迅速接近 $180^\circ$。

如果以后需要加入非零最小角度，一般形式为

```math
\theta(R)
=
2\arcsin
\left[
R\left(1-C_{\min}\right)+C_{\min}
\right]^{1/(\beta+2)},
```

其中

```math
C_{\min}
=
\sin^{\beta+2}\left(\frac{\theta_{\min}}{2}\right).
```

当前程序中 $C_{\min}=0$。

## 5. Ali Rahmati 的能量损失公式

Rahmati 博士论文公式 2.19 给出入射热 O 的相对能量损失：

```math
\boxed{
\frac{\Delta E}{E}
=
\frac{2m_{\mathrm O}m_s}
{(m_{\mathrm O}+m_s)^2}
\left(1-\cos\theta_{\mathrm{COM}}\right)
}
```

这里 $m_{\mathrm O}$ 是入射 O 的质量，$m_s$ 是靶中性粒子的质量。逆变换抽到
的角度本身就是 COM 角，因此可以直接代入该公式，不需要进行 LAB 到 COM 的
角度转换。

对于 O 撞 O，

```math
\frac{\Delta E}{E}
=
\frac12(1-\cos\theta_{\mathrm{COM}})
=
\sin^2\left(\frac{\theta_{\mathrm{COM}}}{2}\right).
```

所以 COM 角为 $90^\circ$ 时损失 50%，为 $180^\circ$ 时损失 100%。

## 6. 从 COM 转换回相对火星静止的 LAB

碰撞前的质心速度为

```math
\mathbf V_{\mathrm{COM}}
=
\frac{m_{\mathrm O}\mathbf v_{\mathrm O}
+m_s\mathbf v_s}
{m_{\mathrm O}+m_s}.
```

当前模型采用静止靶近似：

```math
\mathbf v_s=0.
```

碰撞前相对速度为

```math
\mathbf g
=
\mathbf v_{\mathrm O}-\mathbf v_s.
```

程序按照抽到的 $\theta_{\mathrm{COM}}$ 和均匀方位角 $\phi$ 旋转
$\mathbf g$，得到 $\mathbf g'$。弹性碰撞中

```math
|\mathbf g'|=|\mathbf g|.
```

碰撞后的 LAB 速度为

```math
\mathbf v_{\mathrm O}'
=
\mathbf V_{\mathrm{COM}}
+
\frac{m_s}{m_{\mathrm O}+m_s}\mathbf g',
```

```math
\mathbf v_s'
=
\mathbf V_{\mathrm{COM}}
-
\frac{m_{\mathrm O}}{m_{\mathrm O}+m_s}\mathbf g'.
```

这些关系同时满足线动量守恒和总动能守恒。

## 7. 次级热 O

如果靶粒子是 O，碰撞后的 $\mathbf v_s'$ 属于反冲 O。如果其能量高于最低
追踪能量，程序会把它作为次级热 O 加入粒子队列，并继承相同宏粒子权重。

## 8. 一次碰撞的计算顺序

1. 用中性密度和总截面判断是否发生碰撞。
2. 按 $n_s\sigma_s$ 选择靶成分。
3. 生成均匀随机数 $R$。
4. 用逆 CDF 得到 $\theta_{\mathrm{COM}}$。
5. 独立抽取均匀方位角 $\phi$。
6. 在 COM 中旋转相对速度。
7. 转换回 LAB，得到入射粒子和靶粒子的碰撞后速度。
8. 在 LAB 中继续追踪热 O。
9. 如果靶是 O，则按需要追踪次级热 O。

## 9. 当前近似

1. 不同靶成分共同采用 Rahmati 对 Kharchenko O 与 O 结果的角分布拟合。
2. 靶中性粒子在碰撞前设为静止。
3. 不使用 10°截断，因此会显式模拟数量很多、单次能量损失很小的前向碰撞。
4. 后续获得各碰撞对的能量相关微分截面后，应分别替换共同角分布。

## 参考文献

Rahmati, A. (2016), Oxygen Exosphere of Mars: Evidence from Pickup Ions
Measured by MAVEN, PhD dissertation, University of Kansas, Sections 2.2.3
and 2.2.4.

Kharchenko, V., Dalgarno, A., Zygelman, B., and Yee, J. H. (2000), Energy
transfer in collisions of oxygen atoms in the terrestrial atmosphere,
Journal of Geophysical Research.
