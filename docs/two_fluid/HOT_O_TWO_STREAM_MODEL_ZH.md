# MarsHotO 双流输运模型

## 1. 模型定位

`MarsHotO.jl` 同时包含两套热 O 输运方法：

1. Rahmati Monte Carlo 粒子输运
2. Rahmati 第 2.2.1 节的一维双流输运

这里的 two-stream 指向上和向下两个热 O 通量流，不是等离子体双流体模型。双流方法求解能量和高度空间中的向上通量与向下通量，适合快速计算热 O 通量和光化学逃逸率。

当前模型使用 MGITM 最近日下点剖面，并将该一维剖面当作球对称背景。这是一项模型近似。

## 2. 控制方程

对每个能量区间求解

```math
\pm \langle\mu\rangle
\frac{d\phi^\pm(E,z)}{dz}
=
-\phi^\pm(E,z)\sum_s n_s(z)\sigma_s(E)
+\frac{P(E,z)}{2}
+P_{\mathrm{cascade}}^\pm(E,z)
+P_{\mathrm{secondary}}^\pm(E,z).
```

其中：

* `phi+` 是向上通量
* `phi-` 是向下通量
* `mu = 0.5` 表示半球内各向同性通量的平均俯仰角余弦
* `P/2` 是分配到每个方向的初级热 O 产生率
* `P_cascade` 是高能热 O 碰撞后进入较低能量区间的源
* `P_secondary` 是热 O 与背景 O 碰撞后产生的反冲热 O 源

通量和产生率均对单个能量 bin 积分。因此，内部单位分别为 `m^-2 s^-1` 和 `m^-3 s^-1`。

## 3. 大气和化学输入

默认输入为：

```text
MGITM/MGITM_LS000_F070_150901.dat
```

所需变量包括 CO2、O、N2、CO、O2、O2+ 和电子数密度，以及电子温度。

总热 O 产生率为

```math
Q_{\mathrm{hot\,O}}=2n_en_{\mathrm{O_2^+}}k(T_e).
```

四个解离复合分支和 O2+ 振动态分布从 `data/chemistry/o2plus_dissociative_recombination.toml` 读取。由于两个产物质量相同，每个 O 获得总释放能量和振动能量的一半。当前 MGITM 文件不包含 Ar，因此默认碰撞目标不包含 Ar。

## 4. 碰撞和能量重分布

总截面使用

```math
\sigma_s(E)=\sigma_s(3\,\mathrm{eV})(E/3\,\mathrm{eV})^{-0.2}.
```

`build_two_stream_redistribution` 使用现有 LAB 两体碰撞模块建立能量和方向重分布矩阵。每个入射能量和目标成分都使用固定随机种子采样：

1. 从 MarsASPEN 逆 CDF 表抽取数值，并解释为经验 COM 散射角
2. 在运动学允许范围内对散射角分布进行条件化
3. 抽取均匀方位角
4. 用 LAB 动量和总动能守恒计算碰撞后速度
5. 根据碰撞后垂直速度符号分入同向流或反向流
6. 将主粒子碰撞后能量分入对应能量 bin
7. 当目标为 O 时，将反冲 O 加入次级重分布矩阵

本实现使用完整散射角范围和完整总截面，不设置 10 度最小散射角，不缩放总截面。这与 `MarsHotO.jl` 的统一碰撞约定一致。它与 Rahmati 论文中使用 10 度截断的历史实现存在明确差异。

## 5. 数值积分和边界条件

在一个高度步长内，假定碰撞系数和源项为常数。解析更新为

```math
\phi_{i+1}
=
\phi_i\exp(-\kappa_i\Delta z/\langle\mu\rangle)
+\frac{S_i}{\kappa_i}
\left[1-\exp(-\kappa_i\Delta z/\langle\mu\rangle)\right].
```

默认高度范围为 100 至 300 km，步长为 1 km。边界条件为：

* 100 km：`phi+ = phi-`，表示频繁碰撞造成的各向同性下边界
* 300 km：能量低于局地逃逸能的上行粒子反射为下行粒子
* 300 km：能量高于局地逃逸能的上行粒子离开模型并计入逃逸通量
* 顶部以上的向下源使用 Rahmati 式 2.10 的指数尺度高度近似

碰撞后仍位于同一能量 bin 的粒子会形成 bin 内耦合。求解器对全部能量和两个方向进行迭代，直到相对通量变化满足容差。

完整小角散射使同一能量 bin 内的耦合很强，尤其是低于逃逸能的反射粒子。因此默认最大迭代次数为 2000，相对通量容差为 `1e-3`。提高精度时可减小 `relative_tolerance`，但应同时检查逃逸通量随迭代次数和能量分辨率是否稳定。

## 6. Julia 接口

```julia
using MarsHotO

profile = load_mgitm_subsolar_profile("MGITM/MGITM_LS000_F070_150901.dat")
targets = load_collision_targets(
    "data/cross_sections/rahmati_total_cross_sections.toml"
)
branches = load_reaction_branches(
    "data/chemistry/o2plus_dissociative_recombination.toml"
)

config = TwoStreamConfig(
    altitude_min_m=100e3,
    altitude_max_m=300e3,
    altitude_step_m=1e3,
    energy_edges_eV=collect(range(0.01, 7.01; length=71)),
)

result = run_two_stream(profile, targets, branches; config=config)
```

主要输出包括 `upward_flux_m2_s1`、`downward_flux_m2_s1`、`primary_production_m3_s1`、`escape_energy_eV`、`escape_flux_m2_s1` 和 `converged`。

运行示例：

```text
julia --project=. examples/two_fluid/run_two_stream.jl
```

输出写入 `examples/output/hot_o_two_stream_flux.csv`。

## 7. 验证重点

1. 所有通量和产生率非负
2. 固定随机子的重分布矩阵可重复
3. 每个主粒子重分布概率不超过 1
4. O 目标的反冲次级粒子能量来自守恒的碰撞后速度
5. 下边界满足向上和向下通量相等
6. 低于逃逸能的顶部通量满足反射边界
7. 逃逸通量只包含顶部高于局地逃逸能的向上通量
8. 增加能量和高度分辨率后结果趋于稳定

## 8. 参考文献

Rahmati, A. (2016), *Oxygen Exosphere of Mars: Evidence from Pickup Ions Measured by MAVEN*, PhD dissertation, University of Kansas, Section 2.2.1.
