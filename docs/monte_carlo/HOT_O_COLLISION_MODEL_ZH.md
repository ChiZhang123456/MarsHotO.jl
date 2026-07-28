# MarsHotO Monte Carlo 物理模型总览

MarsHotO 由热 O 源模型和碰撞传输模型两部分组成。

## 1. 热 O 源模型

源模型读取 MGITM 的 $n_e$、$n_{\mathrm{O_2^+}}$、$T_e$ 和 $T_i$，
计算 O2+ 解离复合产生率，抽取四个反应分支和 O2+ 振动态，并在 LAB 中生成
热 O 的初始位置、能量和各向同性速度方向。

详细说明见
[热 O 高度和初生能量分布](HOT_O_SOURCE_MODEL_ZH.md)。

## 2. 碰撞频率和靶成分

中性密度和能量相关总截面共同决定平均自由程、一步内的碰撞概率，以及发生
碰撞时的靶成分。

详细说明见
[热 O 与中性大气的碰撞截面](HOT_O_CROSS_SECTIONS_ZH.md)。

## 3. 散射角和两体运动学

当前模型从 MarsASPEN 的 Kallio 与 Barabash (2001) 查找表抽取角度数值，
并将其作为经验 COM 散射角。碰撞后速度由静止靶两体弹性碰撞的动量和能量
守恒关系计算。

详细说明见
[COM 散射角与两体碰撞](HOT_O_SCATTERING_TWO_BODY_ZH.md)。

## 4. 完整 Monte Carlo 流程

```text
读取 MGITM 大气和等离子体剖面
    -> 计算 O2+ 解离复合产生率 Q_hotO(z)
    -> 抽取源高度、反应分支、振动态和初生速度
    -> 由中性密度和总截面计算平均自由程
    -> 在火星重力下推进粒子
    -> 抽样本步是否发生碰撞
    -> 按 n_s sigma_s 选择靶成分
    -> 从查找表抽取经验 COM 散射角，并抽取均匀方位角
    -> 在 COM 中散射并转换到 LAB 碰撞后速度
    -> 继续追踪主粒子，必要时追踪次级 O
    -> 用驻留时间估计高度和能量分布
```

源粒子按

```math
4\pi r^2 Q(z)\,dz
```

抽样。每个初级宏粒子代表的产生率是总球对称源率除以初级粒子数。次级 O
继承相同宏粒子率。最近日下点剖面被球对称扩展，这是当前全局计算的一项模型
近似。

## 5. 主要输入与代码

| 内容 | 文件 |
|---|---|
| 解离复合与振动态 | `data/chemistry/o2plus_dissociative_recombination.toml` |
| 总碰撞截面 | `data/cross_sections/rahmati_total_cross_sections.toml` |
| 经验 COM 散射角逆 CDF | `data/cross_sections/scattering_angle_distribution.txt` |
| MGITM 大气 | `MGITM/` |
| 初始粒子 | `src/monte_carlo/source_particles.jl` |
| 散射角抽样 | `src/shared/scattering.jl` |
| 两体碰撞 | `src/shared/collision_kinematics.jl` |
| 单粒子传输 | `src/monte_carlo/transport.jl` |
| Monte Carlo 系综 | `src/monte_carlo/ensembles.jl` |
| 完整示例 | `examples/monte_carlo/run_hot_o_corona.jl` |

Rahmati 论文仍用于传输步长、碰撞概率、COM 能量损失关系和终止条件等整体
流程。当前角分布是把 Kallio 查找表的数值经验性地解释为 COM 角。
