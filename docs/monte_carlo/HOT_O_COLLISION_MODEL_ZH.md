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

当前模型采用 Rahmati 对 Kharchenko O 与 O 微分截面的解析拟合，
在 COM 中抽取散射角。默认最小角为零，因此包含 10 度以内的小角散射。
碰撞后速度由静止靶两体弹性碰撞的动量和能量守恒关系计算，再从 COM
转换回静止坐标系。

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
    -> 用解析逆 CDF 抽取 Rahmati COM 散射角，并抽取均匀方位角
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
| Rahmati COM 散射角解析模型 | `src/scattering.jl` |
| MGITM 大气 | `MGITM/` |
| 初始粒子 | `src/source_particles.jl` |
| 散射角抽样 | `src/scattering.jl` |
| 两体碰撞 | `src/collision_kinematics.jl` |
| 单粒子传输 | `src/transport.jl` |
| Monte Carlo 系综 | `src/ensembles.jl` |
| 完整示例 | `examples/run_hot_o_corona.jl` |

Rahmati 论文用于传输步长、碰撞概率、COM 散射角、能量损失关系和终止条件等整体
流程。`data/cross_sections/scattering_angle_distribution.txt` 仅保留为 MarsASPEN
参考数据，不参与当前运行时抽样。
