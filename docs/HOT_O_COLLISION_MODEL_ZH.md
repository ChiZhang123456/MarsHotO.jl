# MarsHotO 物理模型中文说明

MarsHotO 将火星热 O 模型分成两个相互连接的部分。

1. 源模型计算热 O 在什么高度产生，以及刚产生时具有多少能量。
2. 传输模型追踪热 O 的运动、碰撞、能量损失和最终状态。

为了避免把不同物理过程混在一起，详细说明分成以下三个文档。

## 1. 碰撞截面和碰撞概率

[热 O 与中性大气的碰撞截面](HOT_O_CROSS_SECTIONS_ZH.md)

这一部分回答：

* MGITM 提供了哪些中性成分。
* 碰撞截面 $\sigma_s(E)$ 是什么。
* 如何由中性密度和碰撞截面得到平均自由程。
* Monte Carlo 模型如何决定何时碰撞，以及与哪一种粒子碰撞。

## 2. 热 O 的产生率和初始能量

[热 O 高度和初生能量分布](HOT_O_SOURCE_MODEL_ZH.md)

这一部分回答：

* 如何由 $n_e$、$n_{\mathrm{O_2^+}}$ 和 $T_e$ 得到热 O 产生率。
* 四个解离复合分支如何产生四个基本能量峰。
* $T_e$ 和 $T_i$ 如何通过 Monte Carlo 速度抽样展宽能量峰。
* $\mathrm{O_2^+}$ 振动态如何增加反应可用能量。
* 如何得到 $p(E\mid z)$ 和 $Q(E,z)$。

## 3. 散射角和两体碰撞

[LAB、COM、散射角和碰撞能量损失](HOT_O_SCATTERING_TWO_BODY_ZH.md)

这一部分从最基本的概念开始解释：

* LAB 和 COM 分别是什么。
* 文献中的散射角究竟是哪两个方向之间的夹角。
* 微分碰撞截面如何转换成散射角概率。
* 为什么小角散射损失的能量少，大角散射损失的能量多。
* 如何由动量守恒和能量守恒计算碰撞后的两个速度。
* 为什么 O 靶粒子可能成为次级热 O。

## 4. 完整计算流程

```text
读取 MGITM 大气和等离子体剖面
    ↓
计算 O₂⁺ 解离复合产生率 Q_hotO(z)
    ↓
抽取反应分支、振动态、电子速度和 O₂⁺ 速度
    ↓
计算 LAB 系中的热 O 初始位置和速度
    ↓
由中性密度和总碰撞截面计算平均自由程
    ↓
沿火星重力轨迹移动，并抽样是否发生碰撞
    ↓
抽取碰撞成分、COM 散射角和方位角
    ↓
由两体弹性碰撞关系计算碰撞后的 LAB 速度
    ↓
继续追踪主粒子，必要时追踪次级 O
    ↓
累计驻留时间，得到热 O 的 n(E,z)、n(z) 和速度分布
```

源模型的输出 $Q(E,z)$ 是传输模型的输入。碰撞截面决定碰撞频率，散射角分布和两体运动学决定每次碰撞后速度如何变化。

## 5. 代码和输入文件

| 物理内容 | 文件 |
|---|---|
| 解离复合分支和振动态 | `data/chemistry/o2plus_dissociative_recombination.toml` |
| 总碰撞截面 | `data/cross_sections/rahmati_total_cross_sections.toml` |
| 散射角分布 | `data/cross_sections/scattering_angle_distribution.txt` |
| MGITM 大气 | `data/atmosphere/` |
| 热 O 初始粒子 | `src/source_particles.jl` |
| 单粒子传输和碰撞 | `src/transport.jl` |
| Monte Carlo 系综 | `src/ensembles.jl` |
| 完整示例 | `examples/run_hot_o_corona.jl` |

当前模型采用 Rahmati 博士论文列出的 Monte Carlo 传输框架，并使用 Lillis 等人的热 O 源模型思路。具体采用了哪些近似，以及哪些输入还需要改进，分别写在上述三个文档中。
