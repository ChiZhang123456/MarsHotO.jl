# 热 O 与中性大气的碰撞截面

## 1. 总碰撞截面的作用

对于能量为 $E$ 的热 O，目标成分 $s$ 的局地碰撞系数为

```math
\kappa_s(E,z)=n_s(z)\sigma_s(E),
```

总碰撞系数和平均自由程分别为

```math
\kappa_{\mathrm{tot}}(E,z)=\sum_s n_s(z)\sigma_s(E),
```

```math
\lambda(E,z)=\frac{1}{\kappa_{\mathrm{tot}}(E,z)}.
```

其中 $n_s$ 使用 $\mathrm{m^{-3}}$，$\sigma_s$ 使用 $\mathrm{m^2}$，
$\kappa_{\mathrm{tot}}$ 使用 $\mathrm{m^{-1}}$，$\lambda$ 使用 m。

## 2. 当前使用的总截面

参考能量为 3 eV。配置文件
`data/cross_sections/rahmati_total_cross_sections.toml` 包含：

| 目标成分 | $\sigma(3\ \mathrm{eV})$，$\mathrm{cm^2}$ |
|---|---:|
| O | $6.4\times10^{-15}$ |
| CO | $1.8\times10^{-14}$ |
| N2 | $1.8\times10^{-14}$ |
| O2 | $1.8\times10^{-14}$ |
| CO2 | $2.0\times10^{-14}$ |

能量依赖采用

```math
\sigma_s(E)=\sigma_s(3\ \mathrm{eV})
\left(\frac{E}{3\ \mathrm{eV}}\right)^{-0.2}.
```

MGITM 输入和当前碰撞配置都不包含 Ar，因此传输模型不追踪热 O 与 Ar 的
碰撞。

## 3. 每一步如何决定是否碰撞

Rahmati 传输流程的自适应步长为：

1. 当 $\lambda<10\ \mathrm{km}$ 时，$ds=0.1\lambda$。
2. 当 $\lambda\ge10\ \mathrm{km}$ 时，$ds=1\ \mathrm{km}$。

一步内的线性碰撞概率为

```math
P_{\mathrm{coll}}
=
\min\left(ds\,\kappa_{\mathrm{tot}},1\right).
```

抽取均匀随机数 $U$。当 $U<P_{\mathrm{coll}}$ 时，本步发生碰撞。

## 4. 如何选择靶成分

已经确定发生碰撞以后，与成分 $s$ 碰撞的条件概率为

```math
P(s\mid\mathrm{coll})
=
\frac{n_s(z)\sigma_s(E)}
{\sum_j n_j(z)\sigma_j(E)}.
```

所以目标选择同时取决于局地中性密度和该碰撞对的总截面。

## 5. 不再使用角度截断

当前角度由 Kallio 与 Barabash (2001) 的 MarsASPEN 逆 CDF 表抽取。程序使用
表中的完整角度范围，不设置 10° 截断。因此：

1. $\kappa_{\mathrm{tot}}$ 直接使用完整总截面。
2. 不再计算角度保留比例。
3. 不再用角度保留比例缩小平均碰撞频率。
4. 每次被判定的碰撞都从完整查找表抽取一个 LAB 散射角。

角度抽样和碰撞后能量计算见
[LAB 散射角与两体碰撞](HOT_O_SCATTERING_TWO_BODY_ZH.md)。

## 6. 当前近似

总截面的 $E^{-0.2}$ 能量依赖是统一参数化。角度查找表原本描述高能 H 或
H ENA 在火星大气中的散射，目前按用户指定共同用于热 O 与所有中性成分。
因此当前模型适合建立可重复的基线模拟。更高精度的计算需要为每个 O 与中性
成分碰撞对提供各自的能量相关总截面和微分截面。
