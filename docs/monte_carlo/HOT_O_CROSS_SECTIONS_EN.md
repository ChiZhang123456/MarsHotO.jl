# Hot O Collision Cross Sections With the Neutral Atmosphere

## 1. Role of the total collision cross section

For a hot O atom with energy $E$, the local collision coefficient for target species $s$ is

```math
\kappa_s(E,z)=n_s(z)\sigma_s(E).
```

The total collision coefficient and mean free path are

```math
\kappa_{\mathrm{tot}}(E,z)
=\sum_s n_s(z)\sigma_s(E),
```

```math
\lambda(E,z)
=\frac{1}{\kappa_{\mathrm{tot}}(E,z)}.
```

The model uses $\mathrm{m^{-3}}$ for $n_s$, $\mathrm{m^2}$ for $\sigma_s$, $\mathrm{m^{-1}}$ for $\kappa_{\mathrm{tot}}$, and m for $\lambda$.

## 2. Total cross sections currently used

The reference energy is 3 eV. The configuration file `data/cross_sections/rahmati_total_cross_sections.toml` contains:

| Target species | $\sigma(3\ \mathrm{eV})$, $\mathrm{cm^2}$ |
|---|---:|
| O | $6.4\times10^{-15}$ |
| CO | $1.8\times10^{-14}$ |
| N2 | $1.8\times10^{-14}$ |
| O2 | $1.8\times10^{-14}$ |
| CO2 | $2.0\times10^{-14}$ |

The energy dependence is

```math
\sigma_s(E)
=
\sigma_s(3\ \mathrm{eV})
\left(\frac{E}{3\ \mathrm{eV}}\right)^{-0.2}.
```

Neither the MGITM input used here nor the current collision configuration contains Ar. Therefore, the transport model does not include hot O collisions with Ar.

## 3. Collision decision during one step

The adaptive step rule used in the Rahmati transport procedure is:

1. If $\lambda<10\ \mathrm{km}$, use $ds=0.1\lambda$.
2. If $\lambda\ge10\ \mathrm{km}$, use $ds=1\ \mathrm{km}$.

The linear collision probability during one step is

```math
P_{\mathrm{coll}}
=
\min\left(ds\,\kappa_{\mathrm{tot}},1\right).
```

A uniform random number $U$ is sampled. A collision occurs during the step when

```math
U<P_{\mathrm{coll}}.
```

## 4. Target species selection

After a collision has been selected, the conditional probability that the target is species $s$ is

```math
P(s\mid\mathrm{coll})
=
\frac{n_s(z)\sigma_s(E)}
{\sum_j n_j(z)\sigma_j(E)}.
```

Target selection therefore depends on both the local neutral density and the total cross section of the collision pair.

The current angular distribution uses the analytical fit of Rahmati to the Kharchenko O and O differential cross section:

```math
\frac{d\sigma}{d\Omega}
=
\alpha\sin^\beta
\left(\frac{\theta_{\mathrm{COM}}}{2}\right),
\qquad \beta=-1.85.
```

The angular probability density includes the solid angle Jacobian $\sin\theta_{\mathrm{COM}}$. The current model samples the complete interval

```math
0\le\theta_{\mathrm{COM}}\le\pi
```

without a 10 degree cutoff. See [COM scattering angle and two body collisions](HOT_O_SCATTERING_TWO_BODY_EN.md) for the angular sampling and collision kinematics.
