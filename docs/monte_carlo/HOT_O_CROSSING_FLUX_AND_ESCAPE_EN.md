# Hot O Monte Carlo Transport, Directional Flux, and Escape Rate

## 1. Scope

This document describes the complete MarsHotO calculation from MGITM inputs to a hot O escape rate estimate. It covers:

1. Photochemical production and macroparticle weights
2. Sampling of nascent position, energy, and direction
3. Numerical propagation under Martian gravity
4. Neutral collisions, scattering angles, and two body kinematics
5. Production and tracking of secondary hot O
6. Recording of altitude surface crossing events
7. Calculation of upward and downward fluxes
8. Calculation of the escape rate at 300 km

Julia generates the source particles and performs particle transport. Python reads the binary crossing events, calculates directional fluxes, estimates statistical uncertainty, and produces the figures.

## 2. Configuration of the reported simulation

The reported calculation uses:

| Parameter | Value |
|---|---:|
| MGITM file | `MGITM/MGITM_LS000_F070_150901.dat` |
| Atmospheric profile | Nearest subsolar grid columns |
| Source altitude range | 100 to 250 km |
| Source altitude spacing | 1 km |
| Primary particles per altitude in each batch | 500 |
| Number of independent batches | 20 |
| Total primary particles | 1,510,000 |
| Total secondary O particles | 1,496,727 |
| Total tracked particles | 3,006,727 |
| Total crossing and terminal events | 98,810,553 |
| Computational domain | 100 to 2000 km |
| Crossing surface spacing | 10 km |
| Energy range | 0.01 to 7.0 eV |
| Number of energy bins | 140 |
| Random seeds | 20260810 through 20260829 |

The one dimensional nearest subsolar profile is extended as a spherically symmetric atmosphere. This calculation is therefore suitable for examining the physical transport and estimating fluxes, but its absolute global production rate depends on the spherical symmetry approximation.

## 3. Hot O production from MGITM

### 3.1 Dissociative recombination coefficient

The O2+ dissociative recombination coefficient is

```math
k(T_e)
=
1.95\times10^{-7}
\left(\frac{300}{T_e}\right)^{0.70}
\quad \mathrm{cm^3\,s^{-1}},
\qquad
T_e\le1200\ \mathrm{K},
```

```math
k(T_e)
=
7.39\times10^{-8}
\left(\frac{1200}{T_e}\right)^{0.56}
\quad \mathrm{cm^3\,s^{-1}},
\qquad
T_e>1200\ \mathrm{K}.
```

Each reaction produces two O atoms, so the volumetric hot O production rate is

```math
Q_{\mathrm{hot\,O}}(z)
=
2n_e(z)n_{\mathrm{O_2^+}}(z)
k[T_e(z)].
```

The program converts densities, lengths, cross sections, and reaction coefficients to SI units internally. Therefore,

```math
[Q_{\mathrm{hot\,O}}]
=
\mathrm{m^{-3}\,s^{-1}}.
```

### 3.2 Macroparticle weight

For source altitude $z_i$, the associated spherical shell volume is

```math
V_i
=
\frac{4\pi}{3}
\left[
(R_M+z_{i,+})^3
-
(R_M+z_{i,-})^3
\right].
```

The physical hot O production rate in that shell is

```math
\dot N_i
=
Q_{\mathrm{hot\,O}}(z_i)V_i.
```

If $N_i$ Monte Carlo primary macroparticles are generated at that altitude, each particle has rate weight

```math
w_i
=
\frac{\dot N_i}{N_i},
\qquad
[w_i]=\mathrm{s^{-1}}.
```

One simulated particle is not one physical O atom. It represents $w_i$ physical hot O atoms produced per second. A recoil secondary O inherits the same rate weight as its parent.

The 20 batches are independent Monte Carlo estimates of the same physical source. Their fluxes must be averaged. They must not be summed as if the batches represented different physical sources.

## 4. Nascent hot O velocity and energy

### 4.1 Maxwellian reactant distributions

At each source altitude, the electron and O2+ velocities are independently sampled from three dimensional Maxwellian distributions:

```math
f_s(\mathbf v)
=
\left(
\frac{m_s}{2\pi k_BT_s}
\right)^{3/2}
\exp
\left[
-\frac{
m_s|\mathbf v-\mathbf u_s|^2
}{
2k_BT_s
}
\right].
```

The electron and ion bulk velocities are zero:

```math
\mathbf u_e
=
\mathbf u_{\mathrm{O_2^+}}
=(0,0,0).
```

The three velocity components are independent and satisfy

```math
v_j
\sim
\mathcal N
\left(
0,
\frac{k_BT_s}{m_s}
\right).
```

### 4.2 Reactant COM velocity and relative energy

The reactant center of mass velocity is

```math
\mathbf V_{\mathrm{COM}}
=
\frac{
m_e\mathbf v_e+m_i\mathbf v_i
}{
m_e+m_i
}.
```

The reduced mass and relative translational energy are

```math
\mu
=
\frac{m_em_i}{m_e+m_i},
```

```math
E_{\mathrm{rel}}
=
\frac{1}{2}\mu
\left|
\mathbf v_e-\mathbf v_i
\right|^2.
```

### 4.3 Reaction branches and vibration

The four branches are:

| Products | Total released energy | Probability |
|---|---:|---:|
| O(3P) + O(3P) | 6.99 eV | 0.265 |
| O(1D) + O(3P) | 5.02 eV | 0.473 |
| O(1D) + O(1D) | 3.06 eV | 0.204 |
| O(1D) + O(1S) | 0.83 eV | 0.058 |

The program also samples the O2+ vibrational state $v$ from the configured population. Each vibrational quantum contributes 0.23 eV. The total energy available to the two products is

```math
E_{\mathrm{avail}}
=
E_{\mathrm{branch}}
+
E_{\mathrm{rel}}
+
0.23v.
```

For two equal mass O products, each receives half of this energy in the product center of mass frame. The product speed is

```math
u_O
=
\sqrt{
\frac{E_{\mathrm{avail}}}{m_O}
},
```

where the energy is converted to J in the program. The product direction is sampled isotropically. The velocity in the Mars stationary frame is

```math
\mathbf v_{O,\mathrm{LAB}}
=
\mathbf V_{\mathrm{COM}}
+
u_O\hat{\mathbf n}.
```

Here LAB means the Mars stationary reference frame.

## 5. Stepwise hot O transport

### 5.1 Mean free path and step length

At the current particle altitude and energy, the total collision coefficient is

```math
\kappa(E,z)
=
\sum_s n_s(z)\sigma_s(E).
```

The mean free path is

```math
\lambda(E,z)
=
\frac{1}{\kappa(E,z)}.
```

The current model contains O, CO, N2, O2, and CO2. The total cross sections are

```math
\sigma_s(E)
=
\sigma_s(3\ \mathrm{eV})
\left(
\frac{E}{3\ \mathrm{eV}}
\right)^{-0.2}.
```

The Rahmati step rule is

```math
ds
=
\begin{cases}
0.1\lambda,
& \lambda<10\ \mathrm{km},\\
1\ \mathrm{km},
& \lambda\ge10\ \mathrm{km}.
\end{cases}
```

### 5.2 Martian gravity

The gravitational acceleration is

```math
\mathbf a(\mathbf r)
=
-\frac{GM_M}{|\mathbf r|^3}\mathbf r.
```

For each step, the flight time is estimated from $dt=ds/|\mathbf v|$. Position and velocity are then updated with a velocity Verlet form. The particle follows a three dimensional trajectory and exchanges kinetic and gravitational potential energy as its altitude changes.

### 5.3 Collision decision

The collision probability during one numerical step follows the Rahmati procedure:

```math
P_{\mathrm{collision}}
=
\min(ds\,\kappa,1).
```

A uniform random number $R\in[0,1)$ is generated. A collision occurs when

```math
R<P_{\mathrm{collision}}.
```

When a collision occurs, the conditional probability of target species $s$ is

```math
P_s
=
\frac{n_s\sigma_s}
{\sum_jn_j\sigma_j}.
```

### 5.4 COM scattering angle

The current model uses the analytical fit of Rahmati to the Kharchenko O and O differential cross section:

```math
\frac{d\sigma}{d\Omega}
=
\alpha\sin^\beta
\left(
\frac{\theta_{\mathrm{COM}}}{2}
\right),
\qquad
\beta=-1.85.
```

The angular probability includes the solid angle Jacobian. The normalized polar angle probability is proportional to

```math
\frac{d\sigma}{d\Omega}
\sin\theta_{\mathrm{COM}}.
```

The program samples $\theta_{\mathrm{COM}}$ by inverse transform sampling and independently samples an azimuth uniformly over $[0,2\pi)$.

The present version applies this COM angular distribution to all neutral target species. Species specific differential cross sections can replace it when they become available.

### 5.5 Elastic two body collision

Before collision, the neutral target is stationary in the Mars frame. The center of mass velocity is

```math
\mathbf V_{\mathrm{COM}}
=
\frac{
m_O\mathbf v_O+m_s\mathbf v_s
}{
m_O+m_s
},
\qquad
\mathbf v_s=0.
```

The program rotates the relative velocity in COM while preserving its magnitude, and then transforms back to the Mars stationary frame:

```math
\mathbf v'_O
=
\mathbf V_{\mathrm{COM}}
+
\frac{m_s}{m_O+m_s}\mathbf g',
```

```math
\mathbf v'_s
=
\mathbf V_{\mathrm{COM}}
-
\frac{m_O}{m_O+m_s}\mathbf g'.
```

Here $\mathbf g'$ is the scattered relative velocity. These equations conserve total momentum and total kinetic energy.

The fractional energy loss of the incident hot O is

```math
\frac{\Delta E}{E}
=
\frac{2m_Om_s}{(m_O+m_s)^2}
\left(
1-\cos\theta_{\mathrm{COM}}
\right).
```

If the target is O and its recoil energy exceeds 0.01 eV, the recoil O is added to the tracking queue as a secondary hot O.

## 6. Termination conditions and raw events

A particle is terminated when any of the following occurs:

1. Kinetic energy is less than or equal to 0.01 eV.
2. Altitude falls below 100 km.
3. Altitude exceeds 2000 km.
4. The particle exceeds the maximum number of steps.
5. The particle queue exceeds its safety limit.

The directional flux calculation does not retain every integration step. It records only physically useful diagnostic events:

1. Particle birth
2. Crossing of a specified altitude surface
3. Exit through the lower boundary
4. Exit through the upper boundary
5. Thermalization
6. Reaching a numerical safety limit

Each event contains particle ID, parent ID, weight, time, altitude, three dimensional velocity, radial velocity, collision count, event type, and radial direction.

## 7. Directional flux from crossing events

At a spherical surface of radius

```math
r=R_M+z,
```

the crossing events are separated by energy bin and radial velocity:

```math
v_r>0
\quad\text{for upward crossings},
```

```math
v_r<0
\quad\text{for downward crossings}.
```

The directional flux in energy bin $k$ is

```math
\Phi_k(r)
=
\frac{
\sum_{p\in k}w_p
}{
4\pi r^2
}.
```

Its unit is

```math
\mathrm{cm^{-2}\,s^{-1}\ per\ energy\ bin}.
```

The result is not divided by the energy bin width. It is the flux contained in each energy bin, not a spectral density per eV.

The flux is calculated separately for each of the 20 independent batches. The final curve is the batch mean. Shading in the line plot gives the standard error of the mean across the 20 batches. Shading is omitted for bins with very large relative uncertainty.

### 7.1 Altitude and energy snapshots at fixed flight times

In addition to steady state altitude surface crossing events, MarsHotO provides fixed flight time snapshots to examine how an initially released hot O ensemble propagates. All primary hot O particles are released simultaneously at $t=0$. Particle states are recorded at

```math
t=0,\ 10,\ 50,\ 100\ \mathrm{s}.
```

Secondary O produced by collisions is included. The birth time of a secondary O is the actual collision time and is not reset to zero.

The reported snapshot calculation generates 1000 primary particles at every 1 km source altitude from 100 to 250 km, giving

```math
151\times1000
=
151{,}000
```

primary hot O particles. The figure uses 5 km altitude bins and 0.05 eV energy bins. Its displayed ranges are 100 to 1000 km and 0 to 7 eV.

For time $t$, altitude bin $i$, and energy bin $k$, the macroparticle rate weights occupying that bin are first summed:

```math
\dot N_{ik}(t)
=
\sum_{p\in(i,k,t)}w_p,
\qquad
[\dot N_{ik}]=\mathrm{s^{-1}}.
```

This rate is divided by the spherical area at the altitude bin center:

```math
\Phi_{ik}^{\mathrm{snap}}(t)
=
\frac{
\dot N_{ik}(t)
}{
4\pi(R_M+z_i)^2
}.
```

The unit is

```math
\mathrm{cm^{-2}\,s^{-1}\ per\ energy\ bin}.
```

The result is not divided by the energy bin width. The color scale shows

```math
\log_{10}\Phi_{ik}^{\mathrm{snap}}.
```

![Hot O altitude and energy snapshot flux at fixed flight times](../../examples/figures/hot_o_energy_altitude_time_snapshots.png)

Panels a through d show 0, 10, 50, and 100 s, respectively. The initial particles occupy 100 to 250 km. With increasing flight time, energetic particles propagate to higher altitudes and form a clear relation between altitude and energy. Altitude bins containing fewer than 20 particles are assigned the lowest color to suppress low sample noise.

The quantity $\Phi_{ik}^{\mathrm{snap}}$ is the macroparticle production rate occupying an altitude bin at a fixed time, divided by spherical area. It is therefore described as a snapshot flux estimate. It is not separated into upward and downward components and is not the net radial flux through a spherical surface. Escape rates must continue to use the altitude surface crossing flux defined in Section 7.

Run the snapshot calculation from the project root:

```bash
julia --project=. examples/run_hot_o_time_snapshots.jl \
  1000 20260730 examples/output/hot_o_time_snapshots.dat
```

Then create the figure:

```bash
C:\Users\Win\.conda\envs\mars\python.exe \
  examples/plot_hot_o_time_snapshots.py
```

The large local snapshot table remains in `examples/output/` and is not committed to GitHub. GitHub contains the simulation code, plotting code, and final PNG.

## 8. Directional flux from 100 to 300 km

![Directional hot O flux from 100 to 300 km](../../examples/figures/hot_o_directional_flux_100_300km.png)

The left panel shows upward flux and the right panel shows downward flux. Collisions are frequent at low altitudes, so both upward and downward populations are substantial. With increasing altitude, downward high energy particles become rare, while the upward population retains a pronounced high energy tail.

## 9. Energy spectra at 300 km

![Directional hot O energy spectra at 300 km](../../examples/figures/hot_o_directional_flux_spectrum_300km.png)

Summing over all energy bins gives an upward flux at 300 km of

```math
\Phi_{\mathrm{up}}
=
(1.49449\pm0.00506)\times10^8
\ \mathrm{cm^{-2}\,s^{-1}}.
```

The downward flux is

```math
\Phi_{\mathrm{down}}
=
(7.91631\pm0.03441)\times10^7
\ \mathrm{cm^{-2}\,s^{-1}}.
```

## 10. Local escape energy at 300 km

The local escape energy of one O atom at radius $r$ is

```math
E_{\mathrm{esc}}(r)
=
\frac{GM_Mm_O}{r}.
```

Using

```math
R_M=3389.5\ \mathrm{km},
\qquad
z=300\ \mathrm{km},
```

gives

```math
r
=
R_M+z
=
3689.5\ \mathrm{km},
```

```math
E_{\mathrm{esc}}(300\ \mathrm{km})
=
1.92484\ \mathrm{eV}.
```

Summing only upward particles whose energy bin centers satisfy

```math
E_k\ge E_{\mathrm{esc}}
```

gives the upward flux that is energetically capable of escape:

```math
\Phi_{\mathrm{esc}}
=
(2.88003\pm0.01458)\times10^7
\ \mathrm{cm^{-2}\,s^{-1}}.
```

## 11. Escape rate using projected area

The present projected area definition is

```math
A_{\mathrm{proj}}
=
\pi
\left(
R_M+300\ \mathrm{km}
\right)^2.
```

Converting the radius to centimeters gives

```math
r
=
3.6895\times10^8\ \mathrm{cm}.
```

Therefore,

```math
A_{\mathrm{proj}}
=
\pi r^2
=
4.27646\times10^{17}\ \mathrm{cm^2}.
```

If all upward particles are included,

```math
\dot N_{\mathrm{up,proj}}
=
\Phi_{\mathrm{up}}A_{\mathrm{proj}}
=
6.39114\times10^{25}\ \mathrm{s^{-1}}.
```

However, upward particles below the local escape energy remain gravitationally bound and should not be counted directly as escape. The projected area escape rate estimate is therefore

```math
\dot N_{\mathrm{esc,proj}}
=
\Phi_{\mathrm{esc}}A_{\mathrm{proj}},
```

```math
\boxed{
\dot N_{\mathrm{esc,proj}}
=
(1.23163\pm0.00624)\times10^{25}
\ \mathrm{s^{-1}}
}.
```

## 12. Spherically symmetric area comparison

The directional flux itself is defined using the spherical area $4\pi r^2$. If the nearest subsolar profile is interpreted as a globally spherically symmetric atmosphere, the global rate consistent with that geometric assumption is

```math
\dot N_{\mathrm{esc,global}}
=
\Phi_{\mathrm{esc}}4\pi r^2,
```

```math
\dot N_{\mathrm{esc,global}}
=
(4.92653\pm0.02494)\times10^{25}
\ \mathrm{s^{-1}}.
```

This is a spherical extrapolation. It is not a global escape rate obtained by integration through a complete three dimensional MGITM atmosphere.

The criterion $E\ge E_{\mathrm{esc}}$ indicates only that the particle has sufficient local mechanical energy at 300 km. Collisions above 300 km can still alter its final state. The result should therefore be described as an energy criterion escape rate estimate at 300 km.

## 13. Reproduction

### 13.1 Run 20 independent batches

From the project root, execute:

```bash
julia --project=. examples/run_hot_o_crossing_ensemble.jl \
  20 500 20260810 examples/output/run_1p51m_crossings
```

The arguments are:

1. Number of batches
2. Number of primary particles per source altitude in each batch
3. First random seed
4. Local output directory

There are 151 source altitudes. The total number of primary particles is

```math
20\times500\times151
=
1{,}510{,}000.
```

### 13.2 Calculate fluxes and create figures

```bash
C:\Users\Win\.conda\envs\mars\python.exe \
  examples/plot_directional_hot_o_flux.py \
  examples/output/run_1p51m_crossings
```

The raw binary event files are approximately 8.7 GB and remain local. GitHub contains the complete calculation code, fixed inputs, reproduction commands, final PNG files, and a small numerical summary.

The machine readable 300 km result is

```text
examples/results/hot_o_escape_flux_300km.json
```

## 14. Code map

| File | Purpose |
|---|---|
| `src/atmosphere.jl` | Read and interpolate the MGITM atmosphere |
| `src/chemistry.jl` | Dissociative recombination coefficient and hot O production |
| `src/source_particles.jl` | Maxwellian velocities, branches, vibration, and nascent hot O |
| `src/cross_sections.jl` | Total cross sections, collision coefficient, and target selection |
| `src/scattering.jl` | Rahmati COM angular PDF, CDF, and inverse sampling |
| `src/collision_kinematics.jl` | COM elastic collision and Mars frame velocities |
| `src/ensembles.jl` | Gravity propagation, step rule, and residence time estimator |
| `src/crossing_events.jl` | Macroparticle queue, secondary O, and crossing events |
| `examples/run_hot_o_crossing_ensemble.jl` | Entry point for the 20 batch Monte Carlo calculation |
| `examples/plot_directional_hot_o_flux.py` | Event processing, flux, uncertainty, escape rate, and figures |
| `examples/run_hot_o_time_snapshots.jl` | Generate fixed flight time snapshots at 0, 10, 50, and 100 s |
| `examples/plot_hot_o_time_snapshots.py` | Calculate area normalized snapshot flux and create the two by two figure |
| `test/runtests.jl` | Tests of Maxwellians, cross sections, scattering, conservation, and reproducibility |
