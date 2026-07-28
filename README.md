# MarsHotO.jl

Mars Hot Oxygen Transport model, abbreviated MHOT, is a Julia particle-transport package for the production, collisional transport, corona formation, and photochemical escape of hot atomic oxygen at Mars. It supports both direct single-particle propagation and weighted Monte Carlo ensembles.

The current repository contains the Python prototype used to:

1. Read vertical profiles from Mars Global Ionosphere Thermosphere Model, MGITM, output.
2. Calculate hot O production from O2+ dissociative recombination.
3. Sample the four non-negligible reaction branches.
4. Calculate nascent hot O energy distributions using electron and ion thermal velocities.
5. Include the Martian O2+ vibrational distribution from Fox and Hać (1997).
6. Plot MGITM profiles, discrete reaction channels, nascent energy probabilities, and production rates in fixed energy bins.

The Julia core reads MGITM atmospheres, samples O2+ dissociative recombination source particles, samples the Rahmati and Kharchenko analytical COM scattering-angle distribution, performs conservative two-body collision kinematics, and transports particles using collision optical depth and Martian gravity. Python remains responsible for input inspection, validation, analysis, and scientific plotting.

## Package layout

```text
src/                  Julia atmosphere, source, collision, and transport code
data/chemistry/       O2+ dissociative recombination settings
data/cross_sections/  Hot O collision settings
data/atmosphere/      Atmosphere data notes
MGITM/                Packaged MGITM model output
examples/             Python scientific plotting examples
test/                 Julia numerical and physics tests
```

## Collision physics

MarsHotO uses the Rahmati fit to the Kharchenko O and O differential cross section,

```math
\frac{d\sigma}{d\Omega}=\alpha\sin^\beta(\theta_{\mathrm{COM}}/2),
\qquad \beta=-1.85.
```

The probability density includes the solid-angle Jacobian
`2 pi sin(theta_COM)`. The complete interval from 0 to 180 degrees is sampled, with no 10 degree cutoff. The COM angle is generated analytically by inverse transform sampling, the relative velocity is rotated in COM, and the projectile and target velocities are then transformed back to the stationary Mars LAB frame. Post-collision momentum and total kinetic energy are conserved.

物理模型的中文逐步说明见：

* [MarsHotO Monte Carlo 物理模型总览](docs/monte_carlo/HOT_O_COLLISION_MODEL_ZH.md)
* [热 O 与中性大气的碰撞截面](docs/monte_carlo/HOT_O_CROSS_SECTIONS_ZH.md)
* [热 O 高度和初生能量分布](docs/monte_carlo/HOT_O_SOURCE_MODEL_ZH.md)
* [LAB、COM、散射角和碰撞能量损失](docs/monte_carlo/HOT_O_SCATTERING_TWO_BODY_ZH.md)

## Current physics

The Julia package provides `transport_particle!` for direct single-particle propagation and `run_hot_o_corona` for weighted three-dimensional Monte Carlo ensembles.

The total hot O production rate is

```math
Q_{\mathrm{hot\,O}}
=2n_en_{\mathrm{O_2^+}}k(T_e).
```

The four reaction channels use the following total released energies and branching ratios:

| Products | Released energy | Branching ratio |
|---|---:|---:|
| O(3P) + O(3P) | 6.99 eV | 26.5% |
| O(1D) + O(3P) | 5.02 eV | 47.3% |
| O(1D) + O(1D) | 3.06 eV | 20.4% |
| O(1D) + O(1S) | 0.83 eV | 5.8% |

The prototype also supports the Mars exobase O2+ vibrational distribution reported by Fox and Hać (1997). Each vibrational quantum adds approximately 0.23 eV to the total reaction exothermicity.

See [HOT_OXYGEN_MODEL_PLAN.md](HOT_OXYGEN_MODEL_PLAN.md) and [AGENTS.md](AGENTS.md) for the detailed model plan and fixed physical conventions.

## Examples and scripts

### `examples/plot_collision_physics.py`

Plots the Rahmati and Kharchenko analytical inverse CDF, fractional energy loss versus COM scattering angle, and total collision cross sections versus energy.

### `examples/plot_mgitm_profiles.py`

Plots O2+ density, neutral/ion/electron temperatures, and hot O production rate versus altitude.

### `examples/shared/plot_thermal_energy_sampling.py`

Plots normalized 300 K Maxwellian velocity-component and kinetic-energy
probability densities together with the isotropic direction-cosine distribution.

### Complete Rahmati Monte Carlo example

Run the Julia transport model and then plot the residence-time altitude-energy
density:

```bash
julia --project=. examples/run_hot_o_corona.jl 10000
python examples/plot_hot_o_corona.py
```

The particle count is configurable. A production calculation can use
approximately `10000000` primary particles, but should be run on an HPC
system. The packaged MGITM profile ends near 250 km, so the current example
uses a decreasing log-linear density extrapolation above the model top.
The nearest-subsolar column is extended spherically, so its absolute global
source rate is a model approximation rather than a full three-dimensional
MGITM result.

### `plot_mgitm_hot_o_profiles.py`

Plots:

1. O2+ density versus altitude
2. Neutral, ion, and electron temperatures versus altitude
3. Total hot O production rate versus altitude

### `plot_discrete_hot_o_energy_lines.py`

Plots the four discrete hot O source energies without thermal broadening. The left panel shows fixed branching probabilities. The right panel shows the altitude-dependent production rate of each branch.

### `reproduce_hot_o_nascent_energy_map.py`

Calculates and plots:

1. The dimensionless conditional probability in each 0.025 eV energy bin,
   $P_k(z)=N_k(z)/N_{\mathrm{tot}}(z)$
2. The production rate in each 0.025 eV energy bin,
   $Q_k(z)=Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z)P_k(z)$, in
   m$^{-3}$ s$^{-1}$

Neither panel divides by the energy-bin width. Consequently,
$\sum_k P_k(z)=1$ and
$\sum_k Q_k(z)=Q_{\mathrm{hotO}}^{(\mathrm{m^{-3}})}(z)$.
The plotted total production rate is converted from cm$^{-3}$ s$^{-1}$ to
m$^{-3}$ s$^{-1}$ by multiplying by $10^6$.

The calculation currently includes:

1. Zero-bulk three-dimensional Maxwellian electron velocities at $T_e$
2. Zero-bulk three-dimensional Maxwellian O2+ velocities at $T_i$
3. Center-of-mass two-body kinematics
4. Isotropic product directions
5. The Fox and Hać (1997) Mars O2+ vibrational distribution

The default calculation uses $10^6$ dissociative recombination events at each MGITM altitude. Events are processed in batches of $10^5$ to limit memory use. Heatmaps use bilinear display interpolation, similar to MATLAB `shading interp`. The interpolation affects only figure rendering and does not modify the saved source-data grid.

It does not yet include the velocity-dependent dissociative recombination acceptance probability or O2+ rotational energy.

## Input data

The scripts currently use:

```text
MGITM/MGITM_LS000_F070_150901.dat
```

The `MGITM` directory contains all 12 model cases spanning four seasons and three solar activity levels. The current prototype uses the `Ls = 0 degrees`, F070 case as its baseline.

The baseline case uses the two nearest subsolar grid columns:

```text
Longitude: 27.5 degrees
Latitude: -2.5 and 2.5 degrees
Solar zenith angle: 2.73 degrees
```

The two profiles are averaged.

## Python environment

Install NumPy, pandas, and Matplotlib:

```bash
python -m pip install numpy pandas matplotlib
```

Run the scripts:

```bash
python plot_mgitm_hot_o_profiles.py
python examples/plot_mgitm_profiles.py
python examples/plot_collision_physics.py
python plot_discrete_hot_o_energy_lines.py
python reproduce_hot_o_nascent_energy_map.py
```

Generated figures and source-data tables are written to `figures`.

Run the Julia physics tests with:

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

## Current results

### MGITM profiles and total production

![MGITM profiles](examples/figures/mgitm_ls000_f070_profiles.png)

### Hot O collision physics

![Hot O collision physics](examples/figures/hot_o_collision_cross_sections_and_scattering.png)

### Nascent hot O probability and production rate per energy bin

![Hot O energy maps](examples/figures/mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_energy_maps.png)

### Monte Carlo hot O altitude-energy distribution

![Monte Carlo hot O altitude-energy distribution](examples/figures/hot_o_altitude_energy_distribution.png)

## References

1. Fox, J. L., and Hać, A. (1997), Spectrum of hot O at the exobases of the terrestrial planets, Journal of Geophysical Research, 102(A11), 24005-24011, https://doi.org/10.1029/97JA02089.
2. Lillis, R. J., et al. (2017), Photochemical escape of oxygen from Mars: First results from MAVEN in situ data, Journal of Geophysical Research: Space Physics, 122, 3815-3836, https://doi.org/10.1002/2016JA023525.
3. Rahmati, A. (2016), Oxygen Exosphere of Mars: Evidence from Pickup Ions Measured by MAVEN, PhD dissertation, University of Kansas.
4. Kallio, E., and Barabash, S. (2001), Atmospheric effects of precipitating energetic hydrogen atoms on the Martian atmosphere, Journal of Geophysical Research: Space Physics, 106(A1), 165 to 177, https://doi.org/10.1029/2000JA002003.
