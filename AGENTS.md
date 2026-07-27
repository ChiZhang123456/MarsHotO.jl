# Project Instructions

## Scientific Objective

This project studies the Martian hot oxygen corona produced by photochemical reactions.

The primary objective is to develop a hot O transport model that:

1. Uses MGITM thermosphere and ionosphere results as the background atmosphere.
2. Simulates hot O production from O2+ dissociative recombination.
3. Transports hot O atoms through collisions and the Martian gravitational field.
4. Calculates hot O density profiles, energy distributions, escape probabilities, escape fluxes, and global escape rates.

Use Julia for the Monte Carlo transport model and Python for preprocessing, validation, analysis, and plotting.

## Required Scientific References

Use the following local references as the primary methodological sources:

1. `Lillis-Photochemical escape of oxygen from Mar.pdf`
2. `Rahmati_ku_0099D_14448_DATA_1.pdf`
3. `HOT_OXYGEN_MODEL_PLAN.md`

The Lillis model provides the hot O production and escape probability framework. The Rahmati dissertation provides both the two-stream/Liouville method and the three-dimensional Monte Carlo transport method.

## Fixed O2+ Dissociative Recombination Channels

The hot oxygen source must use the following four non-negligible branches:

| Branch | Final state | Released energy | Branching ratio |
|---|---|---:|---:|
| 1 | O(3P) + O(3P) | 6.99 eV | 26.5% |
| 2 | O(1D) + O(3P) | 5.02 eV | 47.3% |
| 3 | O(1D) + O(1D) | 3.06 eV | 20.4% |
| 4 | O(1D) + O(1S) | 0.83 eV | 5.8% |

The branching ratios sum to 100% and must not be replaced without explicit user approval.

For Monte Carlo branch selection, the default cumulative probability intervals are:

```text
Branch 1: [0.000, 0.265)
Branch 2: [0.265, 0.738)
Branch 3: [0.738, 0.942)
Branch 4: [0.942, 1.000]
```

For the initial baseline model, divide the released energy equally between the two identical O atoms in the center-of-mass frame. The corresponding energy of each O atom is:

```text
Branch 1: 3.495 eV
Branch 2: 2.510 eV
Branch 3: 1.530 eV
Branch 4: 0.415 eV
```

The two product O atoms have opposite velocity directions in the center-of-mass frame. A later model may include the thermal velocities of O2+ and electrons using MGITM ion and electron temperatures.

## O2+ Vibrational Distribution

The hot O nascent energy model must support the Martian exobase O2+ vibrational distribution reported by Fox and Hać (1997):

| Vibrational quantum number \(v\) | Published fraction |
|---:|---:|
| 0 | 0.800 |
| 1 | 0.074 |
| 2 | 0.043 |
| 3 | 0.035 |
| 4 | 0.025 |
| 5 | 0.015 |
| 6 | 0.0047 |
| 7 | 0.00027 |
| 8 | 0.00021 |

The published fractions are rounded and do not sum exactly to one. Normalize them before Monte Carlo sampling:

\[
P_v=\frac{f_v}{\sum_j f_j}.
\]

Each vibrational quantum adds approximately 0.23 eV to the total reaction exothermicity:

\[
E_{\mathrm{available}}
=E_{\mathrm{branch}}+0.23v+E_{\mathrm{relative}},
\]

where all energies are in eV. Because the two product O atoms have equal mass, this additional vibrational energy is divided between them through the two-body center-of-mass kinematics.

Use this vibrational distribution as the baseline Mars exobase configuration. Keep it as a named configuration because the distribution may vary with altitude and atmospheric conditions. Do not assume that the fixed exobase distribution is valid globally without a sensitivity test.

Continue to use the user-specified reaction branching ratios of 26.5%, 47.3%, 20.4%, and 5.8%. Do not replace them with the older branching ratios used in Fox and Hać (1997).

## Fixed Dissociative Recombination Rate

Define the O2+ dissociative recombination rate coefficient as

\[
k(T_e)=
\begin{cases}
1.95\times10^{-7}
\left(\dfrac{300}{T_e}\right)^{0.70},
& T_e\leq1200\ \mathrm{K},\\
7.39\times10^{-8}
\left(\dfrac{1200}{T_e}\right)^{0.56},
& T_e>1200\ \mathrm{K}.
\end{cases}
\]

When \(T_e\) is in kelvin, \(k\) has units of \(\mathrm{cm^3\,s^{-1}}\).

The volumetric dissociative recombination rate is

\[
\alpha =
n_e n_{\mathrm{O_2^+}} k(T_e),
\]

with units of \(\mathrm{cm^{-3}\,s^{-1}}\) when both densities are in \(\mathrm{cm^{-3}}\).

Each reaction produces two oxygen atoms. Therefore, the total hot O atom production rate is

\[
Q_{\mathrm{hot\,O}}=2\alpha.
\]

Always distinguish among:

1. \(k(T_e)\), the rate coefficient
2. \(\alpha\), the number of dissociative recombination reactions per unit volume per unit time
3. \(Q_{\mathrm{hot\,O}}\), the number of hot O atoms produced per unit volume per unit time

Do not omit or double count the factor of 2.

## Units

MGITM number densities are stored in \(\mathrm{m^{-3}}\). Use SI units internally in Julia unless a specific validation calculation requires cgs units.

Use the following conversions explicitly:

\[
1\ \mathrm{cm^{-3}}=10^6\ \mathrm{m^{-3}},
\]

\[
1\ \mathrm{cm^3\,s^{-1}}=10^{-6}\ \mathrm{m^3\,s^{-1}}.
\]

Variable names or metadata must record units. Do not mix SI and cgs quantities in the same expression without an explicit conversion.

## MGITM Inputs

The MGITM files are stored in `MGITM`.

They include four seasons:

```text
Ls = 0, 90, 180, and 270 degrees
```

They include three solar activity levels:

```text
F070, F130, and F200
```

Relevant variables include:

```text
Tn, Ti, Te
nCO2, nO, nN2, nCO, nO2
nO2+, nO+, nCO2+, ne
UN, VN, WN
```

The MGITM altitude range is 98.75 to 251.25 km. Any extrapolation beyond the upper boundary must be documented and tested for sensitivity.

## Monte Carlo Requirements

The transport model must include:

1. Isotropic initial directions unless a different source distribution is scientifically justified.
2. Martian gravity between collisions.
3. Collisions of hot O with at least CO2, O, N2, and CO.
4. Energy loss and angular scattering in each collision.
5. Explicit lower, upper, thermalization, and escape boundary criteria.
6. Reproducible random number seeds.
7. Particle weights derived from the physical hot O production rate.
8. Residence-time tallies in spatial bins.

Calculating escape probability alone is insufficient for this project. The O corona density must be calculated from weighted particle residence times:

\[
n_{\mathrm{hot\,O},j}
=\frac{1}{V_j}\sum_i w_i\Delta t_{ij}.
\]

Store ballistic, escaping, and, when included, satellite populations separately.

## Fixed Hot O Collision Cross Sections

For the baseline hot O transport model, use the following total collision cross sections:

| Collision process | Cross section in cgs | Cross section in SI | Reference or assumption |
|---|---:|---:|---|
| O strikes CO2 | \(2.0\times10^{-14}\ \mathrm{cm^2}\) | \(2.0\times10^{-18}\ \mathrm{m^2}\) | Fox and Hać (2014) |
| O strikes O | \(0.6\times10^{-14}\ \mathrm{cm^2}\) | \(6.0\times10^{-19}\ \mathrm{m^2}\) | Kharchenko et al. (2000) |
| O strikes N2 | \(1.8\times10^{-14}\ \mathrm{cm^2}\) | \(1.8\times10^{-18}\ \mathrm{m^2}\) | Balakrishnan et al. (1998) |
| O strikes CO | \(1.8\times10^{-14}\ \mathrm{cm^2}\) | \(1.8\times10^{-18}\ \mathrm{m^2}\) | Assumed equal to O-N2 |

These are fixed, energy-independent total cross sections for the baseline model. Do not replace them with energy-dependent or differential cross sections unless the alternative model is implemented as a separately named configuration and the user explicitly approves the change.

The local inverse mean free path is

\[
\lambda^{-1}(\mathbf r)
=n_{\mathrm{CO_2}}\sigma_{\mathrm{O-CO_2}}
+n_{\mathrm O}\sigma_{\mathrm{O-O}}
+n_{\mathrm{N_2}}\sigma_{\mathrm{O-N_2}}
+n_{\mathrm{CO}}\sigma_{\mathrm{O-CO}}.
\]

The probability that a collision occurs with species \(s\) is

\[
P_s=
\frac{n_s\sigma_s}
{\sum_j n_j\sigma_j}.
\]

All densities and cross sections in these expressions must use a consistent unit system. The Julia baseline implementation must use the SI values in \(\mathrm{m^2}\) together with MGITM densities in \(\mathrm{m^{-3}}\).

## Development Sequence

Develop and validate the model in the following order:

1. One-dimensional spherically symmetric baseline model
2. Two-dimensional solar zenith angle model
3. Full three-dimensional MGITM-coupled model

Start with the subsolar profile from:

```text
MGITM/MGITM_LS000_F070_150901.dat
```

Do not run all 12 MGITM cases until the baseline model passes conservation, collision, convergence, and unit tests.

## Julia and Python Responsibilities

Julia is responsible for:

1. Physical source sampling
2. Particle trajectories
3. Collision sampling
4. Parallel Monte Carlo calculations
5. Residence-time and escape tallies
6. HDF5 or NetCDF output

Python is responsible for:

1. MGITM input inspection
2. Model output validation
3. Scientific plotting
4. Comparisons among seasons and solar activity levels
5. Monte Carlo convergence and uncertainty diagnostics

Prefer the Python environment:

```text
C:\Users\Win\.conda\envs\mars\python.exe
```

## Figure Requirements

Use Arial for all figure text except mathematical formulas. This includes axis labels, tick labels, legends, colorbar labels, titles, annotations, and panel labels.

Plots must show units explicitly and identify the MGITM season, solar activity level, solar zenith angle, and model configuration when applicable.

## Scientific Safety

Do not change reaction branches, branching ratios, reaction rate expressions, collision cross sections, boundary definitions, or unit conventions silently.

When a physical input is uncertain, preserve it as a named configuration parameter and document the tested range. Clearly distinguish literature values, assumptions, numerical approximations, and derived results.
