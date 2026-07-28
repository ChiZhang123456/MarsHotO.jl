# Hot O collision data

The total cross sections follow the values used in the Rahmati transport
model and use

```math
\sigma(E)=\sigma(3\ \mathrm{eV})(E/3\ \mathrm{eV})^{-0.2}.
```

The active scattering model is the analytical Rahmati fit to the Kharchenko
O and O differential cross section. It is implemented in
`src/shared/scattering.jl`, is sampled in COM by inverse transform, and uses
the complete angular interval without a 10 degree cutoff.

`scattering_angle_distribution.txt` is retained only as a MarsASPEN reference
table. It is not used by the active MarsHotO runtime.

Reference:

Kallio, E., and Barabash, S. (2001), Atmospheric effects of precipitating
energetic hydrogen atoms on the Martian atmosphere, Journal of Geophysical
Research: Space Physics, 106(A1), 165 to 177,
https://doi.org/10.1029/2000JA002003.
