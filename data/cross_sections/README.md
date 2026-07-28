# Hot O collision data

The total cross sections follow the values used in the Rahmati transport
model and use

```math
\sigma(E)=\sigma(3\ \mathrm{eV})(E/3\ \mathrm{eV})^{-0.2}.
```

The scattering angles are read from
`scattering_angle_distribution.txt`. This file is copied from MarsASPEN and
is an inverse CDF mapping from a uniform random number to the projectile
scattering angle. The source table labels the angle as LAB and was digitized
from Figure 2 of Kallio and Barabash (2001). MarsHotO intentionally treats
the tabulated angle values as an empirical COM distribution before applying
two-body kinematics.

The complete tabulated angular range is used. MarsHotO applies no additional
minimum angle cutoff and therefore uses the full total collision cross
section.

Reference:

Kallio, E., and Barabash, S. (2001), Atmospheric effects of precipitating
energetic hydrogen atoms on the Martian atmosphere, Journal of Geophysical
Research: Space Physics, 106(A1), 165 to 177,
https://doi.org/10.1029/2000JA002003.
