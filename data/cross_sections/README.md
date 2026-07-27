# Hot O collision data

The total cross sections follow Rahmati's dissertation Table 1 and use
`sigma(E) = sigma(3 eV) (E / 3 eV)^(-0.2)`.

Scattering angles are not copied from MarsASPEN. MarsASPEN's lookup table is
an H and H+ laboratory-frame inverse CDF. MarsHotO samples the continuous
O-O center-of-mass differential cross section fitted by Rahmati:

`d sigma / d Omega = alpha sin(theta_COM / 2)^beta`,

with `alpha = 0.36e-16 cm2 sr-1` and `beta = -1.85`. The polar-angle PDF
includes the solid-angle Jacobian `sin(theta)`. The inverse CDF is evaluated
analytically in `src/scattering.jl`.
