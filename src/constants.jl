"""
Physical constants used throughout MarsHotO.

All internal calculations use SI units. Energies exposed by public particle
and diagnostic interfaces are expressed in electronvolts.
"""
const EV_J = 1.602176634e-19
const BOLTZMANN_J_K = 1.380649e-23
const AMU_KG = 1.66053906892e-27
const ELECTRON_MASS_KG = 9.1093837139e-31
const O_MASS_KG = 15.999 * AMU_KG
const O2P_MASS_KG = 2O_MASS_KG
const MARS_RADIUS_M = 3389.5e3
const MARS_MU_M3_S2 = 4.282837e13
