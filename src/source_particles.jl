"""
Photochemical source-particle sampling.

Electron and O2+ bulk velocities are zero. Their nonnegative thermal energies
use the configured zero-mode half-normal approximation with sigma_E=kB*T.
Thermal directions, reaction branch, vibration, and product direction are
sampled independently.
"""

"""
Sample nonnegative thermal energy from the configured zero-mode half-normal
model, with energy width sigma_E = kB*T.
"""
function _sample_thermal_energy_J(rng, temperature_K)
    temperature_K > 0 ||
        throw(DomainError(temperature_K, "Temperature must be positive"))
    abs(randn(rng)) * BOLTZMANN_J_K * temperature_K
end

"""Sample an isotropic unit vector independently of particle energy."""
@inline function _sample_isotropic_direction(rng)
    mu = 2rand(rng) - 1
    azimuth = 2pi * rand(rng)
    transverse = sqrt(max(1 - mu^2, 0.0))
    (
        transverse * cos(azimuth),
        transverse * sin(azimuth),
        mu,
    )
end

"""
Sample a zero-bulk-velocity particle by drawing its kinetic energy from the
configured zero-mode half-normal model and its direction isotropically.
"""
function _sample_thermal_velocity(rng, temperature_K, mass_kg)
    mass_kg > 0 || throw(DomainError(mass_kg, "Mass must be positive"))
    energy_J = _sample_thermal_energy_J(rng, temperature_K)
    speed_m_s = sqrt(2energy_J / mass_kg)
    _scale(speed_m_s, _sample_isotropic_direction(rng))
end

function _weighted_index(rng, probability)
    u, accumulator = rand(rng), 0.0
    for (i, p) in pairs(probability)
        accumulator += p
        u <= accumulator && return i
    end
    lastindex(probability)
end

"""
Sample one of the two O products of O2+ dissociative recombination. Electron
and O2+ thermal velocities, reaction branch, vibration, and isotropic product
direction are sampled explicitly.
"""
function sample_hot_o_source(rng, position_m, Te_K, Ti_K, branches;
                             vibrational_probability=[1.0],
                             vibrational_quantum_eV=0.23)
    # Both reactants have zero prescribed bulk velocity. Their thermal
    # energies and isotropic directions are sampled independently.
    ve = _sample_thermal_velocity(rng, Te_K, ELECTRON_MASS_KG)
    vi = _sample_thermal_velocity(rng, Ti_K, O2P_MASS_KG)
    total_mass = ELECTRON_MASS_KG + O2P_MASS_KG
    vcom = _scale(1 / total_mass,
        _add(_scale(ELECTRON_MASS_KG, ve), _scale(O2P_MASS_KG, vi)))
    relative = _subtract(ve, vi)
    reduced_mass = ELECTRON_MASS_KG * O2P_MASS_KG / total_mass
    relative_energy_eV = 0.5 * reduced_mass * _dot(relative, relative) / EV_J
    branch = branches[_weighted_index(rng, getfield.(branches, :probability))]
    vibrational_level = _weighted_index(rng, vibrational_probability) - 1
    available_eV = branch.release_energy_eV + relative_energy_eV +
                   vibrational_level * vibrational_quantum_eV
    product_speed = sqrt(available_eV * EV_J / O_MASS_KG)
    direction = _sample_isotropic_direction(rng)
    HotOParticle(position_m, _add(vcom, _scale(product_speed, direction)), 1.0, true, 0)
end
