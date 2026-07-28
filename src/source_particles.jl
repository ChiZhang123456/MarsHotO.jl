"""
Photochemical source-particle sampling.

Electron and O2+ velocities are sampled from normalized three-dimensional
Maxwellian distributions. Their default bulk velocities are zero. Reaction
branch, vibration, and product direction are sampled independently.
"""

"""
Most-probable thermal speed used by the Maxwellian convention.

This equals sqrt(2*kB*T/m), matching the TestParticle.jl Maxwellian
parameterization.
"""
function maxwellian_thermal_speed(temperature_K::Real, mass_kg::Real)
    temperature_K > 0 ||
        throw(DomainError(temperature_K, "Temperature must be positive"))
    mass_kg > 0 || throw(DomainError(mass_kg, "Mass must be positive"))
    sqrt(2BOLTZMANN_J_K * temperature_K / mass_kg)
end

"""
Normalized three-dimensional Maxwellian velocity PDF in SI units.

The returned probability density has units s^3 m^-3 and integrates to one
over all velocity space.
"""
function maxwellian_velocity_pdf(
    velocity_m_s,
    temperature_K::Real,
    mass_kg::Real;
    bulk_velocity_m_s=(0.0, 0.0, 0.0),
)
    length(velocity_m_s) == 3 ||
        throw(DimensionMismatch("Velocity must have three components"))
    length(bulk_velocity_m_s) == 3 ||
        throw(DimensionMismatch("Bulk velocity must have three components"))
    temperature_K > 0 ||
        throw(DomainError(temperature_K, "Temperature must be positive"))
    mass_kg > 0 || throw(DomainError(mass_kg, "Mass must be positive"))

    variance_m2_s2 = BOLTZMANN_J_K * temperature_K / mass_kg
    relative_speed_squared = sum(
        (velocity_m_s[i] - bulk_velocity_m_s[i])^2 for i in 1:3
    )
    exp(-relative_speed_squared / (2variance_m2_s2)) /
        (2pi * variance_m2_s2)^(3 / 2)
end

"""
Sample a three-dimensional Maxwellian velocity.

Each Cartesian component is Gaussian with mean equal to the corresponding
bulk velocity and variance kB*T/m. The resulting direction is isotropic when
the bulk velocity is zero.
"""
function sample_maxwellian_velocity(
    rng,
    temperature_K::Real,
    mass_kg::Real;
    bulk_velocity_m_s=(0.0, 0.0, 0.0),
)
    length(bulk_velocity_m_s) == 3 ||
        throw(DimensionMismatch("Bulk velocity must have three components"))
    thermal_speed = maxwellian_thermal_speed(temperature_K, mass_kg)
    component_sigma = thermal_speed / sqrt(2)
    ntuple(
        i -> Float64(bulk_velocity_m_s[i]) + component_sigma * randn(rng),
        3,
    )
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
    # Both reactants use normalized Maxwellians with zero bulk velocity.
    ve = sample_maxwellian_velocity(rng, Te_K, ELECTRON_MASS_KG)
    vi = sample_maxwellian_velocity(rng, Ti_K, O2P_MASS_KG)
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
