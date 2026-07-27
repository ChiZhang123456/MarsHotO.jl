@inline function _maxwell_velocity(rng, temperature_K, mass_kg)
    sigma = sqrt(BOLTZMANN_J_K * temperature_K / mass_kg)
    (sigma * randn(rng), sigma * randn(rng), sigma * randn(rng))
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
    ve = _maxwell_velocity(rng, Te_K, ELECTRON_MASS_KG)
    vi = _maxwell_velocity(rng, Ti_K, O2P_MASS_KG)
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
    direction = _maxwell_velocity(rng, 1.0, 1.0)
    direction = _scale(1 / _norm(direction), direction)
    HotOParticle(position_m, _add(vcom, _scale(product_speed, direction)), 1.0, true, 0)
end
