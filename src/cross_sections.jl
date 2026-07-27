function load_collision_targets(path::AbstractString)
    input = TOML.parsefile(path)
    [CollisionTarget(
        Symbol(item["species"]),
        item["mass_amu"] * AMU_KG,
        item["sigma_3eV_cm2"] * 1e-4,
    ) for item in input["target"]]
end

"""Rahmati energy-dependent total cross section, in m^2."""
function total_cross_section(target::CollisionTarget, energy_eV::Real)
    energy_eV > 0 || return 0.0
    target.sigma_3eV_m2 * (energy_eV / 3)^(-0.2)
end

function collision_coefficient(targets, density_m3, energy_eV;
                               theta_min_rad=0.0)
    angular_fraction = angular_cross_section_fraction(theta_min_rad)
    sum(get(density_m3, target.species, 0.0) *
        total_cross_section(target, energy_eV) * angular_fraction
        for target in targets)
end

function choose_collision_target(rng, targets, density_m3, energy_eV;
                                 theta_min_rad=0.0)
    angular_fraction = angular_cross_section_fraction(theta_min_rad)
    weights = [get(density_m3, target.species, 0.0) *
               total_cross_section(target, energy_eV) * angular_fraction
               for target in targets]
    total = sum(weights)
    total > 0 || return nothing
    threshold = rand(rng) * total
    accumulator = 0.0
    for (target, weight) in zip(targets, weights)
        accumulator += weight
        threshold <= accumulator && return target
    end
    last(targets)
end
