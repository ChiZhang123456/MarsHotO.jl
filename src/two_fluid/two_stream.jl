function _validate_two_stream_config(config::TwoStreamConfig)
    config.altitude_min_m < config.altitude_max_m ||
        error("The lower two-stream boundary must be below the upper boundary")
    config.altitude_step_m > 0 || error("Altitude step must be positive")
    0 < config.mean_pitch_cosine <= 1 ||
        error("Mean pitch-angle cosine must be in (0, 1]")
    config.redistribution_samples > 0 ||
        error("Redistribution sample count must be positive")
    length(config.energy_edges_eV) >= 2 ||
        error("At least two energy edges are required")
    issorted(config.energy_edges_eV) ||
        error("Energy edges must be sorted")
    all(diff(config.energy_edges_eV) .> 0) ||
        error("Energy edges must be strictly increasing")
end

@inline _energy_bin(energy_eV, edges) =
    energy_eV < first(edges) || energy_eV >= last(edges) ? 0 :
    searchsortedlast(edges, energy_eV)

"""
Construct collision redistribution probabilities for the two-stream solver.

The complete MarsASPEN inverse CDF is interpreted as an empirical COM angle
distribution and the full Rahmati total cross sections are used. No
small-angle cutoff or cross-section rescaling is applied. The incoming stream
is represented by |mu| = mean_pitch_cosine.
"""
function build_two_stream_redistribution(
    targets,
    energy_edges_eV;
    mean_pitch_cosine=0.5,
    samples=20_000,
    seed=73,
    distribution=default_scattering_angle_distribution(),
)
    0 < mean_pitch_cosine <= 1 ||
        error("Mean pitch-angle cosine must be in (0, 1]")
    samples > 0 || error("Sample count must be positive")
    centers = (energy_edges_eV[1:end-1] .+ energy_edges_eV[2:end]) ./ 2
    nt, ne = length(targets), length(centers)
    same = zeros(nt, ne, ne)
    reverse = zeros(nt, ne, ne)
    secondary_same = zeros(nt, ne, ne)
    secondary_reverse = zeros(nt, ne, ne)
    incoming_direction = (
        sqrt(max(1 - mean_pitch_cosine^2, 0.0)), 0.0, mean_pitch_cosine,
    )
    rng = Xoshiro(seed)

    for (itarget, target) in pairs(targets)
        is_oxygen = target.species == :O
        for (isource, energy_eV) in pairs(centers)
            speed = sqrt(2energy_eV * EV_J / O_MASS_KG)
            incoming = _scale(speed, incoming_direction)
            for _ in 1:samples
                theta = sample_scattering_angle(
                    rng, distribution,
                )
                projectile, recoil = elastic_collision(
                    incoming, (0.0, 0.0, 0.0),
                    O_MASS_KG, target.mass_kg,
                    theta, sample_azimuth(rng),
                )
                projectile_energy =
                    0.5O_MASS_KG * _dot(projectile, projectile) / EV_J
                idestination = _energy_bin(projectile_energy, energy_edges_eV)
                if idestination > 0
                    storage = projectile[3] >= 0 ? same : reverse
                    storage[itarget, isource, idestination] += 1 / samples
                end
                if is_oxygen
                    recoil_energy =
                        0.5target.mass_kg * _dot(recoil, recoil) / EV_J
                    isecondary = _energy_bin(recoil_energy, energy_edges_eV)
                    if isecondary > 0
                        storage = recoil[3] >= 0 ?
                            secondary_same : secondary_reverse
                        storage[itarget, isource, isecondary] += 1 / samples
                    end
                end
            end
        end
    end
    TwoStreamRedistribution(same, reverse, secondary_same, secondary_reverse)
end

function _source_spectrum(profile_point, branches, chemistry_path, edges)
    chemistry = TOML.parsefile(chemistry_path)
    vibration = chemistry["vibration"]
    vibration_levels = vibration["level"]
    vibration_probability = vibration["fraction"]
    vibration_probability ./= sum(vibration_probability)
    quantum_eV = vibration["quantum_energy_eV"]
    source = zeros(length(edges) - 1)
    densities = profile_point.density_m3
    total_rate = hot_o_production_rate(
        densities[:e], densities[:O2p], profile_point.Te_K,
    )
    for branch in branches
        for (level, probability) in zip(
            vibration_levels, vibration_probability,
        )
            # The tabulated release and vibrational energies are shared by
            # two equal-mass O products.
            product_energy = (branch.release_energy_eV + level * quantum_eV) / 2
            index = _energy_bin(product_energy, edges)
            index > 0 && (source[index] +=
                total_rate * branch.probability * probability)
        end
    end
    source
end

function _top_scale_height(profile, top_m)
    dz = min(1e3, max((last(profile.altitude_m) - first(profile.altitude_m)) / 20, 1.0))
    lower = interpolate_profile(profile, top_m - dz).density_m3[:O2p]
    upper = interpolate_profile(profile, top_m).density_m3[:O2p]
    slope = (log(max(upper, floatmin(Float64))) -
             log(max(lower, floatmin(Float64)))) / dz
    slope < 0 ? -1 / slope : 70e3
end

@inline function _analytic_flux_step(flux, source, kappa, dz, mu)
    if kappa <= 0
        return flux + source * dz / mu
    end
    attenuation = exp(-kappa * dz / mu)
    flux * attenuation + source / kappa * (1 - attenuation)
end

function _collision_sources!(
    upward_source,
    downward_source,
    upward_flux,
    downward_flux,
    density,
    targets,
    redistribution,
    energy_centers,
)
    fill!(upward_source, 0.0)
    fill!(downward_source, 0.0)
    ne = length(energy_centers)
    for (itarget, target) in pairs(targets)
        target_density = get(density, target.species, 0.0)
        target_density <= 0 && continue
        for source_index in 1:ne
            rate_coefficient = target_density *
                total_cross_section(target, energy_centers[source_index])
            up_collision = rate_coefficient * upward_flux[source_index]
            down_collision = rate_coefficient * downward_flux[source_index]
            for destination_index in 1:source_index
                same_probability =
                    redistribution.same_stream[
                        itarget, source_index, destination_index
                    ] +
                    redistribution.secondary_same_stream[
                        itarget, source_index, destination_index
                    ]
                reverse_probability =
                    redistribution.reverse_stream[
                        itarget, source_index, destination_index
                    ] +
                    redistribution.secondary_reverse_stream[
                        itarget, source_index, destination_index
                    ]
                upward_source[destination_index] +=
                    same_probability * up_collision +
                    reverse_probability * down_collision
                downward_source[destination_index] +=
                    same_probability * down_collision +
                    reverse_probability * up_collision
            end
        end
    end
end

"""
Solve the Rahmati (2016), section 2.2.1, one-dimensional two-stream equations.

Fluxes and production are integrated over each energy bin. Atmosphere inputs
remain in SI units. The lower boundary is isotropic, while the upper boundary
reflects gravitationally bound atoms and lets escaping atoms leave.
"""
function run_two_stream(
    profile::AtmosphereProfile,
    targets,
    branches;
    chemistry_path=normpath(joinpath(
        @__DIR__, "..", "..", "data", "chemistry",
        "o2plus_dissociative_recombination.toml",
    )),
    config=TwoStreamConfig(),
    redistribution=nothing,
)
    _validate_two_stream_config(config)
    altitude = collect(
        config.altitude_min_m:config.altitude_step_m:config.altitude_max_m,
    )
    isapprox(last(altitude), config.altitude_max_m; atol=1e-8) ||
        error("Altitude range must be an integer number of altitude steps")
    edges = config.energy_edges_eV
    centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
    nz, ne = length(altitude), length(centers)
    redistribution = isnothing(redistribution) ?
        build_two_stream_redistribution(
            targets, edges;
            mean_pitch_cosine=config.mean_pitch_cosine,
            samples=config.redistribution_samples,
            seed=config.redistribution_seed,
        ) : redistribution

    density = Vector{Dict{Symbol,Float64}}(undef, nz)
    primary = zeros(nz, ne)
    kappa = zeros(nz, ne)
    for iz in 1:nz
        point = interpolate_profile(profile, altitude[iz])
        density[iz] = point.density_m3
        primary[iz, :] .= _source_spectrum(
            point, branches, chemistry_path, edges,
        )
        for ie in 1:ne
            kappa[iz, ie] =
                collision_coefficient(targets, density[iz], centers[ie])
        end
    end

    upward = zeros(nz, ne)
    downward = zeros(nz, ne)
    previous_upward = similar(upward)
    previous_downward = similar(downward)
    cascade_up = zeros(nz, ne)
    cascade_down = zeros(nz, ne)
    local_up = zeros(ne)
    local_down = zeros(ne)
    scale_height = isnothing(config.top_scale_height_m) ?
        _top_scale_height(profile, last(altitude)) :
        config.top_scale_height_m
    escape_energy =
        MARS_MU_M3_S2 * O_MASS_KG /
        ((MARS_RADIUS_M + last(altitude)) * EV_J)
    converged = false
    iterations = 0

    for iteration in 1:config.maximum_iterations
        iterations = iteration
        copyto!(previous_upward, upward)
        copyto!(previous_downward, downward)
        for iz in 1:nz
            _collision_sources!(
                local_up, local_down,
                view(previous_upward, iz, :),
                view(previous_downward, iz, :),
                density[iz], targets, redistribution, centers,
            )
            cascade_up[iz, :] .= local_up
            cascade_down[iz, :] .= local_down
        end

        for ie in ne:-1:1
            top_external = primary[end, ie] * scale_height / 4
            reflected = centers[ie] < escape_energy ? upward[end, ie] : 0.0
            downward[end, ie] = top_external + reflected
            for iz in nz:-1:2
                source = primary[iz, ie] / 2 + cascade_down[iz, ie]
                downward[iz - 1, ie] = _analytic_flux_step(
                    downward[iz, ie], source, kappa[iz, ie],
                    config.altitude_step_m, config.mean_pitch_cosine,
                )
            end
            upward[1, ie] = downward[1, ie]
            for iz in 1:nz-1
                source = primary[iz, ie] / 2 + cascade_up[iz, ie]
                upward[iz + 1, ie] = _analytic_flux_step(
                    upward[iz, ie], source, kappa[iz, ie],
                    config.altitude_step_m, config.mean_pitch_cosine,
                )
            end
        end
        change = max(
            maximum(abs.(upward .- previous_upward)),
            maximum(abs.(downward .- previous_downward)),
        )
        magnitude = max(maximum(upward), maximum(downward), 1.0)
        if change / magnitude <= config.relative_tolerance
            converged = true
            break
        end
    end

    escaping = centers .>= escape_energy
    escape_flux = sum(view(upward, nz, escaping))
    TwoStreamResult(
        altitude, copy(edges), centers, upward, downward, primary,
        escape_energy, escape_flux, iterations, converged,
    )
end

function write_two_stream_flux(path::AbstractString, result::TwoStreamResult)
    open(path, "w") do io
        println(io, "altitude_km,energy_eV,upward_flux_m-2_s-1,downward_flux_m-2_s-1,primary_production_m-3_s-1")
        for iz in eachindex(result.altitude_m)
            for ie in eachindex(result.energy_centers_eV)
                println(
                    io,
                    result.altitude_m[iz] / 1e3, ",",
                    result.energy_centers_eV[ie], ",",
                    result.upward_flux_m2_s1[iz, ie], ",",
                    result.downward_flux_m2_s1[iz, ie], ",",
                    result.primary_production_m3_s1[iz, ie],
                )
            end
        end
    end
    path
end
