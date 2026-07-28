"""
Weighted particle ensembles and residence-time diagnostics.

This file builds a stratified source distribution, tracks primary and
secondary particles, records termination reasons, and converts residence time
into the spherical density in each altitude and energy bin, in m^-3.
"""
Base.@kwdef struct RahmatiMonteCarloConfig
    particles_per_source_altitude::Int = 10_000
    seed::Int = 20260727
    minimum_energy_eV::Float64 = 0.01
    maximum_altitude_m::Float64 = 110_000e3
    maximum_steps_per_particle::Int = 2_000_000
    maximum_total_particles::Int = 50_000_000
    source_altitudes_km::Vector{Float64} = collect(100.0:1.0:250.0)
    altitude_edges_km::Vector{Float64} = collect(100.0:5.0:1000.0)
    energy_edges_eV::Vector{Float64} = collect(range(0.01, 7.0; length=141))
end

struct HotOCoronaResult
    altitude_edges_km::Vector{Float64}
    energy_edges_eV::Vector{Float64}
    density_m3_per_bin::Matrix{Float64}
    upward_density_m3_per_bin::Matrix{Float64}
    downward_density_m3_per_bin::Matrix{Float64}
    source_altitudes_km::Vector{Float64}
    source_particle_weights_s1::Vector{Float64}
    particles_per_source_altitude::Int
    primary_particles::Int
    secondary_particles::Int
    total_source_rate_s1::Float64
    stop_counts::Dict{Symbol,Int}
end

"""Rahmati step rule: 0.1 mfp below 10 km, otherwise 1 km."""
rahmati_step_length(mean_free_path_m::Real) =
    mean_free_path_m < 10_000 ? 0.1 * mean_free_path_m : 1000.0

function _shell_edges_m(altitude_m)
    length(altitude_m) >= 2 || error("At least two source altitudes are required")
    edges = Vector{Float64}(undef, length(altitude_m) + 1)
    edges[2:end-1] .= (altitude_m[1:end-1] .+ altitude_m[2:end]) ./ 2
    edges[1] = altitude_m[1] - (edges[2] - altitude_m[1])
    edges[end] = altitude_m[end] + (altitude_m[end] - edges[end-1])
    edges
end

function _load_vibration(path)
    input = TOML.parsefile(path)["vibration"]
    probability = Float64.(input["fraction"])
    probability ./= sum(probability)
    probability, Float64(input["quantum_energy_eV"])
end

function _accumulate_residence!(
    residence_particles, particle_weight_s1, dt, altitude_m, energy_eV,
    altitude_edges_km, energy_edges_eV,
)
    ia = searchsortedlast(altitude_edges_km, altitude_m / 1000)
    ie = searchsortedlast(energy_edges_eV, energy_eV)
    if 1 <= ia < length(altitude_edges_km) &&
       1 <= ie < length(energy_edges_eV)
        residence_particles[ia, ie] += particle_weight_s1 * dt
    end
end

function _advance_gravity(position, velocity, ds)
    radius = _norm(position)
    speed = _norm(velocity)
    dt = ds / max(speed, 1.0)
    acceleration0 = _scale(-MARS_MU_M3_S2 / radius^3, position)
    position1 = _add(
        position,
        _add(_scale(dt, velocity), _scale(0.5dt^2, acceleration0)),
    )
    radius1 = _norm(position1)
    acceleration1 = _scale(-MARS_MU_M3_S2 / radius1^3, position1)
    velocity1 = _add(
        velocity, _scale(0.5dt, _add(acceleration0, acceleration1)),
    )
    position1, velocity1, dt
end

"""
Run the spherically extended Rahmati Monte Carlo model.

The nearest-subsolar MGITM profile is extended spherically. The model launches
the configured number of particles at every source altitude. A source particle
at altitude i represents Q_hotO(z_i) * V_shell,i / N_i physical O atoms per
second. The returned height-energy density uses a residence-time estimator and
is reported in m^-3 per energy bin.
"""
function run_hot_o_corona(
    atmosphere::AtmosphereProfile, targets, branches;
    chemistry_path::AbstractString,
    config::RahmatiMonteCarloConfig=RahmatiMonteCarloConfig(),
)
    config.particles_per_source_altitude > 0 ||
        error("particles_per_source_altitude must be positive")
    length(config.source_altitudes_km) >= 2 ||
        error("At least two source altitudes are required")
    issorted(config.source_altitudes_km) ||
        error("source_altitudes_km must be sorted")
    all(diff(config.source_altitudes_km) .> 0) ||
        error("source_altitudes_km must be strictly increasing")
    rng = Xoshiro(config.seed)
    source_altitudes_m = 1000 .* config.source_altitudes_km
    all(source_altitudes_m .> atmosphere.altitude_m[1]) ||
        error("Source altitudes must be above the atmospheric lower boundary")
    all(source_altitudes_m .< config.maximum_altitude_m) ||
        error("Source altitudes must be below maximum_altitude_m")
    source_edges_m = _shell_edges_m(source_altitudes_m)
    shell_volume_m3 = (4pi / 3) .* (
        (MARS_RADIUS_M .+ source_edges_m[2:end]).^3 .-
        (MARS_RADIUS_M .+ source_edges_m[1:end-1]).^3
    )
    source_states = [
        interpolate_profile(atmosphere, altitude_m)
        for altitude_m in source_altitudes_m
    ]
    source_q_m3_s1 = [
        hot_o_production_rate(
            state.density_m3[:e],
            state.density_m3[:O2p],
            state.Te_K,
        ) for state in source_states
    ]
    source_rate_s1 = source_q_m3_s1 .* shell_volume_m3
    source_particle_weights_s1 =
        source_rate_s1 ./ config.particles_per_source_altitude
    total_source_rate_s1 = sum(source_rate_s1)
    total_source_rate_s1 > 0 || error("The atmosphere has zero hot O source")
    primary_particles =
        length(source_altitudes_m) * config.particles_per_source_altitude
    primary_particles <= config.maximum_total_particles ||
        error("Primary particles exceed maximum_total_particles")
    vibration_probability, vibration_quantum_eV =
        _load_vibration(chemistry_path)

    queue = HotOParticle[]
    sizehint!(queue, min(primary_particles, 1_000_000))
    for iz in eachindex(source_altitudes_m)
        altitude_m = source_altitudes_m[iz]
        position = (MARS_RADIUS_M + altitude_m, 0.0, 0.0)
        state = source_states[iz]
        for _ in 1:config.particles_per_source_altitude
            push!(queue, sample_hot_o_source(
                rng, position, state.Te_K, state.Ti_K, branches;
                vibrational_probability=vibration_probability,
                vibrational_quantum_eV=vibration_quantum_eV,
                weight_s1=source_particle_weights_s1[iz],
            ))
        end
    end

    upward_residence_particles = zeros(
        length(config.altitude_edges_km) - 1,
        length(config.energy_edges_eV) - 1,
    )
    downward_residence_particles = zeros(
        length(config.altitude_edges_km) - 1,
        length(config.energy_edges_eV) - 1,
    )
    stops = Dict{Symbol,Int}()
    secondary_count = 0
    next_particle = 1
    while next_particle <= length(queue)
        particle = queue[next_particle]
        next_particle += 1
        reason = :maximum_steps
        for _ in 1:config.maximum_steps_per_particle
            altitude_m = _norm(particle.position_m) - MARS_RADIUS_M
            energy_eV = kinetic_energy_eV(particle.velocity_m_s)
            if energy_eV <= config.minimum_energy_eV
                reason = :thermalized
                break
            elseif altitude_m >= config.maximum_altitude_m
                reason = :escaped
                break
            elseif altitude_m <= atmosphere.altitude_m[1]
                reason = :lower_boundary
                break
            end

            local_state = interpolate_profile(atmosphere, altitude_m)
            kappa = collision_coefficient(
                targets, local_state.density_m3, energy_eV,
            )
            mfp = kappa > 0 ? inv(kappa) : Inf
            ds = rahmati_step_length(mfp)
            position1, velocity1, dt =
                _advance_gravity(particle.position_m, particle.velocity_m_s, ds)
            radial_velocity_m_s =
                _dot(particle.position_m, particle.velocity_m_s) /
                _norm(particle.position_m)
            directional_residence_particles =
                radial_velocity_m_s >= 0 ?
                upward_residence_particles : downward_residence_particles
            _accumulate_residence!(
                directional_residence_particles, particle.weight_s1, dt,
                altitude_m, energy_eV,
                config.altitude_edges_km, config.energy_edges_eV,
            )
            particle.position_m = position1
            particle.velocity_m_s = velocity1

            if kappa > 0 && rand(rng) < min(ds * kappa, 1.0)
                target = choose_collision_target(
                    rng, targets, local_state.density_m3, energy_eV,
                )
                theta_com = sample_scattering_angle(rng)
                projectile_after, target_after = elastic_collision(
                    particle.velocity_m_s, (0.0, 0.0, 0.0),
                    O_MASS_KG, target.mass_kg,
                    theta_com, sample_azimuth(rng),
                )
                particle.velocity_m_s = projectile_after
                particle.collisions += 1
                if target.species == :O &&
                   kinetic_energy_eV(target_after) > config.minimum_energy_eV
                    if length(queue) >= config.maximum_total_particles
                        reason = :maximum_particles
                        break
                    end
                    push!(queue, HotOParticle(
                        particle.position_m, target_after,
                        particle.weight_s1, true, 0,
                    ))
                    secondary_count += 1
                end
            end
        end
        stops[reason] = get(stops, reason, 0) + 1
        reason == :maximum_particles && break
    end

    altitude_radius_m = MARS_RADIUS_M .+
        1000 .* config.altitude_edges_km
    diagnostic_volume_m3 = (4pi / 3) .* (
        altitude_radius_m[2:end].^3 .- altitude_radius_m[1:end-1].^3
    )
    volume_column_m3 = reshape(diagnostic_volume_m3, :, 1)
    upward_density_m3_per_bin =
        upward_residence_particles ./ volume_column_m3
    downward_density_m3_per_bin =
        downward_residence_particles ./ volume_column_m3
    density_m3_per_bin =
        upward_density_m3_per_bin .+ downward_density_m3_per_bin
    HotOCoronaResult(
        copy(config.altitude_edges_km), copy(config.energy_edges_eV),
        density_m3_per_bin,
        upward_density_m3_per_bin,
        downward_density_m3_per_bin,
        copy(config.source_altitudes_km),
        source_particle_weights_s1,
        config.particles_per_source_altitude,
        primary_particles, secondary_count,
        total_source_rate_s1, stops,
    )
end

"""
Write total, upward, and downward residence-time density diagnostics.

The upward population has positive Mars-centered radial velocity and the
downward population has negative radial velocity. Densities are in m^-3 per
energy bin.
"""
function write_directional_corona_distribution(
    path::AbstractString, result::HotOCoronaResult,
)
    open(path, "w") do io
        println(
            io,
            "# altitude_km energy_eV total_density_m-3_per_bin ",
            "upward_density_m-3_per_bin downward_density_m-3_per_bin",
        )
        for ia in axes(result.density_m3_per_bin, 1)
            altitude = (
                result.altitude_edges_km[ia] +
                result.altitude_edges_km[ia + 1]
            ) / 2
            for ie in axes(result.density_m3_per_bin, 2)
                energy = (
                    result.energy_edges_eV[ie] +
                    result.energy_edges_eV[ie + 1]
                ) / 2
                println(
                    io,
                    altitude, ' ', energy, ' ',
                    result.density_m3_per_bin[ia, ie], ' ',
                    result.upward_density_m3_per_bin[ia, ie], ' ',
                    result.downward_density_m3_per_bin[ia, ie],
                )
            end
        end
    end
    path
end

"""Write a compact long-form text table for Python plotting."""
function write_corona_distribution(path::AbstractString, result::HotOCoronaResult)
    open(path, "w") do io
        println(io, "# altitude_km energy_eV density_m-3_per_bin")
        for ia in axes(result.density_m3_per_bin, 1)
            altitude = (result.altitude_edges_km[ia] +
                        result.altitude_edges_km[ia + 1]) / 2
            for ie in axes(result.density_m3_per_bin, 2)
                energy = (result.energy_edges_eV[ie] +
                          result.energy_edges_eV[ie + 1]) / 2
                println(io, altitude, ' ', energy, ' ',
                        result.density_m3_per_bin[ia, ie])
            end
        end
    end
    path
end
