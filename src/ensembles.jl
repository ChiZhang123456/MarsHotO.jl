Base.@kwdef struct RahmatiMonteCarloConfig
    primary_particles::Int = 10_000
    seed::Int = 20260727
    scattering_angle_path::String = DEFAULT_SCATTERING_ANGLE_PATH
    minimum_energy_eV::Float64 = 0.01
    maximum_altitude_m::Float64 = 110_000e3
    maximum_steps_per_particle::Int = 2_000_000
    maximum_total_particles::Int = 50_000_000
    source_minimum_altitude_m::Float64 = 150e3
    source_maximum_altitude_m::Float64 = 500e3
    altitude_edges_km::Vector{Float64} = collect(100.0:5.0:1000.0)
    energy_edges_eV::Vector{Float64} = collect(range(0.01, 7.0; length=141))
end

struct HotOCoronaResult
    altitude_edges_km::Vector{Float64}
    energy_edges_eV::Vector{Float64}
    density_cm3_eV1::Matrix{Float64}
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

function _sample_weighted_index(rng, cumulative_weight)
    searchsortedfirst(cumulative_weight, rand(rng) * cumulative_weight[end])
end

function _load_vibration(path)
    input = TOML.parsefile(path)["vibration"]
    probability = Float64.(input["fraction"])
    probability ./= sum(probability)
    probability, Float64(input["quantum_energy_eV"])
end

function _accumulate_residence!(
    residence_s, particle_rate_s1, dt, altitude_m, energy_eV,
    altitude_edges_km, energy_edges_eV,
)
    ia = searchsortedlast(altitude_edges_km, altitude_m / 1000)
    ie = searchsortedlast(energy_edges_eV, energy_eV)
    if 1 <= ia < length(altitude_edges_km) &&
       1 <= ie < length(energy_edges_eV)
        residence_s[ia, ie] += particle_rate_s1 * dt
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
Run the spherical-column Rahmati Monte Carlo model.

The nearest-subsolar MGITM profile is extended spherically. Source particles
are sampled in proportion to Q_hotO(z) times shell volume. The returned
height-energy density uses a residence-time estimator.
"""
function run_hot_o_corona(
    atmosphere::AtmosphereProfile, targets, branches;
    chemistry_path::AbstractString,
    config::RahmatiMonteCarloConfig=RahmatiMonteCarloConfig(),
)
    config.primary_particles > 0 || error("primary_particles must be positive")
    rng = Xoshiro(config.seed)
    source_edges_m = _shell_edges_m(atmosphere.altitude_m)
    shell_volume_m3 = (4pi / 3) .* (
        (MARS_RADIUS_M .+ source_edges_m[2:end]).^3 .-
        (MARS_RADIUS_M .+ source_edges_m[1:end-1]).^3
    )
    source_q_m3_s1 = [
        hot_o_production_rate(
            atmosphere.density_m3[:e][i],
            atmosphere.density_m3[:O2p][i],
            atmosphere.Te_K[i],
        ) for i in eachindex(atmosphere.altitude_m)
    ]
    source_weight = source_q_m3_s1 .* shell_volume_m3
    source_weight .*= (
        (atmosphere.altitude_m .>= config.source_minimum_altitude_m) .&
        (atmosphere.altitude_m .<= config.source_maximum_altitude_m)
    )
    cumulative_source = cumsum(source_weight)
    total_source_rate_s1 = cumulative_source[end]
    total_source_rate_s1 > 0 || error("The atmosphere has zero hot O source")
    macro_rate_s1 = total_source_rate_s1 / config.primary_particles
    vibration_probability, vibration_quantum_eV =
        _load_vibration(chemistry_path)
    scattering_distribution = load_scattering_angle_distribution(
        config.scattering_angle_path,
    )

    queue = HotOParticle[]
    sizehint!(queue, min(config.primary_particles, 1_000_000))
    for _ in 1:config.primary_particles
        iz = _sample_weighted_index(rng, cumulative_source)
        zlo, zhi = source_edges_m[iz], source_edges_m[iz + 1]
        altitude_m = zlo + rand(rng) * (zhi - zlo)
        position = (MARS_RADIUS_M + altitude_m, 0.0, 0.0)
        push!(queue, sample_hot_o_source(
            rng, position, atmosphere.Te_K[iz], atmosphere.Ti_K[iz], branches;
            vibrational_probability=vibration_probability,
            vibrational_quantum_eV=vibration_quantum_eV,
        ))
    end

    residence = zeros(
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
                reason = :upper_boundary
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
            _accumulate_residence!(
                residence, macro_rate_s1, dt, altitude_m, energy_eV,
                config.altitude_edges_km, config.energy_edges_eV,
            )
            particle.position_m = position1
            particle.velocity_m_s = velocity1

            if kappa > 0 && rand(rng) < min(ds * kappa, 1.0)
                target = choose_collision_target(
                    rng, targets, local_state.density_m3, energy_eV,
                )
                theta_lab = sample_scattering_angle(
                    rng, scattering_distribution,
                    O_MASS_KG, target.mass_kg,
                )
                projectile_after, target_after = elastic_collision_lab(
                    particle.velocity_m_s,
                    O_MASS_KG, target.mass_kg,
                    theta_lab, sample_azimuth(rng),
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
                        particle.weight, true, 0,
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
    energy_width_eV = diff(config.energy_edges_eV)
    density_m3_eV1 = residence ./
        (diagnostic_volume_m3 .* transpose(energy_width_eV))
    HotOCoronaResult(
        copy(config.altitude_edges_km), copy(config.energy_edges_eV),
        density_m3_eV1 .* 1e-6,
        config.primary_particles, secondary_count,
        total_source_rate_s1, stops,
    )
end

"""Write a compact long-form text table for Python plotting."""
function write_corona_distribution(path::AbstractString, result::HotOCoronaResult)
    open(path, "w") do io
        println(io, "# altitude_km energy_eV density_cm-3_eV-1")
        for ia in axes(result.density_cm3_eV1, 1)
            altitude = (result.altitude_edges_km[ia] +
                        result.altitude_edges_km[ia + 1]) / 2
            for ie in axes(result.density_cm3_eV1, 2)
                energy = (result.energy_edges_eV[ie] +
                          result.energy_edges_eV[ie + 1]) / 2
                println(io, altitude, ' ', energy, ' ',
                        result.density_cm3_eV1[ia, ie])
            end
        end
    end
    path
end
