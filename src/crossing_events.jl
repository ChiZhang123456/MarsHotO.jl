"""
Particle-level altitude-crossing event output.

The binary stream is intended as the primary raw output for directional flux
diagnostics. Each record stores one particle birth, altitude-surface crossing,
or terminal event. Downstream analysis can reconstruct upward and downward
energy spectra without retaining every numerical integration step.
"""

const HOT_O_EVENT_MAGIC = codeunits("MHOTE001")
const HOT_O_EVENT_HEADER_BYTES = 64

const EVENT_BIRTH = Int8(1)
const EVENT_CROSSING = Int8(2)
const EVENT_EXIT_LOWER = Int8(3)
const EVENT_EXIT_UPPER = Int8(4)
const EVENT_THERMALIZED = Int8(5)
const EVENT_MAXIMUM_STEPS = Int8(6)
const EVENT_MAXIMUM_PARTICLES = Int8(7)

Base.@kwdef struct HotOCrossingConfig
    particles_per_source_altitude::Int = 10_000
    seed::Int = 20260810
    minimum_energy_eV::Float64 = 0.01
    domain_minimum_altitude_km::Float64 = 100.0
    domain_maximum_altitude_km::Float64 = 2000.0
    maximum_steps_per_particle::Int = 2_000_000
    maximum_total_particles::Int = 50_000_000
    source_altitudes_km::Vector{Float64} = collect(100.0:1.0:250.0)
    crossing_altitudes_km::Vector{Float64} = collect(100.0:10.0:2000.0)
end

"""
Fixed-width native little-endian event record.

The final four bytes are structure padding on 64-bit Julia. Python readers
must use the documented field offsets and an item size equal to
`sizeof(HotOEventRecord)`.
"""
struct HotOEventRecord
    particle_id::Int64
    parent_id::Int64
    weight_s1::Float64
    time_s::Float64
    altitude_km::Float64
    velocity_x_m_s::Float64
    velocity_y_m_s::Float64
    velocity_z_m_s::Float64
    radial_velocity_m_s::Float64
    event_index::Int32
    collisions::Int32
    surface_index::Int16
    event_code::Int8
    direction::Int8
end

struct HotOCrossingRunResult
    event_records::Int64
    primary_particles::Int
    secondary_particles::Int
    total_source_rate_s1::Float64
    source_particle_weights_s1::Vector{Float64}
    stop_counts::Dict{Symbol,Int}
    output_path::String
end

mutable struct _EventTrackedParticle
    particle::HotOParticle
    particle_id::Int64
    parent_id::Int64
    time_s::Float64
    event_index::Int32
end

@inline function _radial_velocity(position_m, velocity_m_s)
    _dot(position_m, velocity_m_s) / _norm(position_m)
end

@inline _direction_code(radial_velocity_m_s) =
    radial_velocity_m_s >= 0 ? Int8(1) : Int8(-1)

function _event_record(
    tracked::_EventTrackedParticle,
    event_code::Int8;
    surface_index::Int16=Int16(-1),
    altitude_km=nothing,
    position_m=tracked.particle.position_m,
    velocity_m_s=tracked.particle.velocity_m_s,
)
    radial_velocity_m_s = _radial_velocity(position_m, velocity_m_s)
    event_altitude_km = isnothing(altitude_km) ?
        (_norm(position_m) - MARS_RADIUS_M) / 1000 : altitude_km
    HotOEventRecord(
        tracked.particle_id,
        tracked.parent_id,
        tracked.particle.weight_s1,
        tracked.time_s,
        event_altitude_km,
        velocity_m_s[1],
        velocity_m_s[2],
        velocity_m_s[3],
        radial_velocity_m_s,
        tracked.event_index,
        Int32(tracked.particle.collisions),
        surface_index,
        event_code,
        _direction_code(radial_velocity_m_s),
    )
end

function _write_event!(
    io, event_count::Base.RefValue{Int64},
    tracked::_EventTrackedParticle, event_code::Int8; kwargs...,
)
    record = _event_record(tracked, event_code; kwargs...)
    write(io, Ref(record))
    event_count[] += 1
    nothing
end

function _interpolate_crossing_state(
    position0, velocity0, position1, velocity1,
    altitude0_km, altitude1_km, crossing_altitude_km,
)
    fraction = (
        crossing_altitude_km - altitude0_km
    ) / (altitude1_km - altitude0_km)
    position = _add(
        position0,
        _scale(fraction, _subtract(position1, position0)),
    )
    velocity = _add(
        velocity0,
        _scale(fraction, _subtract(velocity1, velocity0)),
    )
    position, velocity, fraction
end

function _write_crossings!(
    io, event_count, tracked,
    position0, velocity0, position1, velocity1,
    elapsed_step_s, crossing_altitudes_km,
)
    altitude0_km = (_norm(position0) - MARS_RADIUS_M) / 1000
    altitude1_km = (_norm(position1) - MARS_RADIUS_M) / 1000
    altitude1_km == altitude0_km && return

    indices = if altitude1_km > altitude0_km
        first_index = searchsortedlast(crossing_altitudes_km, altitude0_km) + 1
        last_index = searchsortedlast(crossing_altitudes_km, altitude1_km)
        first_index:last_index
    else
        first_index = searchsortedfirst(crossing_altitudes_km, altitude1_km)
        last_index = searchsortedfirst(crossing_altitudes_km, altitude0_km) - 1
        last_index:-1:first_index
    end

    time_at_step_start = tracked.time_s - elapsed_step_s
    for surface_index in indices
        1 <= surface_index <= length(crossing_altitudes_km) || continue
        crossing_altitude_km = crossing_altitudes_km[surface_index]
        position, velocity, fraction = _interpolate_crossing_state(
            position0, velocity0, position1, velocity1,
            altitude0_km, altitude1_km, crossing_altitude_km,
        )
        tracked.event_index += 1
        original_time_s = tracked.time_s
        tracked.time_s = time_at_step_start + fraction * elapsed_step_s
        _write_event!(
            io, event_count, tracked, EVENT_CROSSING;
            surface_index=Int16(surface_index),
            altitude_km=crossing_altitude_km,
            position_m=position,
            velocity_m_s=velocity,
        )
        tracked.time_s = original_time_s
    end
    nothing
end

function _write_event_header(io)
    write(io, HOT_O_EVENT_MAGIC)
    write(io, UInt32(1))
    write(io, UInt32(sizeof(HotOEventRecord)))
    write(io, UInt64(0))
    write(io, UInt64(0))
    write(io, UInt64(0))
    write(io, UInt64(0))
    write(io, zeros(UInt8, HOT_O_EVENT_HEADER_BYTES - position(io)))
    nothing
end

function _finalize_event_header(
    io, event_records, primary_particles, secondary_particles,
)
    final_position = position(io)
    seek(io, 16)
    write(io, UInt64(event_records))
    write(io, UInt64(primary_particles))
    write(io, UInt64(secondary_particles))
    write(io, UInt64(primary_particles + secondary_particles))
    seek(io, final_position)
    nothing
end

"""
Run a weighted hot O ensemble and stream particle-level crossing events.

The computational domain is closed at the configured lower and upper
altitudes. A particle is terminated immediately after crossing either domain
boundary. Event records are written incrementally, so memory use does not
scale with the number of recorded crossings.
"""
function run_hot_o_crossing_events(
    atmosphere::AtmosphereProfile, targets, branches;
    chemistry_path::AbstractString,
    output_path::AbstractString,
    config::HotOCrossingConfig=HotOCrossingConfig(),
)
    config.particles_per_source_altitude > 0 ||
        error("particles_per_source_altitude must be positive")
    issorted(config.source_altitudes_km) ||
        error("source_altitudes_km must be sorted")
    issorted(config.crossing_altitudes_km) ||
        error("crossing_altitudes_km must be sorted")
    all(diff(config.crossing_altitudes_km) .> 0) ||
        error("crossing_altitudes_km must be strictly increasing")
    first(config.crossing_altitudes_km) ==
        config.domain_minimum_altitude_km ||
        error("First crossing altitude must equal the lower domain boundary")
    last(config.crossing_altitudes_km) ==
        config.domain_maximum_altitude_km ||
        error("Last crossing altitude must equal the upper domain boundary")
    all(
        (config.domain_minimum_altitude_km .<=
         config.source_altitudes_km) .&
        (config.source_altitudes_km .<=
         config.domain_maximum_altitude_km)
    ) || error("Source altitudes must lie inside the domain")

    rng = Xoshiro(config.seed)
    source_altitudes_m = 1000 .* config.source_altitudes_km
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
    primary_particles =
        length(source_altitudes_m) * config.particles_per_source_altitude
    primary_particles <= config.maximum_total_particles ||
        error("Primary particles exceed maximum_total_particles")
    vibration_probability, vibration_quantum_eV =
        _load_vibration(chemistry_path)

    mkpath(dirname(abspath(output_path)))
    event_count = Ref{Int64}(0)
    result = open(output_path, "w+") do io
        _write_event_header(io)
        queue = _EventTrackedParticle[]
        sizehint!(queue, min(primary_particles, 1_000_000))
        next_id = Int64(0)
        for iz in eachindex(source_altitudes_m)
            altitude_m = source_altitudes_m[iz]
            position_m = (MARS_RADIUS_M + altitude_m, 0.0, 0.0)
            state = source_states[iz]
            for _ in 1:config.particles_per_source_altitude
                next_id += 1
                particle = sample_hot_o_source(
                    rng, position_m, state.Te_K, state.Ti_K, branches;
                    vibrational_probability=vibration_probability,
                    vibrational_quantum_eV=vibration_quantum_eV,
                    weight_s1=source_particle_weights_s1[iz],
                )
                tracked = _EventTrackedParticle(
                    particle, next_id, Int64(0), 0.0, Int32(0),
                )
                push!(queue, tracked)
                _write_event!(io, event_count, tracked, EVENT_BIRTH)
            end
        end

        stops = Dict{Symbol,Int}()
        secondary_count = 0
        next_particle = 1
        maximum_particles_reached = false
        while next_particle <= length(queue)
            tracked = queue[next_particle]
            next_particle += 1
            particle = tracked.particle
            reason = :maximum_steps
            terminal_written = false
            for _ in 1:config.maximum_steps_per_particle
                altitude_km =
                    (_norm(particle.position_m) - MARS_RADIUS_M) / 1000
                energy_eV = kinetic_energy_eV(particle.velocity_m_s)
                if energy_eV <= config.minimum_energy_eV
                    reason = :thermalized
                    tracked.event_index += 1
                    _write_event!(
                        io, event_count, tracked, EVENT_THERMALIZED,
                    )
                    terminal_written = true
                    break
                elseif altitude_km < config.domain_minimum_altitude_km
                    reason = :lower_boundary
                    tracked.event_index += 1
                    _write_event!(
                        io, event_count, tracked, EVENT_EXIT_LOWER,
                    )
                    terminal_written = true
                    break
                elseif altitude_km > config.domain_maximum_altitude_km
                    reason = :upper_boundary
                    tracked.event_index += 1
                    _write_event!(
                        io, event_count, tracked, EVENT_EXIT_UPPER,
                    )
                    terminal_written = true
                    break
                end

                local_state = interpolate_profile(
                    atmosphere, 1000 * altitude_km,
                )
                kappa = collision_coefficient(
                    targets, local_state.density_m3, energy_eV,
                )
                mfp = kappa > 0 ? inv(kappa) : Inf
                ds = rahmati_step_length(mfp)
                position0 = particle.position_m
                velocity0 = particle.velocity_m_s
                position1, velocity1, dt = _advance_gravity(
                    position0, velocity0, ds,
                )
                particle.position_m = position1
                particle.velocity_m_s = velocity1
                tracked.time_s += dt
                _write_crossings!(
                    io, event_count, tracked,
                    position0, velocity0, position1, velocity1,
                    dt, config.crossing_altitudes_km,
                )

                altitude1_km =
                    (_norm(position1) - MARS_RADIUS_M) / 1000
                if altitude1_km < config.domain_minimum_altitude_km
                    reason = :lower_boundary
                    tracked.event_index += 1
                    _write_event!(
                        io, event_count, tracked, EVENT_EXIT_LOWER,
                    )
                    terminal_written = true
                    break
                elseif altitude1_km > config.domain_maximum_altitude_km
                    reason = :upper_boundary
                    tracked.event_index += 1
                    _write_event!(
                        io, event_count, tracked, EVENT_EXIT_UPPER,
                    )
                    terminal_written = true
                    break
                end

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
                       kinetic_energy_eV(target_after) >
                       config.minimum_energy_eV
                        if length(queue) >= config.maximum_total_particles
                            reason = :maximum_particles
                            tracked.event_index += 1
                            _write_event!(
                                io, event_count, tracked,
                                EVENT_MAXIMUM_PARTICLES,
                            )
                            terminal_written = true
                            maximum_particles_reached = true
                            break
                        end
                        next_id += 1
                        secondary = _EventTrackedParticle(
                            HotOParticle(
                                particle.position_m,
                                target_after,
                                particle.weight_s1,
                                true,
                                0,
                            ),
                            next_id,
                            tracked.particle_id,
                            0.0,
                            Int32(0),
                        )
                        push!(queue, secondary)
                        secondary_count += 1
                        _write_event!(
                            io, event_count, secondary, EVENT_BIRTH,
                        )
                    end
                end
            end
            if !terminal_written
                tracked.event_index += 1
                _write_event!(
                    io, event_count, tracked, EVENT_MAXIMUM_STEPS,
                )
            end
            stops[reason] = get(stops, reason, 0) + 1
            maximum_particles_reached && break
        end
        _finalize_event_header(
            io, event_count[], primary_particles, secondary_count,
        )
        HotOCrossingRunResult(
            event_count[],
            primary_particles,
            secondary_count,
            total_source_rate_s1,
            source_particle_weights_s1,
            stops,
            abspath(output_path),
        )
    end
    result
end
