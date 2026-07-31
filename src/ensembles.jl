"""
Weighted particle ensembles and residence-time diagnostics.

This file builds a stratified source distribution, tracks primary and
secondary particles, records termination reasons, and converts residence time
into the spherical density in each altitude and energy bin, in m^-3.
"""
Base.@kwdef struct RahmatiMonteCarloConfig
    events_per_source_altitude::Int = 5_000
    seed::Int = 20260727
    minimum_energy_eV::Float64 = 0.01
    maximum_altitude_m::Float64 = 110_000e3
    maximum_steps_per_particle::Int = 2_000_000
    maximum_total_particles::Int = 50_000_000
    show_progress::Bool = true
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
    source_event_weights_s1::Vector{Float64}
    events_per_source_altitude::Int
    primary_events::Int
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

function _simulate_source_altitude!(
    upward_residence_particles, downward_residence_particles,
    rng, atmosphere, targets, branches, position, source_state,
    event_count, event_weight_s1, vibration_probability,
    vibration_quantum_eV, config, global_particle_count,
)
    queue = HotOParticle[]
    sizehint!(queue, 2event_count)
    for _ in 1:event_count
        event = sample_dissociative_recombination_event(
            rng, position, source_state.Te_K, source_state.Ti_K, branches;
            vibrational_probability=vibration_probability,
            vibrational_quantum_eV=vibration_quantum_eV,
            weight_s1=event_weight_s1,
            plasma_bulk_velocity_m_s=(0.0, 0.0, 0.0),
        )
        push!(queue, event.products...)
    end

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

            step = advance_hot_o_step!(
                rng, particle, atmosphere, targets;
                minimum_secondary_energy_eV=config.minimum_energy_eV,
            )
            radial_velocity_m_s =
                _dot(step.position_before_m, step.velocity_before_m_s) /
                _norm(step.position_before_m)
            directional_residence_particles = radial_velocity_m_s >= 0 ?
                upward_residence_particles : downward_residence_particles
            _accumulate_residence!(
                directional_residence_particles, particle.weight_s1,
                step.dt_s, altitude_m, energy_eV,
                config.altitude_edges_km, config.energy_edges_eV,
            )
            if !isnothing(step.secondary)
                previous_count = Threads.atomic_add!(global_particle_count, 1)
                if previous_count >= config.maximum_total_particles
                    Threads.atomic_sub!(global_particle_count, 1)
                    reason = :maximum_particles
                    break
                end
                push!(queue, step.secondary)
                secondary_count += 1
            end
        end
        stops[reason] = get(stops, reason, 0) + 1
        reason == :maximum_particles && break
    end
    secondary_count, stops
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
    config.events_per_source_altitude > 0 ||
        error("events_per_source_altitude must be positive")
    length(config.source_altitudes_km) >= 2 ||
        error("At least two source altitudes are required")
    issorted(config.source_altitudes_km) ||
        error("source_altitudes_km must be sorted")
    all(diff(config.source_altitudes_km) .> 0) ||
        error("source_altitudes_km must be strictly increasing")
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
    source_event_rate_m3_s1 = [
        dissociative_recombination_event_rate(
            state.density_m3[:e],
            state.density_m3[:O2p],
            state.Te_K,
        ) for state in source_states
    ]
    source_event_rate_s1 = source_event_rate_m3_s1 .* shell_volume_m3
    source_event_weights_s1 =
        source_event_rate_s1 ./ config.events_per_source_altitude
    total_source_rate_s1 = 2sum(source_event_rate_s1)
    total_source_rate_s1 > 0 || error("The atmosphere has zero hot O source")
    primary_events =
        length(source_altitudes_m) * config.events_per_source_altitude
    primary_particles = 2primary_events
    primary_particles <= config.maximum_total_particles ||
        error("Primary particles exceed maximum_total_particles")
    vibration_probability, vibration_quantum_eV =
        _load_vibration(chemistry_path)

    altitude_bin_count = length(config.altitude_edges_km) - 1
    energy_bin_count = length(config.energy_edges_eV) - 1
    # Julia may assign work from a nondefault thread pool whose thread ID is
    # larger than Threads.nthreads(:default). Allocate by the largest possible
    # thread ID so every worker has a private reduction slot.
    thread_count = Threads.maxthreadid()
    upward_by_thread = zeros(
        altitude_bin_count, energy_bin_count, thread_count,
    )
    downward_by_thread = zeros(
        altitude_bin_count, energy_bin_count, thread_count,
    )
    stops_by_thread = [Dict{Symbol,Int}() for _ in 1:thread_count]
    secondary_by_thread = zeros(Int, thread_count)
    global_particle_count = Threads.Atomic{Int}(primary_particles)
    completed_altitudes = Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()

    Threads.@threads :dynamic for iz in eachindex(source_altitudes_m)
        thread_id = Threads.threadid()
        rng = Xoshiro(config.seed + iz - 1)
        altitude_m = source_altitudes_m[iz]
        position = (MARS_RADIUS_M + altitude_m, 0.0, 0.0)
        secondary_count, local_stops = _simulate_source_altitude!(
            @view(upward_by_thread[:, :, thread_id]),
            @view(downward_by_thread[:, :, thread_id]),
            rng, atmosphere, targets, branches, position, source_states[iz],
            config.events_per_source_altitude,
            source_event_weights_s1[iz], vibration_probability,
            vibration_quantum_eV, config, global_particle_count,
        )
        secondary_by_thread[thread_id] += secondary_count
        for (reason, count) in local_stops
            stops_by_thread[thread_id][reason] =
                get(stops_by_thread[thread_id], reason, 0) + count
        end
        completed = Threads.atomic_add!(completed_altitudes, 1) + 1
        if config.show_progress
            lock(progress_lock) do
                println(
                    "Completed source altitude ",
                    config.source_altitudes_km[iz], " km, ", completed, "/",
                    length(source_altitudes_m), ", thread=", thread_id,
                    ", secondaries=", secondary_count,
                )
                flush(stdout)
            end
        end
    end

    upward_residence_particles = dropdims(
        sum(upward_by_thread; dims=3); dims=3,
    )
    downward_residence_particles = dropdims(
        sum(downward_by_thread; dims=3); dims=3,
    )
    stops = Dict{Symbol,Int}()
    for local_stops in stops_by_thread, (reason, count) in local_stops
        stops[reason] = get(stops, reason, 0) + count
    end
    secondary_count = sum(secondary_by_thread)

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
        source_event_weights_s1,
        config.events_per_source_altitude,
        primary_events, primary_particles, secondary_count,
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

"""Write the complete compact output set for one steady-state corona run."""
function write_corona_outputs(
    output_directory::AbstractString, result::HotOCoronaResult,
)
    mkpath(output_directory)
    total_path = joinpath(output_directory, "hot_o_density_total.dat")
    directional_path =
        joinpath(output_directory, "hot_o_density_directional.dat")
    source_path = joinpath(output_directory, "hot_o_source_events.dat")
    summary_path = joinpath(output_directory, "hot_o_run_summary.toml")
    write_corona_distribution(total_path, result)
    write_directional_corona_distribution(directional_path, result)
    open(source_path, "w") do io
        println(
            io,
            "# source_altitude_km event_weight_s-1 events primary_O",
        )
        for (altitude_km, weight_s1) in zip(
            result.source_altitudes_km, result.source_event_weights_s1,
        )
            println(
                io, altitude_km, ' ', weight_s1, ' ',
                result.events_per_source_altitude, ' ',
                2result.events_per_source_altitude,
            )
        end
    end
    summary = Dict{String,Any}(
        "events_per_source_altitude" => result.events_per_source_altitude,
        "source_altitude_count" => length(result.source_altitudes_km),
        "primary_events" => result.primary_events,
        "primary_particles" => result.primary_particles,
        "secondary_particles" => result.secondary_particles,
        "tracked_particles" =>
            result.primary_particles + result.secondary_particles,
        "total_hot_o_source_rate_s1" => result.total_source_rate_s1,
        "density_unit" => "m^-3 per energy bin",
        "energy_bin_unit" => "eV",
        "altitude_unit" => "km",
        "spherical_atmosphere_approximation" => true,
        "stop_counts" => Dict(
            String(reason) => count for (reason, count) in result.stop_counts
        ),
        "files" => Dict(
            "total_density" => basename(total_path),
            "directional_density" => basename(directional_path),
            "source_events" => basename(source_path),
        ),
    )
    open(summary_path, "w") do io
        TOML.print(io, summary; sorted=true)
    end
    (
        total_density=total_path,
        directional_density=directional_path,
        source_events=source_path,
        summary=summary_path,
    )
end
