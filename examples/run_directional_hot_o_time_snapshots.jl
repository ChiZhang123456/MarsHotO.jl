using MarsHotO
using Random

const ROOT = normpath(joinpath(@__DIR__, ".."))
const EVENTS_PER_ALTITUDE = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260731
const SNAPSHOT_TIMES_S = collect(0.0:5.0:100.0)
const SOURCE_ALTITUDES_KM = collect(100.0:1.0:250.0)
const ALTITUDE_EDGES_KM = collect(100.0:5.0:400.0)
const ENERGY_EDGES_EV = collect(0.0:0.1:7.0)
const MINIMUM_ENERGY_EV = 0.01
const MAXIMUM_STEPS = 2_000_000
const MAXIMUM_PARTICLES_PER_ALTITUDE = 200_000
const OUTPUT = joinpath(
    ROOT, "examples", "output", "directional_hot_o_time_snapshots.dat",
)

mutable struct TimedParticle
    particle::HotOParticle
    time_s::Float64
    next_snapshot::Int
end

function accumulate!(histogram, snapshot_index, particle, position, velocity)
    altitude_km = (MarsHotO._norm(position) - MARS_RADIUS_M) / 1000
    energy_eV = MarsHotO.kinetic_energy_eV(velocity)
    ia = searchsortedlast(ALTITUDE_EDGES_KM, altitude_km)
    ie = searchsortedlast(ENERGY_EDGES_EV, energy_eV)
    if 1 <= ia < length(ALTITUDE_EDGES_KM) &&
       1 <= ie < length(ENERGY_EDGES_EV)
        radial_velocity = MarsHotO._dot(position, velocity) /
                          MarsHotO._norm(position)
        direction_index = radial_velocity >= 0 ? 1 : 2
        histogram[snapshot_index, ia, ie, direction_index] +=
            particle.weight_s1
    end
end

function interpolate_state(position0, velocity0, position1, velocity1, fraction)
    position = ntuple(
        i -> position0[i] + fraction * (position1[i] - position0[i]), 3,
    )
    velocity = ntuple(
        i -> velocity0[i] + fraction * (velocity1[i] - velocity0[i]), 3,
    )
    position, velocity
end

function simulate_altitude!(
    histogram, rng, atmosphere, targets, branches, source_altitude_km,
    source_state, event_weight_s1, vibration_probability,
    vibration_quantum_eV,
)
    position = (MARS_RADIUS_M + 1000source_altitude_km, 0.0, 0.0)
    queue = TimedParticle[]
    sizehint!(queue, 2EVENTS_PER_ALTITUDE)
    for _ in 1:EVENTS_PER_ALTITUDE
        event = sample_dissociative_recombination_event(
            rng, position, source_state.Te_K, source_state.Ti_K, branches;
            vibrational_probability=vibration_probability,
            vibrational_quantum_eV=vibration_quantum_eV,
            weight_s1=event_weight_s1,
            plasma_bulk_velocity_m_s=(0.0, 0.0, 0.0),
        )
        for particle in event.products
            push!(queue, TimedParticle(particle, 0.0, 1))
        end
    end

    next_particle = 1
    secondaries = 0
    while next_particle <= length(queue)
        tracked = queue[next_particle]
        next_particle += 1
        particle = tracked.particle
        for _ in 1:MAXIMUM_STEPS
            while tracked.next_snapshot <= length(SNAPSHOT_TIMES_S) &&
                  SNAPSHOT_TIMES_S[tracked.next_snapshot] == tracked.time_s
                accumulate!(
                    histogram, tracked.next_snapshot, particle,
                    particle.position_m, particle.velocity_m_s,
                )
                tracked.next_snapshot += 1
            end
            tracked.next_snapshot > length(SNAPSHOT_TIMES_S) && break
            altitude_m = MarsHotO._norm(particle.position_m) - MARS_RADIUS_M
            energy_eV = MarsHotO.kinetic_energy_eV(particle.velocity_m_s)
            if energy_eV <= MINIMUM_ENERGY_EV || altitude_m <= 100e3 ||
               altitude_m >= 400e3
                break
            end

            position0 = particle.position_m
            velocity0 = particle.velocity_m_s
            time0 = tracked.time_s
            step = advance_hot_o_step!(
                rng, particle, atmosphere, targets;
                minimum_secondary_energy_eV=MINIMUM_ENERGY_EV,
            )
            time1 = time0 + step.dt_s
            while tracked.next_snapshot <= length(SNAPSHOT_TIMES_S) &&
                  SNAPSHOT_TIMES_S[tracked.next_snapshot] <= time1
                snapshot_time = SNAPSHOT_TIMES_S[tracked.next_snapshot]
                fraction = (snapshot_time - time0) / step.dt_s
                snapshot_position, snapshot_velocity = interpolate_state(
                    position0, velocity0, step.position_after_m,
                    step.ballistic_velocity_after_m_s, fraction,
                )
                accumulate!(
                    histogram, tracked.next_snapshot, particle,
                    snapshot_position, snapshot_velocity,
                )
                tracked.next_snapshot += 1
            end
            tracked.time_s = time1
            if !isnothing(step.secondary)
                length(queue) < MAXIMUM_PARTICLES_PER_ALTITUDE || break
                next_snapshot = searchsortedfirst(SNAPSHOT_TIMES_S, time1)
                push!(queue, TimedParticle(
                    step.secondary, time1, next_snapshot,
                ))
                secondaries += 1
            end
        end
    end
    secondaries
end

function main()
    atmosphere = load_mgitm_subsolar_profile(joinpath(
        ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
    ))
    chemistry_path = joinpath(
        ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
    )
    branches = load_reaction_branches(chemistry_path)
    targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
    ))
    vibration_probability, vibration_quantum_eV =
        MarsHotO._load_vibration(chemistry_path)
    source_altitudes_m = 1000 .* SOURCE_ALTITUDES_KM
    source_edges_m = MarsHotO._shell_edges_m(source_altitudes_m)
    shell_volumes = (4pi / 3) .* (
        (MARS_RADIUS_M .+ source_edges_m[2:end]).^3 .-
        (MARS_RADIUS_M .+ source_edges_m[1:end-1]).^3
    )
    source_states = [
        interpolate_profile(atmosphere, altitude_m)
        for altitude_m in source_altitudes_m
    ]
    event_weights = [
        dissociative_recombination_event_rate(
            state.density_m3[:e], state.density_m3[:O2p], state.Te_K,
        ) * shell_volumes[i] / EVENTS_PER_ALTITUDE
        for (i, state) in enumerate(source_states)
    ]

    thread_count = Threads.maxthreadid()
    histograms = zeros(
        Float64, length(SNAPSHOT_TIMES_S),
        length(ALTITUDE_EDGES_KM) - 1, length(ENERGY_EDGES_EV) - 1,
        2, thread_count,
    )
    completed = Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()
    Threads.@threads :dynamic for i in eachindex(SOURCE_ALTITUDES_KM)
        thread_id = Threads.threadid()
        secondaries = simulate_altitude!(
            @view(histograms[:, :, :, :, thread_id]),
            Xoshiro(SEED + i - 1), atmosphere, targets, branches,
            SOURCE_ALTITUDES_KM[i], source_states[i], event_weights[i],
            vibration_probability, vibration_quantum_eV,
        )
        count = Threads.atomic_add!(completed, 1) + 1
        lock(progress_lock) do
            println(
                "Completed ", SOURCE_ALTITUDES_KM[i], " km, ", count, "/",
                length(SOURCE_ALTITUDES_KM), ", thread=", thread_id,
                ", secondaries=", secondaries,
            )
            flush(stdout)
        end
    end
    weighted_rate = dropdims(sum(histograms; dims=5); dims=5)

    altitude_radii = MARS_RADIUS_M .+ 1000 .* ALTITUDE_EDGES_KM
    volumes = (4pi / 3) .* (
        altitude_radii[2:end].^3 .- altitude_radii[1:end-1].^3
    )
    energy_widths = diff(ENERGY_EDGES_EV)
    density_cm3_eV1 = weighted_rate ./
        reshape(volumes, 1, :, 1, 1) ./
        reshape(energy_widths, 1, 1, :, 1) ./ 1e6

    mkpath(dirname(OUTPUT))
    open(OUTPUT, "w") do io
        println(io, "# time_s altitude_km energy_eV upward_cm-3_eV-1 downward_cm-3_eV-1")
        for it in eachindex(SNAPSHOT_TIMES_S),
            ia in 1:(length(ALTITUDE_EDGES_KM) - 1),
            ie in 1:(length(ENERGY_EDGES_EV) - 1)
            altitude = (ALTITUDE_EDGES_KM[ia] + ALTITUDE_EDGES_KM[ia + 1]) / 2
            energy = (ENERGY_EDGES_EV[ie] + ENERGY_EDGES_EV[ie + 1]) / 2
            println(
                io, SNAPSHOT_TIMES_S[it], ' ', altitude, ' ', energy, ' ',
                density_cm3_eV1[it, ia, ie, 1], ' ',
                density_cm3_eV1[it, ia, ie, 2],
            )
        end
    end
    println("Output: ", OUTPUT)
end

main()
