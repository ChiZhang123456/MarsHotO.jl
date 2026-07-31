using MarsHotO
using Random

const ROOT = normpath(joinpath(@__DIR__, ".."))
const EVENTS_PER_ALTITUDE =
    length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 500
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260730
const OUTPUT_PATH = length(ARGS) >= 3 ? abspath(ARGS[3]) : joinpath(
    ROOT, "examples", "output", "hot_o_time_snapshots.dat",
)

const SNAPSHOT_INTERVAL_S =
    length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : nothing
const MAXIMUM_SNAPSHOT_TIME_S =
    length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 100.0
const SNAPSHOT_TIMES_S = if isnothing(SNAPSHOT_INTERVAL_S)
    [0.0, 10.0, 50.0, 100.0]
else
    SNAPSHOT_INTERVAL_S > 0 ||
        error("SNAPSHOT_INTERVAL_S must be positive")
    collect(0.0:SNAPSHOT_INTERVAL_S:MAXIMUM_SNAPSHOT_TIME_S)
end
const SOURCE_ALTITUDES_KM = collect(100.0:1.0:250.0)
const ALTITUDE_EDGES_KM = collect(100.0:5.0:1000.0)
const ENERGY_EDGES_EV = collect(0.0:0.05:7.0)
const DOMAIN_MINIMUM_ALTITUDE_M = 100e3
const DOMAIN_MAXIMUM_ALTITUDE_M = 2000e3
const MINIMUM_ENERGY_EV = 0.01
const MAXIMUM_TOTAL_PARTICLES = 5_000_000
const MAXIMUM_STEPS_PER_PARTICLE = 2_000_000

mutable struct SnapshotTrackedParticle
    particle::HotOParticle
    time_s::Float64
    next_snapshot::Int
end

function accumulate_snapshot!(
    weighted_histogram, count_histogram, snapshot_index,
    particle_weight_s1, position_m, velocity_m_s,
)
    altitude_km =
        (MarsHotO._norm(position_m) - MARS_RADIUS_M) / 1000
    energy_eV = MarsHotO.kinetic_energy_eV(velocity_m_s)
    altitude_bin = searchsortedlast(ALTITUDE_EDGES_KM, altitude_km)
    energy_bin = searchsortedlast(ENERGY_EDGES_EV, energy_eV)
    if 1 <= altitude_bin < length(ALTITUDE_EDGES_KM) &&
       1 <= energy_bin < length(ENERGY_EDGES_EV)
        weighted_histogram[snapshot_index, altitude_bin, energy_bin] +=
            particle_weight_s1
        count_histogram[snapshot_index, altitude_bin, energy_bin] += 1
    end
end

function interpolate_state(position0, velocity0, position1, velocity1, fraction)
    position = ntuple(
        index -> position0[index] +
                 fraction * (position1[index] - position0[index]),
        3,
    )
    velocity = ntuple(
        index -> velocity0[index] +
                 fraction * (velocity1[index] - velocity0[index]),
        3,
    )
    position, velocity
end

function main()
    EVENTS_PER_ALTITUDE > 0 ||
        error("EVENTS_PER_ALTITUDE must be positive")
    atmosphere = load_mgitm_subsolar_profile(joinpath(
        ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
    ))
    chemistry_path = joinpath(
        ROOT, "data", "chemistry",
        "o2plus_dissociative_recombination.toml",
    )
    branches = load_reaction_branches(chemistry_path)
    targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections",
        "rahmati_total_cross_sections.toml",
    ))
    vibration_probability, vibration_quantum_eV =
        MarsHotO._load_vibration(chemistry_path)

    source_altitudes_m = 1000 .* SOURCE_ALTITUDES_KM
    source_edges_m = MarsHotO._shell_edges_m(source_altitudes_m)
    shell_volume_m3 = (4pi / 3) .* (
        (MARS_RADIUS_M .+ source_edges_m[2:end]).^3 .-
        (MARS_RADIUS_M .+ source_edges_m[1:end-1]).^3
    )
    source_states = [
        interpolate_profile(atmosphere, altitude_m)
        for altitude_m in source_altitudes_m
    ]
    source_event_rate_s1 = [
        dissociative_recombination_event_rate(
            state.density_m3[:e],
            state.density_m3[:O2p],
            state.Te_K,
        ) for state in source_states
    ] .* shell_volume_m3
    source_event_weights_s1 =
        source_event_rate_s1 ./ EVENTS_PER_ALTITUDE

    rng = Xoshiro(SEED)
    queue = SnapshotTrackedParticle[]
    primary_events = length(SOURCE_ALTITUDES_KM) * EVENTS_PER_ALTITUDE
    primary_particles = 2primary_events
    sizehint!(queue, min(primary_particles, MAXIMUM_TOTAL_PARTICLES))
    for source_index in eachindex(source_altitudes_m)
        altitude_m = source_altitudes_m[source_index]
        position_m = (MARS_RADIUS_M + altitude_m, 0.0, 0.0)
        state = source_states[source_index]
        for _ in 1:EVENTS_PER_ALTITUDE
            event = sample_dissociative_recombination_event(
                rng, position_m, state.Te_K, state.Ti_K, branches;
                vibrational_probability=vibration_probability,
                vibrational_quantum_eV=vibration_quantum_eV,
                weight_s1=source_event_weights_s1[source_index],
            )
            for particle in event.products
                push!(queue, SnapshotTrackedParticle(particle, 0.0, 1))
            end
        end
    end

    weighted_histogram = zeros(
        length(SNAPSHOT_TIMES_S),
        length(ALTITUDE_EDGES_KM) - 1,
        length(ENERGY_EDGES_EV) - 1,
    )
    count_histogram = zeros(
        Int64,
        size(weighted_histogram),
    )
    secondary_particles = 0
    stop_counts = Dict{Symbol,Int}()
    next_particle = 1
    maximum_time_s = last(SNAPSHOT_TIMES_S)

    while next_particle <= length(queue)
        tracked = queue[next_particle]
        next_particle += 1
        particle = tracked.particle
        reason = :maximum_steps

        for _ in 1:MAXIMUM_STEPS_PER_PARTICLE
            altitude_m =
                MarsHotO._norm(particle.position_m) - MARS_RADIUS_M
            energy_eV = MarsHotO.kinetic_energy_eV(particle.velocity_m_s)

            while tracked.next_snapshot <= length(SNAPSHOT_TIMES_S) &&
                  SNAPSHOT_TIMES_S[tracked.next_snapshot] == tracked.time_s
                accumulate_snapshot!(
                    weighted_histogram, count_histogram,
                    tracked.next_snapshot, particle.weight_s1,
                    particle.position_m, particle.velocity_m_s,
                )
                tracked.next_snapshot += 1
            end

            if tracked.time_s >= maximum_time_s
                reason = :maximum_snapshot_time
                break
            elseif energy_eV <= MINIMUM_ENERGY_EV
                reason = :thermalized
                break
            elseif altitude_m < DOMAIN_MINIMUM_ALTITUDE_M
                reason = :lower_boundary
                break
            elseif altitude_m > DOMAIN_MAXIMUM_ALTITUDE_M
                reason = :upper_boundary
                break
            end

            position0 = particle.position_m
            velocity0 = particle.velocity_m_s
            time0_s = tracked.time_s
            step = advance_hot_o_step!(
                rng, particle, atmosphere, targets;
                minimum_secondary_energy_eV=MINIMUM_ENERGY_EV,
            )
            position1 = step.position_after_m
            velocity1 = step.ballistic_velocity_after_m_s
            time1_s = time0_s + step.dt_s

            while tracked.next_snapshot <= length(SNAPSHOT_TIMES_S)
                snapshot_time_s = SNAPSHOT_TIMES_S[tracked.next_snapshot]
                snapshot_time_s <= time1_s || break
                fraction = (snapshot_time_s - time0_s) / step.dt_s
                snapshot_position, snapshot_velocity = interpolate_state(
                    position0, velocity0, position1, velocity1, fraction,
                )
                accumulate_snapshot!(
                    weighted_histogram, count_histogram,
                    tracked.next_snapshot, particle.weight_s1,
                    snapshot_position, snapshot_velocity,
                )
                tracked.next_snapshot += 1
            end

            tracked.time_s = time1_s

            if !isnothing(step.secondary)
                    length(queue) < MAXIMUM_TOTAL_PARTICLES ||
                        error("Maximum particle queue size reached")
                    next_snapshot = searchsortedfirst(
                        SNAPSHOT_TIMES_S, tracked.time_s,
                    )
                    push!(
                        queue,
                        SnapshotTrackedParticle(
                            step.secondary, tracked.time_s, next_snapshot,
                        ),
                    )
                    secondary_particles += 1
            end
        end
        stop_counts[reason] = get(stop_counts, reason, 0) + 1
    end

    weighted_probability = zeros(size(weighted_histogram))
    for snapshot_index in axes(weighted_histogram, 1)
        for altitude_index in axes(weighted_histogram, 2)
            total_weight = sum(@view weighted_histogram[
                snapshot_index, altitude_index, :,
            ])
            if total_weight > 0
                weighted_probability[
                    snapshot_index, altitude_index, :,
                ] .= @view(weighted_histogram[
                    snapshot_index, altitude_index, :,
                ]) ./ total_weight
            end
        end
    end

    mkpath(dirname(OUTPUT_PATH))
    open(OUTPUT_PATH, "w") do io
        println(
            io,
            "# time_s altitude_km energy_eV weighted_probability ",
            "particle_count weighted_rate_s-1",
        )
        for snapshot_index in eachindex(SNAPSHOT_TIMES_S)
            for altitude_index in 1:(length(ALTITUDE_EDGES_KM) - 1)
                altitude_km = (
                    ALTITUDE_EDGES_KM[altitude_index] +
                    ALTITUDE_EDGES_KM[altitude_index + 1]
                ) / 2
                for energy_index in 1:(length(ENERGY_EDGES_EV) - 1)
                    energy_eV = (
                        ENERGY_EDGES_EV[energy_index] +
                        ENERGY_EDGES_EV[energy_index + 1]
                    ) / 2
                    println(
                        io,
                        SNAPSHOT_TIMES_S[snapshot_index], ' ',
                        altitude_km, ' ',
                        energy_eV, ' ',
                        weighted_probability[
                            snapshot_index, altitude_index, energy_index,
                        ], ' ',
                        count_histogram[
                            snapshot_index, altitude_index, energy_index,
                        ], ' ',
                        weighted_histogram[
                            snapshot_index, altitude_index, energy_index,
                        ],
                    )
                end
            end
        end
    end

    println("Primary particles: ", primary_particles)
    println("Secondary particles: ", secondary_particles)
    println("Tracked particles: ", length(queue))
    println("Stops: ", stop_counts)
    for snapshot_index in eachindex(SNAPSHOT_TIMES_S)
        println(
            "Snapshot ", SNAPSHOT_TIMES_S[snapshot_index],
            " s, particle count=",
            sum(@view count_histogram[snapshot_index, :, :]),
            ", nonempty altitude bins=",
            count(
                sum(@view(
                    count_histogram[snapshot_index, :, :]
                ); dims=2) .> 0
            ),
        )
    end
    println("Output: ", OUTPUT_PATH)
end

main()
