using MarsHotO
using Random

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTPUT_DIR = joinpath(ROOT, "examples", "output")
const TRAJECTORY_OUTPUT = joinpath(
    OUTPUT_DIR, "two_opposite_hot_o_trajectories.dat",
)
const COLLISION_OUTPUT = joinpath(
    OUTPUT_DIR, "two_opposite_hot_o_collisions.dat",
)

const INITIAL_ALTITUDE_KM = 140.0
const RANDOM_SEED = 20260729
const MINIMUM_ENERGY_EV = 0.01
const MAXIMUM_ALTITUDE_KM = 500.0
const MAXIMUM_STEPS = 200_000

const SPECIES_CODE = Dict(
    :O => 1,
    :CO => 2,
    :N2 => 3,
    :O2 => 4,
    :CO2 => 5,
)

function trace_particle(
    particle_id, initial_velocity_m_s, rng, atmosphere, targets,
)
    position_m = (
        MarsHotO.MARS_RADIUS_M + 1000 * INITIAL_ALTITUDE_KM,
        0.0,
        0.0,
    )
    particle = HotOParticle(
        position_m, initial_velocity_m_s, 1.0, true, 0,
    )
    trajectory = NamedTuple[]
    collisions = NamedTuple[]
    elapsed_time_s = 0.0
    path_length_km = 0.0
    collision_index = 0
    reason = :maximum_steps

    function record_trajectory!(step)
        altitude_km = (
            MarsHotO._norm(particle.position_m) - MarsHotO.MARS_RADIUS_M
        ) / 1000
        push!(trajectory, (
            particle_id=particle_id,
            step=step,
            time_s=elapsed_time_s,
            path_km=path_length_km,
            east_km=particle.position_m[2] / 1000,
            north_km=particle.position_m[3] / 1000,
            altitude_km=altitude_km,
            energy_eV=MarsHotO.kinetic_energy_eV(particle.velocity_m_s),
            collision_index=collision_index,
        ))
    end

    record_trajectory!(0)
    for step in 1:MAXIMUM_STEPS
        altitude_m =
            MarsHotO._norm(particle.position_m) - MarsHotO.MARS_RADIUS_M
        energy_before_step_eV =
            MarsHotO.kinetic_energy_eV(particle.velocity_m_s)
        if energy_before_step_eV <= MINIMUM_ENERGY_EV
            reason = :thermalized
            break
        elseif altitude_m >= 1000 * MAXIMUM_ALTITUDE_KM
            reason = :upper_boundary
            break
        elseif altitude_m <= atmosphere.altitude_m[1]
            reason = :lower_boundary
            break
        end

        energy_before_collision_eV = energy_before_step_eV
        step_result = advance_hot_o_step!(
            rng, particle, atmosphere, targets;
            minimum_secondary_energy_eV=MINIMUM_ENERGY_EV,
        )
        elapsed_time_s += step_result.dt_s
        path_length_km += step_result.ds_m / 1000

        if !isnothing(step_result.target)
            target = step_result.target
            theta_com_rad = step_result.scattering_angle_com_rad
            collision_index += 1
            energy_after_collision_eV =
                MarsHotO.kinetic_energy_eV(particle.velocity_m_s)
            altitude_collision_km = (
                MarsHotO._norm(particle.position_m) - MarsHotO.MARS_RADIUS_M
            ) / 1000
            push!(collisions, (
                particle_id=particle_id,
                collision_index=collision_index,
                time_s=elapsed_time_s,
                path_km=path_length_km,
                east_km=particle.position_m[2] / 1000,
                north_km=particle.position_m[3] / 1000,
                altitude_km=altitude_collision_km,
                species_code=SPECIES_CODE[target.species],
                theta_com_deg=rad2deg(theta_com_rad),
                energy_before_eV=energy_before_collision_eV,
                energy_after_eV=energy_after_collision_eV,
                fractional_energy_loss=max(
                    1 - energy_after_collision_eV /
                        energy_before_collision_eV,
                    0.0,
                ),
            ))
        end
        record_trajectory!(step)
    end
    trajectory, collisions, reason
end

function write_trajectory(path, rows)
    open(path, "w") do io
        println(
            io,
            "# particle_id step time_s path_km east_km north_km ",
            "altitude_km energy_eV collision_index",
        )
        for row in rows
            println(
                io,
                row.particle_id, ' ', row.step, ' ', row.time_s, ' ',
                row.path_km, ' ', row.east_km, ' ', row.north_km, ' ',
                row.altitude_km, ' ', row.energy_eV, ' ',
                row.collision_index,
            )
        end
    end
end

function write_collisions(path, rows)
    open(path, "w") do io
        println(io, "# species_code: 1=O 2=CO 3=N2 4=O2 5=CO2")
        println(
            io,
            "# particle_id collision_index time_s path_km east_km ",
            "north_km altitude_km species_code theta_com_deg ",
            "energy_before_eV energy_after_eV fractional_energy_loss",
        )
        for row in rows
            println(
                io,
                row.particle_id, ' ', row.collision_index, ' ',
                row.time_s, ' ', row.path_km, ' ', row.east_km, ' ',
                row.north_km, ' ', row.altitude_km, ' ',
                row.species_code, ' ', row.theta_com_deg, ' ',
                row.energy_before_eV, ' ', row.energy_after_eV, ' ',
                row.fractional_energy_loss,
            )
        end
    end
end

function main()
    atmosphere = load_mgitm_subsolar_profile(
        joinpath(ROOT, "MGITM", "MGITM_LS000_F070_150901.dat"),
    )
    targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
    ))
    chemistry_path = joinpath(
        ROOT, "data", "chemistry",
        "o2plus_dissociative_recombination.toml",
    )
    branches = load_reaction_branches(chemistry_path)
    vibration_probability, vibration_quantum_eV =
        MarsHotO._load_vibration(chemistry_path)
    source_state = interpolate_profile(
        atmosphere, 1000INITIAL_ALTITUDE_KM,
    )
    source_position_m = (
        MARS_RADIUS_M + 1000INITIAL_ALTITUDE_KM, 0.0, 0.0,
    )
    source_event = sample_dissociative_recombination_event(
        Xoshiro(RANDOM_SEED), source_position_m,
        source_state.Te_K, source_state.Ti_K, branches;
        vibrational_probability=vibration_probability,
        vibrational_quantum_eV=vibration_quantum_eV,
        plasma_bulk_velocity_m_s=(0.0, 0.0, 0.0),
    )
    velocity_1 = source_event.products[1].velocity_m_s
    velocity_2 = source_event.products[2].velocity_m_s

    all_trajectories = NamedTuple[]
    all_collisions = NamedTuple[]
    reasons = Symbol[]
    for (particle_id, velocity) in enumerate((velocity_1, velocity_2))
        trajectory, collisions, reason = trace_particle(
            particle_id,
            velocity,
            Xoshiro(RANDOM_SEED + particle_id),
            atmosphere,
            targets,
        )
        append!(all_trajectories, trajectory)
        append!(all_collisions, collisions)
        push!(reasons, reason)
    end

    mkpath(OUTPUT_DIR)
    write_trajectory(TRAJECTORY_OUTPUT, all_trajectories)
    write_collisions(COLLISION_OUTPUT, all_collisions)

    initial_momentum = MarsHotO._add(
        MarsHotO._scale(O_MASS_KG, velocity_1),
        MarsHotO._scale(O_MASS_KG, velocity_2),
    )
    println("Initial altitude: ", INITIAL_ALTITUDE_KM, " km")
    println("Reaction branch: ", source_event.branch.products)
    println("Vibrational level: ", source_event.vibrational_level)
    println("Available energy: ", source_event.available_energy_eV, " eV")
    println("Event COM velocity: ", source_event.com_velocity_m_s, " m s^-1")
    println("Event COM speed: ", MarsHotO._norm(
        source_event.com_velocity_m_s,
    ), " m s^-1")
    println("Initial O energies: ",
        MarsHotO.kinetic_energy_eV(velocity_1), ", ",
        MarsHotO.kinetic_energy_eV(velocity_2), " eV")
    println(
        "Initial total momentum magnitude: ",
        MarsHotO._norm(initial_momentum),
        " kg m s^-1",
    )
    for particle_id in 1:2
        selected = filter(
            row -> row.particle_id == particle_id,
            all_collisions,
        )
        final_row = last(filter(
            row -> row.particle_id == particle_id,
            all_trajectories,
        ))
        println(
            "O", particle_id,
            ": stop=", reasons[particle_id],
            ", collisions=", length(selected),
            ", final altitude=", final_row.altitude_km, " km",
            ", final energy=", final_row.energy_eV, " eV",
        )
        for event in selected
            species = first(
                key for (key, value) in SPECIES_CODE
                if value == event.species_code
            )
            println(
                "  collision ", event.collision_index,
                ": target=", species,
                ", altitude=", event.altitude_km, " km",
                ", theta_COM=", event.theta_com_deg, " deg",
                ", E=", event.energy_before_eV,
                " -> ", event.energy_after_eV, " eV",
            )
        end
    end
    println("Trajectory output: ", TRAJECTORY_OUTPUT)
    println("Collision output: ", COLLISION_OUTPUT)
end

main()
