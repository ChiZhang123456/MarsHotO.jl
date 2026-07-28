using MarsHotO

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BATCH_COUNT = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
const PARTICLES_PER_ALTITUDE =
    length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 500
const BASE_SEED = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20260810
const OUTPUT_DIRECTORY = length(ARGS) >= 4 ? abspath(ARGS[4]) : joinpath(
    ROOT, "examples", "output", "run_1p51m_crossings",
)

BATCH_COUNT > 0 || error("BATCH_COUNT must be positive")
PARTICLES_PER_ALTITUDE > 0 ||
    error("PARTICLES_PER_ALTITUDE must be positive")

atmosphere = load_mgitm_subsolar_profile(
    joinpath(ROOT, "MGITM", "MGITM_LS000_F070_150901.dat"),
)
chemistry_path = joinpath(
    ROOT, "data", "chemistry",
    "o2plus_dissociative_recombination.toml",
)
branches = load_reaction_branches(chemistry_path)
targets = load_collision_targets(joinpath(
    ROOT, "data", "cross_sections",
    "rahmati_total_cross_sections.toml",
))

mkpath(OUTPUT_DIRECTORY)
total_primary = 0
total_secondary = 0
total_events = Int64(0)
total_start_time = time()

for batch_index in 1:BATCH_COUNT
    seed = BASE_SEED + batch_index - 1
    output_path = joinpath(
        OUTPUT_DIRECTORY,
        "batch_$(lpad(batch_index, 2, '0')).bin",
    )
    isfile(output_path) && error(
        "Output already exists: $output_path. Choose a new output directory.",
    )
    config = HotOCrossingConfig(
        particles_per_source_altitude=PARTICLES_PER_ALTITUDE,
        seed=seed,
        source_altitudes_km=collect(100.0:1.0:250.0),
        crossing_altitudes_km=collect(100.0:10.0:2000.0),
        domain_minimum_altitude_km=100.0,
        domain_maximum_altitude_km=2000.0,
    )

    println(
        "Batch ", batch_index, "/", BATCH_COUNT,
        ", particles per altitude=", PARTICLES_PER_ALTITUDE,
        ", seed=", seed,
    )
    batch_start_time = time()
    result = run_hot_o_crossing_events(
        atmosphere, targets, branches;
        chemistry_path=chemistry_path,
        output_path=output_path,
        config=config,
    )
    total_primary += result.primary_particles
    total_secondary += result.secondary_particles
    total_events += result.event_records
    println(
        "Completed batch ", batch_index,
        " in ", round(time() - batch_start_time; digits=1), " s",
        ", primary=", result.primary_particles,
        ", secondary=", result.secondary_particles,
        ", events=", result.event_records,
    )
end

println("All batches completed")
println("Elapsed seconds: ", time() - total_start_time)
println("Primary particles: ", total_primary)
println("Secondary particles: ", total_secondary)
println("Tracked particles: ", total_primary + total_secondary)
println("Event records: ", total_events)
println("Output directory: ", OUTPUT_DIRECTORY)
