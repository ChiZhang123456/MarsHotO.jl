using MarsHotO

const ROOT = normpath(joinpath(@__DIR__, ".."))
const BATCH_COUNT = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
const EVENTS_PER_ALTITUDE =
    length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 250
const BASE_SEED = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 20260810
const OUTPUT_DIRECTORY = length(ARGS) >= 4 ? abspath(ARGS[4]) : joinpath(
    ROOT, "examples", "output", "run_paired_dr_crossings",
)

BATCH_COUNT > 0 || error("BATCH_COUNT must be positive")
EVENTS_PER_ALTITUDE > 0 ||
    error("EVENTS_PER_ALTITUDE must be positive")

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
primary_by_batch = zeros(Int, BATCH_COUNT)
secondary_by_batch = zeros(Int, BATCH_COUNT)
events_by_batch = zeros(Int64, BATCH_COUNT)
total_start_time = time()
progress_lock = ReentrantLock()

Threads.@threads :dynamic for batch_index in 1:BATCH_COUNT
    seed = BASE_SEED + batch_index - 1
    output_path = joinpath(
        OUTPUT_DIRECTORY,
        "batch_$(lpad(batch_index, 2, '0')).bin",
    )
    isfile(output_path) && error(
        "Output already exists: $output_path. Choose a new output directory.",
    )
    config = HotOCrossingConfig(
        events_per_source_altitude=EVENTS_PER_ALTITUDE,
        seed=seed,
        source_altitudes_km=collect(100.0:1.0:250.0),
        crossing_altitudes_km=collect(100.0:10.0:2000.0),
        domain_minimum_altitude_km=100.0,
        domain_maximum_altitude_km=2000.0,
    )

    lock(progress_lock) do
        println(
            "Batch ", batch_index, "/", BATCH_COUNT,
            ", DR events per altitude=", EVENTS_PER_ALTITUDE,
            ", seed=", seed, ", thread=", Threads.threadid(),
        )
        flush(stdout)
    end
    batch_start_time = time()
    result = run_hot_o_crossing_events(
        atmosphere, targets, branches;
        chemistry_path=chemistry_path,
        output_path=output_path,
        config=config,
    )
    primary_by_batch[batch_index] = result.primary_particles
    secondary_by_batch[batch_index] = result.secondary_particles
    events_by_batch[batch_index] = result.event_records
    lock(progress_lock) do
        println(
            "Completed batch ", batch_index,
            " in ", round(time() - batch_start_time; digits=1), " s",
            ", primary=", result.primary_particles,
            ", secondary=", result.secondary_particles,
            ", events=", result.event_records,
        )
        flush(stdout)
    end
end

println("All batches completed")
println("Elapsed seconds: ", time() - total_start_time)
println("Primary particles: ", sum(primary_by_batch))
println("Secondary particles: ", sum(secondary_by_batch))
println("Tracked particles: ", sum(primary_by_batch) + sum(secondary_by_batch))
println("Event records: ", sum(events_by_batch))
println("Output directory: ", OUTPUT_DIRECTORY)
