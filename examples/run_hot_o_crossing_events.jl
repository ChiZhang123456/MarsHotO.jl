using MarsHotO

const ROOT = normpath(joinpath(@__DIR__, ".."))
const EVENTS_PER_ALTITUDE =
    length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 250
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260810
const OUTPUT_PATH = length(ARGS) >= 3 ? abspath(ARGS[3]) : joinpath(
    ROOT, "examples", "output", "run_paired_dr_crossings",
    "hot_o_crossing_events.bin",
)

atmosphere = load_mgitm_subsolar_profile(
    joinpath(ROOT, "MGITM", "MGITM_LS000_F070_150901.dat"),
)
chemistry_path = joinpath(
    ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
)
branches = load_reaction_branches(chemistry_path)
targets = load_collision_targets(joinpath(
    ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
))

config = HotOCrossingConfig(
    events_per_source_altitude=EVENTS_PER_ALTITUDE,
    seed=SEED,
    source_altitudes_km=collect(100.0:1.0:250.0),
    crossing_altitudes_km=collect(100.0:10.0:2000.0),
    domain_minimum_altitude_km=100.0,
    domain_maximum_altitude_km=2000.0,
)

println(
    "Running crossing-event model with ",
    EVENTS_PER_ALTITUDE,
    " DR events per source altitude, seed=",
    SEED,
)
start_time = time()
result = run_hot_o_crossing_events(
    atmosphere, targets, branches;
    chemistry_path=chemistry_path,
    output_path=OUTPUT_PATH,
    config=config,
)

println("Elapsed seconds: ", time() - start_time)
println("Event records: ", result.event_records)
println("Primary particles: ", result.primary_particles)
println("Secondary particles: ", result.secondary_particles)
println("Total source rate: ", result.total_source_rate_s1, " s^-1")
println("Stops: ", result.stop_counts)
println("Output: ", result.output_path)
