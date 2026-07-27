using MarsHotO

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PARTICLES = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 10_000
const SEED = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260727

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

config = RahmatiMonteCarloConfig(
    primary_particles=PARTICLES,
    seed=SEED,
    altitude_edges_km=collect(100.0:10.0:2000.0),
    energy_edges_eV=collect(range(0.01, 7.0; length=141)),
)

println("Running MarsHotO with $PARTICLES primary particles, seed=$SEED")
result = run_hot_o_corona(
    atmosphere, targets, branches;
    chemistry_path=chemistry_path,
    config=config,
)

output_dir = joinpath(@__DIR__, "output")
mkpath(output_dir)
output_path = joinpath(output_dir, "hot_o_altitude_energy_distribution.dat")
write_corona_distribution(output_path, result)

println("Primary particles: ", result.primary_particles)
println("Secondary particles: ", result.secondary_particles)
println("Total source rate: ", result.total_source_rate_s1, " s^-1")
println("Stops: ", result.stop_counts)
println("Output: ", output_path)
