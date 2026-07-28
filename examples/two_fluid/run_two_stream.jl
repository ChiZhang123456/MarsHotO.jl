using MarsHotO

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

profile = load_mgitm_subsolar_profile(joinpath(
    ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
))
targets = load_collision_targets(joinpath(
    ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
))
branches = load_reaction_branches(joinpath(
    ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
))

config = TwoStreamConfig(
    altitude_min_m=100e3,
    altitude_max_m=300e3,
    altitude_step_m=1e3,
    energy_edges_eV=collect(range(0.01, 7.01; length=71)),
    redistribution_samples=20_000,
    redistribution_seed=73,
    maximum_iterations=2_000,
    relative_tolerance=1e-3,
)
result = run_two_stream(
    profile, targets, branches;
    chemistry_path=joinpath(
        ROOT, "data", "chemistry",
        "o2plus_dissociative_recombination.toml",
    ),
    config=config,
)

output_directory = joinpath(ROOT, "examples", "output")
mkpath(output_directory)
output_path = joinpath(output_directory, "hot_o_two_stream_flux.csv")
write_two_stream_flux(output_path, result)

println("Converged: ", result.converged)
println("Iterations: ", result.iterations)
println("Escape energy at model top: ", result.escape_energy_eV, " eV")
println("Escape flux: ", result.escape_flux_m2_s1, " m^-2 s^-1")
println("Wrote: ", output_path)
