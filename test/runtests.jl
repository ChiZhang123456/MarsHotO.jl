using MarsHotO
using Random
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "Dissociative recombination" begin
    @test dissociative_recombination_coefficient(300.0) ≈ 1.95e-13
    @test dissociative_recombination_coefficient(1200.0) ≈
          1.95e-13 * (300 / 1200)^0.70
    @test dissociative_recombination_coefficient(1200.001) ≈
          7.39e-14 * (1200 / 1200.001)^0.56
    branches = load_reaction_branches(joinpath(
        ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
    ))
    @test getfield.(branches, :probability) == [0.265, 0.473, 0.204, 0.058]
end

@testset "Cross sections" begin
    targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
    ))
    oxygen = only(filter(x -> x.species == :O, targets))
    @test total_cross_section(oxygen, 3.0) ≈ 6.4e-19
    @test total_cross_section(oxygen, 12.0) <
          total_cross_section(oxygen, 3.0)
    density_m3 = Dict(:O => 1.0e15)
    @test collision_coefficient(targets, density_m3, 3.0) ≈
          1.0e15 * 6.4e-19
end

@testset "Two-stream transport" begin
    profile = load_mgitm_subsolar_profile(joinpath(
        ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
    ))
    branches = load_reaction_branches(joinpath(
        ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
    ))
    all_targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
    ))
    targets = filter(target -> target.species in (:O, :CO2), all_targets)
    edges = collect(0.01:0.5:7.01)
    redistribution1 = build_two_stream_redistribution(
        targets, edges; samples=300, seed=73,
    )
    redistribution2 = build_two_stream_redistribution(
        targets, edges; samples=300, seed=73,
    )
    @test redistribution1.same_stream == redistribution2.same_stream
    @test all(redistribution1.same_stream .>= 0)
    @test all(redistribution1.reverse_stream .>= 0)
    for target_index in eachindex(targets), source_index in 1:length(edges)-1
        projectile_probability = sum(
            redistribution1.same_stream[target_index, source_index, :] +
            redistribution1.reverse_stream[target_index, source_index, :],
        )
        @test projectile_probability <= 1 + 1e-12
    end
    oxygen_index = findfirst(target -> target.species == :O, targets)
    @test sum(redistribution1.secondary_same_stream[oxygen_index, :, :]) > 0

    config = TwoStreamConfig(
        altitude_min_m=100e3,
        altitude_max_m=110e3,
        altitude_step_m=2e3,
        energy_edges_eV=edges,
        redistribution_samples=300,
        maximum_iterations=100,
        relative_tolerance=1e-5,
    )
    result = run_two_stream(
        profile, targets, branches;
        chemistry_path=joinpath(
            ROOT, "data", "chemistry",
            "o2plus_dissociative_recombination.toml",
        ),
        config=config,
        redistribution=redistribution1,
    )
    @test size(result.upward_flux_m2_s1) ==
          (length(result.altitude_m), length(edges) - 1)
    @test all(isfinite, result.upward_flux_m2_s1)
    @test all(result.upward_flux_m2_s1 .>= 0)
    @test all(result.downward_flux_m2_s1 .>= 0)
    @test result.upward_flux_m2_s1[1, :] ≈
          result.downward_flux_m2_s1[1, :]
    @test result.escape_flux_m2_s1 >= 0
end

@testset "Scattering distribution" begin
    beta = -1.85
    exponent = beta + 2
    @test scattering_angle_cdf(0.0) == 0.0
    @test scattering_angle_cdf(pi) == 1.0
    median_angle = 2asin(0.5^(1 / exponent))
    @test scattering_angle_cdf(median_angle) ≈ 0.5
    for probability in range(0.0, 1.0; length=101)
        theta = 2asin(probability^(1 / exponent))
        @test scattering_angle_cdf(theta) ≈ probability atol=1e-13
    end
    rng = Xoshiro(73)
    samples = [sample_scattering_angle(rng) for _ in 1:100_000]
    @test all(0.0 .<= samples .<= pi)
    empirical = sum(x <= deg2rad(10.0) for x in samples) / length(samples)
    expected = scattering_angle_cdf(deg2rad(10.0))
    @test abs(empirical - expected) < 0.005
    @test minimum(samples) < deg2rad(0.01)
end

@testset "Elastic collision" begin
    v1 = (6000.0, 0.0, 0.0)
    v2 = (0.0, 0.0, 0.0)
    m1, m2 = 16.0, 44.0
    before_energy = 0.5m1 * sum(abs2, v1) + 0.5m2 * sum(abs2, v2)
    before_momentum = ntuple(i -> m1 * v1[i] + m2 * v2[i], 3)
    v1_after, v2_after = elastic_collision(v1, v2, m1, m2, 1.1, 2.2)
    after_energy = 0.5m1 * sum(abs2, v1_after) +
                   0.5m2 * sum(abs2, v2_after)
    after_momentum = ntuple(i -> m1 * v1_after[i] + m2 * v2_after[i], 3)
    @test after_energy ≈ before_energy rtol=1e-13
    @test collect(after_momentum) ≈ collect(before_momentum)
    measured_loss = 1 - sum(abs2, v1_after) / sum(abs2, v1)
    @test measured_loss ≈ fractional_energy_loss(1.1, m1, m2)

    theta_lab = deg2rad(42.0)
    v1_after_lab, v2_after_lab = elastic_collision_lab(
        v1, m1, m2, theta_lab, 2.2,
    )
    after_energy_lab = 0.5m1 * sum(abs2, v1_after_lab) +
                       0.5m2 * sum(abs2, v2_after_lab)
    after_momentum_lab =
        ntuple(i -> m1 * v1_after_lab[i] + m2 * v2_after_lab[i], 3)
    @test after_energy_lab ≈ before_energy rtol=1e-13
    @test collect(after_momentum_lab) ≈ collect(before_momentum)
    measured_loss_lab = 1 - sum(abs2, v1_after_lab) / sum(abs2, v1)
    @test measured_loss_lab ≈
          fractional_energy_loss_lab(theta_lab, m1, m2)
    @test maximum_lab_scattering_angle(16.0, 16.0) ≈ pi / 2
    @test maximum_lab_scattering_angle(16.0, 44.0) ≈ pi
    @test MarsHotO.lab_to_com_scattering_angle(
        deg2rad(60.0), 16.0, 16.0,
    ) ≈ deg2rad(120.0)
    theta_com_from_lab = MarsHotO.lab_to_com_scattering_angle(
        theta_lab, m1, m2,
    )
    @test fractional_energy_loss_lab(theta_lab, m1, m2) ≈
          fractional_energy_loss(theta_com_from_lab, m1, m2)
    @test_throws DomainError fractional_energy_loss_lab(
        deg2rad(91.0), 16.0, 16.0,
    )
end

@testset "MGITM atmosphere" begin
    profile = load_mgitm_subsolar_profile(joinpath(
        ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
    ))
    @test issorted(profile.altitude_m)
    @test profile.latitude_deg == [-2.5, 2.5]
    @test profile.longitude_deg == [27.5, 27.5]
    @test profile.sza_deg ≈ 2.73115946 atol=1e-7
    @test all(profile.density_m3[:O2p] .> 0)
    @test interpolate_profile(profile, 500e3).density_m3[:CO2] <
          profile.density_m3[:CO2][end]

    branches = load_reaction_branches(joinpath(
        ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
    ))
    targets = load_collision_targets(joinpath(
        ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
    ))
    rng = Xoshiro(20260727)
    source = sample_hot_o_source(
        rng, (MARS_RADIUS_M + 200e3, 0.0, 0.0), 1500.0, 300.0, branches,
    )
    @test source.weight == 1.0
    @test sum(abs2, source.velocity_m_s) > 0
    result = transport_particle!(
        rng, source, profile, targets; max_steps=10, step_m=100.0,
    )
    @test result.reason == :maximum_steps
    @test result.steps == 10

    @test rahmati_step_length(5_000.0) == 500.0
    @test rahmati_step_length(10_000.0) == 1000.0
    corona = run_hot_o_corona(
        profile, targets, branches;
        chemistry_path=joinpath(
            ROOT, "data", "chemistry",
            "o2plus_dissociative_recombination.toml",
        ),
        config=RahmatiMonteCarloConfig(
            primary_particles=20,
            seed=73,
            maximum_altitude_m=260e3,
            maximum_steps_per_particle=10_000,
            maximum_total_particles=10_000,
            altitude_edges_km=collect(100.0:10.0:260.0),
            energy_edges_eV=collect(range(0.01, 7.0; length=29)),
        ),
    )
    @test corona.primary_particles == 20
    @test all(corona.density_cm3_eV1 .>= 0)
    @test sum(corona.density_cm3_eV1) > 0
end
