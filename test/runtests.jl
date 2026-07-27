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
        ROOT, "data", "cross_sections", "ali_total_cross_sections.toml",
    ))
    oxygen = only(filter(x -> x.species == :O, targets))
    @test total_cross_section(oxygen, 3.0) ≈ 6.4e-19
    @test total_cross_section(oxygen, 12.0) <
          total_cross_section(oxygen, 3.0)
end

@testset "Scattering distribution" begin
    theta_min = deg2rad(10.0)
    @test scattering_angle_cdf(theta_min; theta_min_rad=theta_min) == 0.0
    @test scattering_angle_cdf(pi; theta_min_rad=theta_min) == 1.0
    @test 0.25 < angular_cross_section_fraction(theta_min) < 0.35
    theta = range(theta_min, pi; length=100_001)
    pdf = scattering_angle_pdf.(theta; theta_min_rad=theta_min)
    integral = sum((pdf[1:end-1] .+ pdf[2:end]) .* diff(theta) ./ 2)
    @test integral ≈ 1.0 atol=2e-6

    rng = Xoshiro(73)
    samples = [sample_scattering_angle(rng; theta_min_rad=theta_min)
               for _ in 1:100_000]
    @test all(theta_min .<= samples .<= pi)
    @test abs(sum(x <= pi / 2 for x in samples) / length(samples) -
              scattering_angle_cdf(pi / 2; theta_min_rad=theta_min)) < 0.005
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
        ROOT, "data", "cross_sections", "ali_total_cross_sections.toml",
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

    @test ali_step_length(5_000.0) == 500.0
    @test ali_step_length(10_000.0) == 1000.0
    corona = run_hot_o_corona(
        profile, targets, branches;
        chemistry_path=joinpath(
            ROOT, "data", "chemistry",
            "o2plus_dissociative_recombination.toml",
        ),
        config=AliMonteCarloConfig(
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
