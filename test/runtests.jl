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

@testset "Maxwellian thermal source sampling" begin
    temperature_K = 300.0
    mass_kg = O_MASS_KG
    thermal_speed = maxwellian_thermal_speed(temperature_K, mass_kg)
    component_variance =
        MarsHotO.BOLTZMANN_J_K * temperature_K / mass_kg

    speed_grid = range(0.0, 10thermal_speed; length=200_001)
    speed_pdf = [
        4pi * speed^2 * maxwellian_velocity_pdf(
            (speed, 0.0, 0.0), temperature_K, mass_kg,
        )
        for speed in speed_grid
    ]
    normalization = sum(
        (speed_pdf[i] + speed_pdf[i + 1]) / 2
        for i in 1:(length(speed_pdf) - 1)
    ) * step(speed_grid)
    @test normalization ≈ 1.0 atol=1e-10

    bulk_velocity = (120.0, -45.0, 30.0)
    rng = Xoshiro(73)
    sample_count = 200_000
    velocities = [
        sample_maxwellian_velocity(
            rng, temperature_K, mass_kg;
            bulk_velocity_m_s=bulk_velocity,
        )
        for _ in 1:sample_count
    ]
    mean_velocity = [
        sum(velocity[component] for velocity in velocities) / sample_count
        for component in 1:3
    ]
    component_sigma = sqrt(component_variance)
    @test all(
        abs.(mean_velocity .- collect(bulk_velocity)) .<
        0.01component_sigma
    )
    sampled_component_variance = [
        sum(
            (velocity[component] - bulk_velocity[component])^2
            for velocity in velocities
        ) / sample_count
        for component in 1:3
    ]
    @test all(
        isapprox.(sampled_component_variance, component_variance; rtol=0.01)
    )
    sampled_mean_energy = sum(
        0.5mass_kg * sum(
            (velocity[component] - bulk_velocity[component])^2
            for component in 1:3
        )
        for velocity in velocities
    ) / sample_count
    expected_mean_energy =
        1.5 * MarsHotO.BOLTZMANN_J_K * temperature_K
    @test sampled_mean_energy ≈ expected_mean_energy rtol=0.01

    rng1 = Xoshiro(20260728)
    rng2 = Xoshiro(20260728)
    @test sample_maxwellian_velocity(rng1, temperature_K, mass_kg) ==
          sample_maxwellian_velocity(rng2, temperature_K, mass_kg)
    @test_throws DomainError maxwellian_thermal_speed(0.0, mass_kg)
    @test_throws DomainError maxwellian_velocity_pdf(
        (0.0, 0.0, 0.0), temperature_K, -mass_kg,
    )
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
        rng, (MARS_RADIUS_M + 200e3, 0.0, 0.0), 1500.0, 300.0, branches;
        weight_s1=2.5e20,
    )
    @test source.weight_s1 == 2.5e20
    @test sum(abs2, source.velocity_m_s) > 0
    result = transport_particle!(
        rng, source, profile, targets; max_steps=10, step_m=100.0,
    )
    @test result.reason == :maximum_steps
    @test result.steps == 10

    @test rahmati_step_length(5_000.0) == 500.0
    @test rahmati_step_length(10_000.0) == 1000.0
    corona_config = RahmatiMonteCarloConfig(
        particles_per_source_altitude=10,
        source_altitudes_km=[150.0, 160.0],
        seed=73,
        maximum_altitude_m=260e3,
        maximum_steps_per_particle=10_000,
        maximum_total_particles=10_000,
        altitude_edges_km=collect(100.0:10.0:260.0),
        energy_edges_eV=collect(range(0.01, 7.0; length=29)),
    )
    corona = run_hot_o_corona(
        profile, targets, branches;
        chemistry_path=joinpath(
            ROOT, "data", "chemistry",
            "o2plus_dissociative_recombination.toml",
        ),
        config=corona_config,
    )
    @test corona.primary_particles == 20
    @test corona.particles_per_source_altitude == 10
    @test corona.source_altitudes_km == [150.0, 160.0]
    @test all(corona.source_particle_weights_s1 .> 0)
    @test sum(corona.source_particle_weights_s1) * 10 ≈
          corona.total_source_rate_s1
    @test all(corona.density_m3_per_bin .>= 0)
    @test sum(corona.density_m3_per_bin) > 0

    source_altitudes_m = 1000 .* corona.source_altitudes_km
    source_edges_m = MarsHotO._shell_edges_m(source_altitudes_m)
    source_shell_volume_m3 = (4pi / 3) .* (
        (MARS_RADIUS_M .+ source_edges_m[2:end]).^3 .-
        (MARS_RADIUS_M .+ source_edges_m[1:end-1]).^3
    )
    for i in eachindex(source_altitudes_m)
        state = interpolate_profile(profile, source_altitudes_m[i])
        expected_weight_s1 = hot_o_production_rate(
            state.density_m3[:e], state.density_m3[:O2p], state.Te_K,
        ) * source_shell_volume_m3[i] /
            corona.particles_per_source_altitude
        @test corona.source_particle_weights_s1[i] ≈ expected_weight_s1
    end

    corona_repeat = run_hot_o_corona(
        profile, targets, branches;
        chemistry_path=joinpath(
            ROOT, "data", "chemistry",
            "o2plus_dissociative_recombination.toml",
        ),
        config=corona_config,
    )
    @test corona_repeat.source_particle_weights_s1 ==
          corona.source_particle_weights_s1
    @test corona_repeat.density_m3_per_bin == corona.density_m3_per_bin
    @test corona_repeat.stop_counts == corona.stop_counts
end
