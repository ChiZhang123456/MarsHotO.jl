using MarsHotO
using Random

const ROOT = normpath(joinpath(@__DIR__, ".."))
const SEED = 20260731
const SOURCE_ALTITUDE_KM = 140.0

atmosphere = load_mgitm_subsolar_profile(joinpath(
    ROOT, "MGITM", "MGITM_LS000_F070_150901.dat",
))
chemistry_path = joinpath(
    ROOT, "data", "chemistry", "o2plus_dissociative_recombination.toml",
)
branches = load_reaction_branches(chemistry_path)
targets = load_collision_targets(joinpath(
    ROOT, "data", "cross_sections", "rahmati_total_cross_sections.toml",
))
vibration_probability, vibration_quantum_eV =
    MarsHotO._load_vibration(chemistry_path)

position_m = (MARS_RADIUS_M + 1000SOURCE_ALTITUDE_KM, 0.0, 0.0)
state = interpolate_profile(atmosphere, 1000SOURCE_ALTITUDE_KM)
event_rate_m3_s1 = dissociative_recombination_event_rate(
    state.density_m3[:e], state.density_m3[:O2p], state.Te_K,
)
radius_inner_m = MARS_RADIUS_M + 1000(SOURCE_ALTITUDE_KM - 0.5)
radius_outer_m = MARS_RADIUS_M + 1000(SOURCE_ALTITUDE_KM + 0.5)
shell_volume_m3 = (4pi / 3) * (radius_outer_m^3 - radius_inner_m^3)
event_weight_s1 = event_rate_m3_s1 * shell_volume_m3
rng = Xoshiro(SEED)

event = sample_dissociative_recombination_event(
    rng, position_m, state.Te_K, state.Ti_K, branches;
    vibrational_probability=vibration_probability,
    vibrational_quantum_eV=vibration_quantum_eV,
    weight_s1=event_weight_s1,
    plasma_bulk_velocity_m_s=(0.0, 0.0, 0.0),
)

println("seed = ", SEED)
println("source altitude = ", SOURCE_ALTITUDE_KM, " km")
println("Te = ", state.Te_K, " K, Ti = ", state.Ti_K, " K")
println("ne = ", state.density_m3[:e], " m^-3")
println("nO2+ = ", state.density_m3[:O2p], " m^-3")
println("DR event rate = ", event_rate_m3_s1, " m^-3 s^-1")
println("event weight for one event in the 1 km shell = ",
        event.weight_s1, " s^-1")
println("electron velocity = ", event.electron_velocity_m_s, " m/s")
println("O2+ velocity = ", event.o2p_velocity_m_s, " m/s")
println("event COM velocity = ", event.com_velocity_m_s, " m/s")
println("event COM speed = ", MarsHotO._norm(event.com_velocity_m_s), " m/s")
println("branch = ", event.branch.products)
println("vibrational level = ", event.vibrational_level)
println("relative energy = ", event.relative_energy_eV, " eV")
println("available energy = ", event.available_energy_eV, " eV")
println("O1 initial velocity = ", event.products[1].velocity_m_s, " m/s")
println("O2 initial velocity = ", event.products[2].velocity_m_s, " m/s")

o1_com = MarsHotO._subtract(
    event.products[1].velocity_m_s, event.com_velocity_m_s,
)
o2_com = MarsHotO._subtract(
    event.products[2].velocity_m_s, event.com_velocity_m_s,
)
println("O1 velocity in event COM = ", o1_com, " m/s")
println("O2 velocity in event COM = ", o2_com, " m/s")
println("COM opposition residual |v1_COM + v2_COM| = ",
        MarsHotO._norm(MarsHotO._add(o1_com, o2_com)), " m/s")

particle = event.products[1]
step = advance_hot_o_step!(rng, particle, atmosphere, targets)
println("\nOne advance_hot_o_step! result")
println("ds = ", step.ds_m, " m, dt = ", step.dt_s, " s")
println("collision coefficient = ", step.collision_coefficient_m1, " m^-1")
println("collision probability = ",
        min(step.ds_m * step.collision_coefficient_m1, 1.0))
println("collision uniform = ", step.collision_uniform)
println("collision occurred = ", !isnothing(step.target))
println("target uniform = ", step.target_uniform)
println("target = ", isnothing(step.target) ? nothing : step.target.species)
println("scattering uniform = ", step.scattering_uniform)
println("COM scattering angle = ", step.scattering_angle_com_rad, " rad")
println("azimuth uniform = ", step.azimuth_uniform)
println("final position = ", particle.position_m, " m")
println("final velocity = ", particle.velocity_m_s, " m/s")
