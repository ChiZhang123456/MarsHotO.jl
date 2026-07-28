"""
MarsHotO models photochemical hot oxygen at Mars.

The package exposes low-level atmosphere, chemistry, scattering, and
two-body-kinematics functions, a direct single-particle propagator, and a
weighted particle-ensemble driver with residence-time density diagnostics.
All internal transport calculations use SI units unless a field or function
name explicitly states eV, km, cm^-3, or another output unit.
"""
module MarsHotO

using DelimitedFiles
using Random
using TOML

include("constants.jl")
include("types.jl")
include("atmosphere.jl")
include("chemistry.jl")
include("cross_sections.jl")
include("scattering.jl")
include("collision_kinematics.jl")
include("source_particles.jl")
include("transport.jl")
include("ensembles.jl")

export AtmosphereProfile, CollisionTarget, DRBranch, HotOParticle
export MARS_RADIUS_M, O_MASS_KG, EV_J
export load_mgitm_subsolar_profile, interpolate_profile
export dissociative_recombination_coefficient, hot_o_production_rate
export load_reaction_branches, load_collision_targets
export total_cross_section, collision_coefficient, choose_collision_target
export differential_cross_section, scattering_angle_pdf
export scattering_angle_cdf
export sample_scattering_angle, sample_azimuth
export fractional_energy_loss, fractional_energy_loss_lab
export maximum_lab_scattering_angle
export elastic_collision, elastic_collision_lab
export maxwellian_thermal_speed, maxwellian_velocity_pdf
export sample_maxwellian_velocity
export sample_hot_o_source, transport_particle!
export RahmatiMonteCarloConfig, HotOCoronaResult, rahmati_step_length
export run_hot_o_corona, write_corona_distribution

end
