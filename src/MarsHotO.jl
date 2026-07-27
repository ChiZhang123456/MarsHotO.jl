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
export ScatteringAngleDistribution
export MARS_RADIUS_M, O_MASS_KG, EV_J
export load_mgitm_subsolar_profile, interpolate_profile
export dissociative_recombination_coefficient, hot_o_production_rate
export load_reaction_branches, load_collision_targets
export total_cross_section, collision_coefficient, choose_collision_target
export load_scattering_angle_distribution, scattering_angle_cdf
export sample_scattering_angle, sample_azimuth
export fractional_energy_loss, fractional_energy_loss_lab
export elastic_collision, elastic_collision_lab
export sample_hot_o_source, transport_particle!
export RahmatiMonteCarloConfig, HotOCoronaResult, rahmati_step_length
export run_hot_o_corona, write_corona_distribution

end
