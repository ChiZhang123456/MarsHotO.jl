module MarsHotO

using DelimitedFiles
using Random
using TOML

include("shared/constants.jl")
include("shared/types.jl")
include("shared/atmosphere.jl")
include("shared/chemistry.jl")
include("shared/cross_sections.jl")
include("shared/scattering.jl")
include("shared/collision_kinematics.jl")

include("monte_carlo/types.jl")
include("monte_carlo/source_particles.jl")
include("monte_carlo/transport.jl")
include("monte_carlo/ensembles.jl")

include("two_fluid/types.jl")
include("two_fluid/two_stream.jl")

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
export maximum_lab_scattering_angle
export elastic_collision, elastic_collision_lab
export sample_hot_o_source, transport_particle!
export RahmatiMonteCarloConfig, HotOCoronaResult, rahmati_step_length
export run_hot_o_corona, write_corona_distribution
export TwoStreamConfig, TwoStreamRedistribution, TwoStreamResult
export build_two_stream_redistribution, run_two_stream, write_two_stream_flux

end
