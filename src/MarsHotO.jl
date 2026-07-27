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

export AtmosphereProfile, CollisionTarget, DRBranch, HotOParticle
export MARS_RADIUS_M, O_MASS_KG, EV_J
export load_mgitm_subsolar_profile, interpolate_profile
export dissociative_recombination_coefficient, hot_o_production_rate
export load_reaction_branches, load_collision_targets
export total_cross_section, collision_coefficient, choose_collision_target
export differential_cross_section, scattering_angle_pdf
export angular_cross_section_fraction, scattering_angle_cdf
export sample_scattering_angle, sample_azimuth
export fractional_energy_loss, elastic_collision
export sample_hot_o_source, transport_particle!

end
