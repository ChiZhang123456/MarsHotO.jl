"""
Core atmosphere, chemistry, and collision data structures.

Field names include units wherever practical so that transport inputs remain
traceable across file reading, interpolation, and collision calculations.
"""
struct AtmosphereProfile
    altitude_m::Vector{Float64}
    Tn_K::Vector{Float64}
    Ti_K::Vector{Float64}
    Te_K::Vector{Float64}
    density_m3::Dict{Symbol,Vector{Float64}}
    longitude_deg::Vector{Float64}
    latitude_deg::Vector{Float64}
    sza_deg::Float64
end

struct DRBranch
    products::String
    release_energy_eV::Float64
    probability::Float64
end

struct CollisionTarget
    species::Symbol
    mass_kg::Float64
    sigma_3eV_m2::Float64
end

"""
One weighted hot O test particle propagated in the stationary Mars frame.

`weight_s1` is the physical hot O production rate represented by this
macroparticle, in s^-1.
"""
mutable struct HotOParticle
    position_m::NTuple{3,Float64}
    velocity_m_s::NTuple{3,Float64}
    weight_s1::Float64
    alive::Bool
    collisions::Int
end

"""One weighted O2+ dissociative-recombination event and its paired products."""
struct DissociativeRecombinationEvent
    position_m::NTuple{3,Float64}
    weight_s1::Float64
    electron_velocity_m_s::NTuple{3,Float64}
    o2p_velocity_m_s::NTuple{3,Float64}
    com_velocity_m_s::NTuple{3,Float64}
    relative_energy_eV::Float64
    branch::DRBranch
    vibrational_level::Int
    available_energy_eV::Float64
    products::NTuple{2,HotOParticle}
end

"""Random choices and outcome of one numerical hot O transport step."""
struct HotOTransportStep
    position_before_m::NTuple{3,Float64}
    velocity_before_m_s::NTuple{3,Float64}
    position_after_m::NTuple{3,Float64}
    ballistic_velocity_after_m_s::NTuple{3,Float64}
    velocity_after_m_s::NTuple{3,Float64}
    dt_s::Float64
    ds_m::Float64
    collision_coefficient_m1::Float64
    collision_uniform::Float64
    target_uniform::Union{Nothing,Float64}
    scattering_uniform::Union{Nothing,Float64}
    azimuth_uniform::Union{Nothing,Float64}
    target::Union{Nothing,CollisionTarget}
    scattering_angle_com_rad::Union{Nothing,Float64}
    secondary::Union{Nothing,HotOParticle}
end
