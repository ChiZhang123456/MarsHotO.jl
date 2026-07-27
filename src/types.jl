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

struct ScatteringAngleDistribution
    random_number::Vector{Float64}
    theta_lab_rad::Vector{Float64}
end

mutable struct HotOParticle
    position_m::NTuple{3,Float64}
    velocity_m_s::NTuple{3,Float64}
    weight::Float64
    alive::Bool
    collisions::Int
end
