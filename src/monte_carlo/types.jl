mutable struct HotOParticle
    position_m::NTuple{3,Float64}
    velocity_m_s::NTuple{3,Float64}
    weight::Float64
    alive::Bool
    collisions::Int
end
