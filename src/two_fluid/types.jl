Base.@kwdef struct TwoStreamConfig
    altitude_min_m::Float64 = 100e3
    altitude_max_m::Float64 = 300e3
    altitude_step_m::Float64 = 1e3
    energy_edges_eV::Vector{Float64} = collect(range(0.01, 7.01; length=71))
    mean_pitch_cosine::Float64 = 0.5
    redistribution_samples::Int = 20_000
    redistribution_seed::Int = 73
    maximum_iterations::Int = 2_000
    relative_tolerance::Float64 = 1e-3
    top_scale_height_m::Union{Nothing,Float64} = nothing
end

struct TwoStreamRedistribution
    same_stream::Array{Float64,3}
    reverse_stream::Array{Float64,3}
    secondary_same_stream::Array{Float64,3}
    secondary_reverse_stream::Array{Float64,3}
end

struct TwoStreamResult
    altitude_m::Vector{Float64}
    energy_edges_eV::Vector{Float64}
    energy_centers_eV::Vector{Float64}
    upward_flux_m2_s1::Matrix{Float64}
    downward_flux_m2_s1::Matrix{Float64}
    primary_production_m3_s1::Matrix{Float64}
    escape_energy_eV::Float64
    escape_flux_m2_s1::Float64
    iterations::Int
    converged::Bool
end
