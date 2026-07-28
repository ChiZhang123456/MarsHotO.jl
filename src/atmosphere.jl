"""
MGITM atmosphere input and interpolation.

The nearest-subsolar grid columns are averaged at each altitude. Number
densities remain in m^-3, temperatures in K, and altitude is converted to m.
"""
const MGITM_COLUMNS = (
    :longitude_deg, :latitude_deg, :altitude_km, :Tn_K, :Ti_K, :Te_K,
    :CO2, :O, :N2, :CO, :O2, :O2p, :Op, :CO2p, :e, :UN, :VN, :WN,
)

function _subsolar_longitude(path::AbstractString)
    for line in eachline(path)
        match_value = match(r"subsolar_longitude:\s*([+-]?\d+(?:\.\d+)?)", line)
        isnothing(match_value) || return parse(Float64, match_value.captures[1])
    end
    error("Subsolar longitude was not found in $path")
end

@inline _angular_distance_deg(a, b) = mod(a - b + 180, 360) - 180

"""
Read the MGITM columns nearest the subsolar point and average tied latitude
columns. Number densities remain in m^-3 and altitude is converted to metres.
"""
function load_mgitm_subsolar_profile(path::AbstractString)
    raw = readdlm(path, comments=true, comment_char='#')
    size(raw, 2) == length(MGITM_COLUMNS) ||
        error("Expected $(length(MGITM_COLUMNS)) MGITM columns, found $(size(raw, 2))")
    subsolar_lon = _subsolar_longitude(path)
    lon, lat = raw[:, 1], raw[:, 2]
    grid_pairs = unique(collect(zip(lon, lat)))
    sza(pair) = rad2deg(acos(clamp(
        cosd(pair[2]) * cosd(_angular_distance_deg(pair[1], subsolar_lon)),
        -1, 1,
    )))
    minimum_sza = minimum(sza.(grid_pairs))
    selected = filter(
        pair -> isapprox(sza(pair), minimum_sza; atol=1e-10), grid_pairs,
    )
    keep = [any(pair -> row[1] == pair[1] && row[2] == pair[2], selected)
            for row in eachrow(raw)]
    chosen = raw[keep, :]
    altitudes = sort(unique(chosen[:, 3]))
    averaged = zeros(length(altitudes), size(chosen, 2))
    for (i, altitude) in pairs(altitudes)
        rows = chosen[chosen[:, 3] .== altitude, :]
        averaged[i, :] .= vec(sum(rows; dims=1)) ./ size(rows, 1)
    end
    density = Dict{Symbol,Vector{Float64}}(
        species => vec(averaged[:, column])
        for (species, column) in
        ((:CO2, 7), (:O, 8), (:N2, 9), (:CO, 10), (:O2, 11),
         (:O2p, 12), (:Op, 13), (:CO2p, 14), (:e, 15))
    )
    AtmosphereProfile(
        vec(averaged[:, 3]) .* 1000,
        vec(averaged[:, 4]), vec(averaged[:, 5]), vec(averaged[:, 6]),
        density, first.(selected), last.(selected), minimum_sza,
    )
end

function _linear_interp(x::Real, xp, fp)
    x <= xp[1] && return fp[1]
    x >= xp[end] && return fp[end]
    i = searchsortedlast(xp, x)
    w = (x - xp[i]) / (xp[i + 1] - xp[i])
    return muladd(w, fp[i + 1] - fp[i], fp[i])
end

function _log_density_interp(x::Real, xp, fp)
    positive = max.(fp, floatmin(Float64))
    if x <= xp[1]
        return positive[1]
    elseif x >= xp[end]
        slope = (log(positive[end]) - log(positive[end-1])) /
                (xp[end] - xp[end-1])
        return exp(log(positive[end]) + min(slope, 0.0) * (x - xp[end]))
    end
    exp(_linear_interp(x, xp, log.(positive)))
end

"""Interpolate temperatures and densities at altitude in metres."""
function interpolate_profile(profile::AtmosphereProfile, altitude_m::Real)
    densities = Dict(k => _log_density_interp(altitude_m, profile.altitude_m, v)
                     for (k, v) in profile.density_m3)
    return (
        Tn_K=_linear_interp(altitude_m, profile.altitude_m, profile.Tn_K),
        Ti_K=_linear_interp(altitude_m, profile.altitude_m, profile.Ti_K),
        Te_K=_linear_interp(altitude_m, profile.altitude_m, profile.Te_K),
        density_m3=densities,
    )
end
