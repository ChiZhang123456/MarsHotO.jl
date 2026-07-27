const DEFAULT_SCATTERING_ANGLE_PATH = normpath(joinpath(
    @__DIR__, "..", "data", "cross_sections",
    "scattering_angle_distribution.txt",
))
const _DEFAULT_SCATTERING_ANGLE_DISTRIBUTION =
    Ref{Union{Nothing,ScatteringAngleDistribution}}(nothing)

"""
Read the Kallio and Barabash (2001) inverse-CDF lookup table.

The first column is a uniform random number and the second column is the
projectile scattering angle in the LAB frame, in degrees.
"""
function load_scattering_angle_distribution(path::AbstractString)
    random_number = Float64[]
    theta_lab_rad = Float64[]
    for line in eachline(path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, '#') ||
         startswith(stripped, "random_number")) && continue
        fields = split(stripped)
        length(fields) >= 2 || continue
        push!(random_number, parse(Float64, fields[1]))
        push!(theta_lab_rad, deg2rad(parse(Float64, fields[2])))
    end
    length(random_number) >= 2 ||
        error("Scattering-angle table must contain at least two rows")
    random_number[1] = 0.0
    issorted(random_number) ||
        error("Scattering-angle random-number column must be sorted")
    issorted(theta_lab_rad) ||
        error("Scattering-angle column must be sorted")
    isapprox(random_number[1], 0.0; atol=1e-12) ||
        error("Scattering-angle random-number grid must start at zero")
    isapprox(random_number[end], 1.0; atol=1e-12) ||
        error("Scattering-angle random-number grid must end at one")
    ScatteringAngleDistribution(random_number, theta_lab_rad)
end

function default_scattering_angle_distribution()
    if isnothing(_DEFAULT_SCATTERING_ANGLE_DISTRIBUTION[])
        _DEFAULT_SCATTERING_ANGLE_DISTRIBUTION[] =
            load_scattering_angle_distribution(DEFAULT_SCATTERING_ANGLE_PATH)
    end
    _DEFAULT_SCATTERING_ANGLE_DISTRIBUTION[]::ScatteringAngleDistribution
end

@inline function _linear_lookup(x, xp, fp)
    x <= xp[1] && return fp[1]
    x >= xp[end] && return fp[end]
    index = searchsortedlast(xp, x)
    weight = (x - xp[index]) / (xp[index + 1] - xp[index])
    muladd(weight, fp[index + 1] - fp[index], fp[index])
end

"""Draw one LAB scattering angle from the tabulated inverse CDF."""
sample_scattering_angle(
    rng=Random.default_rng(),
    distribution::ScatteringAngleDistribution=
        default_scattering_angle_distribution(),
) = _linear_lookup(
    rand(rng), distribution.random_number, distribution.theta_lab_rad,
)

"""
Draw a kinematically accessible LAB projectile angle for a stationary target.
For equal or lighter targets, the inverse CDF is conditioned on the maximum
LAB angle permitted by elastic two-body kinematics.
"""
function sample_scattering_angle(
    rng,
    distribution::ScatteringAngleDistribution,
    projectile_mass_kg::Real,
    target_mass_kg::Real,
)
    theta_max = maximum_lab_scattering_angle(
        projectile_mass_kg, target_mass_kg,
    )
    cumulative_max = scattering_angle_cdf(theta_max, distribution)
    _linear_lookup(
        rand(rng) * cumulative_max,
        distribution.random_number,
        distribution.theta_lab_rad,
    )
end

"""Return the tabulated cumulative probability at a LAB scattering angle."""
scattering_angle_cdf(
    theta_lab_rad::Real,
    distribution::ScatteringAngleDistribution=
        default_scattering_angle_distribution(),
) = _linear_lookup(
    theta_lab_rad, distribution.theta_lab_rad, distribution.random_number,
)

sample_azimuth(rng=Random.default_rng()) = 2pi * rand(rng)
