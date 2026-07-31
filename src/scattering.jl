"""
Rahmati analytical COM scattering-angle distribution.

The probability density includes the solid-angle Jacobian sin(theta). The
complete interval from 0 to pi is sampled without a minimum-angle cutoff.
"""
const RAHMATI_DCS_ALPHA_M2_SR = 0.36e-16 * 1e-4
const RAHMATI_DCS_BETA = -1.85

"""Rahmati fit to the Kharchenko O-O COM differential cross section."""
function differential_cross_section(
    theta_com_rad::Real;
    alpha_m2_sr=RAHMATI_DCS_ALPHA_M2_SR,
    beta=RAHMATI_DCS_BETA,
)
    0 < theta_com_rad <= pi ||
        throw(DomainError(theta_com_rad, "Require 0 < theta <= pi"))
    alpha_m2_sr * sin(theta_com_rad / 2)^beta
end

"""
Normalized COM polar-angle PDF including the solid-angle Jacobian.

The default theta_min is zero, so scattering below 10 degrees is retained.
"""
function scattering_angle_pdf(
    theta_com_rad::Real;
    theta_min_rad=0.0,
    beta=RAHMATI_DCS_BETA,
)
    theta_min_rad <= theta_com_rad <= pi || return 0.0
    exponent = beta + 2
    exponent > 0 ||
        throw(DomainError(beta, "Angular distribution is not integrable"))
    lower = sin(theta_min_rad / 2)^exponent
    exponent * sin(theta_com_rad / 2)^(beta + 1) *
        cos(theta_com_rad / 2) / (2 * (1 - lower))
end

"""CDF of the Rahmati COM scattering-angle distribution."""
function scattering_angle_cdf(
    theta_com_rad::Real;
    theta_min_rad=0.0,
    beta=RAHMATI_DCS_BETA,
)
    0 <= theta_min_rad < pi ||
        throw(DomainError(theta_min_rad, "Require 0 <= theta_min < pi"))
    theta_com_rad <= theta_min_rad && return 0.0
    theta_com_rad >= pi && return 1.0
    exponent = beta + 2
    exponent > 0 ||
        throw(DomainError(beta, "Angular distribution is not integrable"))
    lower = sin(theta_min_rad / 2)^exponent
    (sin(theta_com_rad / 2)^exponent - lower) / (1 - lower)
end

"""
Inverse-CDF sample of the Rahmati COM scattering angle.

By default theta_min = 0, retaining the complete forward peak.
"""
function sample_scattering_angle(
    rng=Random.default_rng();
    theta_min_rad=0.0,
    beta=RAHMATI_DCS_BETA,
)
    0 <= theta_min_rad < pi ||
        throw(DomainError(theta_min_rad, "Require 0 <= theta_min < pi"))
    exponent = beta + 2
    exponent > 0 ||
        throw(DomainError(beta, "Angular distribution is not integrable"))
    lower = sin(theta_min_rad / 2)^exponent
    2asin((lower + rand(rng) * (1 - lower))^(1 / exponent))
end

"""Inverse-CDF scattering angle for an explicitly supplied uniform variate."""
function sample_scattering_angle_from_uniform(
    u::Real; theta_min_rad=0.0, beta=RAHMATI_DCS_BETA,
)
    0 <= u < 1 ||
        throw(DomainError(u, "Scattering uniform must be in [0, 1)"))
    0 <= theta_min_rad < pi ||
        throw(DomainError(theta_min_rad, "Require 0 <= theta_min < pi"))
    exponent = beta + 2
    exponent > 0 ||
        throw(DomainError(beta, "Angular distribution is not integrable"))
    lower = sin(theta_min_rad / 2)^exponent
    2asin((lower + u * (1 - lower))^(1 / exponent))
end

sample_azimuth(rng=Random.default_rng()) = 2pi * rand(rng)
