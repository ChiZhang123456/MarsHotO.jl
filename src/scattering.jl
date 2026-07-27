const RAHMATI_DCS_ALPHA_M2_SR = 0.36e-16 * 1e-4
const RAHMATI_DCS_BETA = -1.85

"""Fraction of the fitted angular cross section retained above theta_min."""
function angular_cross_section_fraction(theta_min_rad::Real;
                                        beta=RAHMATI_DCS_BETA)
    0 <= theta_min_rad < pi ||
        throw(DomainError(theta_min_rad, "Require 0 <= theta_min < pi"))
    1 - sin(theta_min_rad / 2)^(beta + 2)
end

"""Kharchenko O-O differential cross section fitted by Rahmati."""
function differential_cross_section(theta_rad::Real;
                                    alpha_m2_sr=RAHMATI_DCS_ALPHA_M2_SR,
                                    beta=RAHMATI_DCS_BETA)
    0 < theta_rad <= pi || throw(DomainError(theta_rad, "Require 0 < theta <= pi"))
    alpha_m2_sr * sin(theta_rad / 2)^beta
end

"""Normalized polar-angle PDF, including the sin(theta) solid-angle Jacobian."""
function scattering_angle_pdf(theta_rad::Real; theta_min_rad=deg2rad(10.0),
                              beta=RAHMATI_DCS_BETA)
    theta_min_rad <= theta_rad <= pi || return 0.0
    p = beta + 2
    denominator = 1 - sin(theta_min_rad / 2)^p
    p * sin(theta_rad / 2)^(beta + 1) * cos(theta_rad / 2) /
        (2denominator)
end

function scattering_angle_cdf(theta_rad::Real; theta_min_rad=deg2rad(10.0),
                              beta=RAHMATI_DCS_BETA)
    theta_rad <= theta_min_rad && return 0.0
    theta_rad >= pi && return 1.0
    p = beta + 2
    lower = sin(theta_min_rad / 2)^p
    (sin(theta_rad / 2)^p - lower) / (1 - lower)
end

"""Inverse-CDF sample of the COM polar angle."""
function sample_scattering_angle(rng=Random.default_rng();
                                 theta_min_rad=deg2rad(10.0),
                                 beta=RAHMATI_DCS_BETA)
    0 <= theta_min_rad < pi ||
        throw(DomainError(theta_min_rad, "Require 0 <= theta_min < pi"))
    p = beta + 2
    p > 0 || throw(DomainError(beta, "The angular distribution is not integrable"))
    lower = sin(theta_min_rad / 2)^p
    2asin((lower + rand(rng) * (1 - lower))^(1 / p))
end

sample_azimuth(rng=Random.default_rng()) = 2pi * rand(rng)
