@inline _add(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _subtract(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _scale(c, a) = (c * a[1], c * a[2], c * a[3])
@inline _dot(a, b) = a[1]b[1] + a[2]b[2] + a[3]b[3]
@inline _norm(a) = sqrt(_dot(a, a))
@inline _cross(a, b) = (
    a[2]b[3] - a[3]b[2], a[3]b[1] - a[1]b[3], a[1]b[2] - a[2]b[1],
)

fractional_energy_loss(theta_com_rad::Real, projectile_mass_kg::Real,
                       target_mass_kg::Real) =
    2projectile_mass_kg * target_mass_kg /
    (projectile_mass_kg + target_mass_kg)^2 * (1 - cos(theta_com_rad))

function _rotate_relative_velocity(relative_velocity, theta, phi)
    e0 = _scale(1 / _norm(relative_velocity), relative_velocity)
    helper = abs(e0[3]) > 0.9 ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0)
    e1raw = _cross(helper, e0)
    e1 = _scale(1 / _norm(e1raw), e1raw)
    e2 = _cross(e0, e1)
    speed = _norm(relative_velocity)
    _scale(speed, _add(_scale(cos(theta), e0),
        _scale(sin(theta), _add(_scale(cos(phi), e1), _scale(sin(phi), e2)))))
end

"""
Elastic two-body collision. Angles are defined in the COM frame. Returns the
post-collision projectile and target LAB velocities.
"""
function elastic_collision(projectile_velocity, target_velocity,
                           projectile_mass_kg, target_mass_kg,
                           theta_com_rad, phi_rad)
    total_mass = projectile_mass_kg + target_mass_kg
    vcom = _scale(1 / total_mass,
        _add(_scale(projectile_mass_kg, projectile_velocity),
             _scale(target_mass_kg, target_velocity)))
    relative_before = _subtract(projectile_velocity, target_velocity)
    relative_after = _rotate_relative_velocity(relative_before, theta_com_rad, phi_rad)
    projectile_after = _add(vcom, _scale(target_mass_kg / total_mass, relative_after))
    target_after = _subtract(vcom, _scale(projectile_mass_kg / total_mass, relative_after))
    projectile_after, target_after
end
