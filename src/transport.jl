@inline kinetic_energy_eV(velocity, mass=O_MASS_KG) =
    0.5mass * _dot(velocity, velocity) / EV_J

"""
Transport one hot O particle in a spherically symmetric MGITM atmosphere.
Collisions are sampled by optical depth. Gravity is updated at every spatial
step. Returns a named tuple with the final state and stopping reason.
"""
function transport_particle!(rng, particle::HotOParticle,
                             atmosphere::AtmosphereProfile, targets;
                             scattering_distribution=
                                 default_scattering_angle_distribution(),
                             step_m=1000.0, lower_altitude_m=80e3,
                             upper_altitude_m=1000e3, max_steps=1_000_000)
    optical_depth_to_collision = -log(rand(rng))
    accumulated_depth = 0.0
    for step in 1:max_steps
        position, velocity = particle.position_m, particle.velocity_m_s
        radius = _norm(position)
        altitude = radius - MARS_RADIUS_M
        if altitude <= lower_altitude_m
            particle.alive = false
            return (particle=particle, reason=:lower_boundary, steps=step)
        elseif altitude >= upper_altitude_m
            particle.alive = false
            specific_energy = 0.5 * _dot(velocity, velocity) - MARS_MU_M3_S2 / radius
            reason = specific_energy >= 0 ? :escaped : :upper_boundary_bound
            return (particle=particle, reason=reason, steps=step)
        end
        local_state = interpolate_profile(atmosphere, altitude)
        energy_eV = kinetic_energy_eV(velocity)
        coefficient = collision_coefficient(
            targets, local_state.density_m3, energy_eV,
        )
        speed = _norm(velocity)
        dt = step_m / max(speed, 1.0)
        acceleration = _scale(-MARS_MU_M3_S2 / radius^3, position)
        velocity = _add(velocity, _scale(dt, acceleration))
        direction = _scale(1 / _norm(velocity), velocity)
        particle.position_m = _add(position, _scale(step_m, direction))
        particle.velocity_m_s = velocity
        accumulated_depth += coefficient * step_m
        if coefficient > 0 && accumulated_depth >= optical_depth_to_collision
            target = choose_collision_target(
                rng, targets, local_state.density_m3, energy_eV,
            )
            theta = sample_scattering_angle(
                rng, scattering_distribution, O_MASS_KG, target.mass_kg,
            )
            phi = sample_azimuth(rng)
            particle.velocity_m_s, _ = elastic_collision_lab(
                particle.velocity_m_s,
                O_MASS_KG, target.mass_kg, theta, phi,
            )
            particle.collisions += 1
            accumulated_depth = 0.0
            optical_depth_to_collision = -log(rand(rng))
        end
    end
    particle.alive = false
    (particle=particle, reason=:maximum_steps, steps=max_steps)
end
