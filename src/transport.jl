"""
Propagation of one hot O particle through the spherical Mars atmosphere.

This file supports direct single-particle calculations. It advances a
particle under gravity, evaluates collision probability, samples collision
physics, and returns any recoil O secondary for optional further tracking.
"""
@inline kinetic_energy_eV(velocity, mass=O_MASS_KG) =
    0.5mass * _dot(velocity, velocity) / EV_J

function _advance_gravity(position, velocity, ds)
    radius = _norm(position)
    speed = _norm(velocity)
    dt = ds / max(speed, 1.0)
    acceleration0 = _scale(-MARS_MU_M3_S2 / radius^3, position)
    position1 = _add(
        position,
        _add(_scale(dt, velocity), _scale(0.5dt^2, acceleration0)),
    )
    radius1 = _norm(position1)
    acceleration1 = _scale(-MARS_MU_M3_S2 / radius1^3, position1)
    velocity1 = _add(
        velocity, _scale(0.5dt, _add(acceleration0, acceleration1)),
    )
    position1, velocity1, dt
end

"""
Advance one known O particle by one Rahmati step.

The supplied RNG is the entire stochastic state of the trajectory. Every step
draws a collision variate. If a collision occurs, separate variates select the
neutral target, COM scattering angle, and azimuth. The returned record exposes
these values for exact trajectory diagnostics and reproducibility.
"""
function advance_hot_o_step!(
    rng, particle::HotOParticle, atmosphere::AtmosphereProfile, targets;
    step_m=nothing, minimum_secondary_energy_eV=0.01,
)
    position0 = particle.position_m
    velocity0 = particle.velocity_m_s
    altitude_m = _norm(position0) - MARS_RADIUS_M
    local_state = interpolate_profile(atmosphere, altitude_m)
    energy_eV = kinetic_energy_eV(velocity0)
    kappa = collision_coefficient(targets, local_state.density_m3, energy_eV)
    mfp = kappa > 0 ? inv(kappa) : Inf
    ds = isnothing(step_m) ? rahmati_step_length(mfp) : Float64(step_m)
    ds > 0 || throw(DomainError(ds, "Step length must be positive"))
    position1, velocity1, dt = _advance_gravity(position0, velocity0, ds)
    particle.position_m = position1
    particle.velocity_m_s = velocity1

    collision_u = rand(rng)
    target_u = nothing
    scattering_u = nothing
    azimuth_u = nothing
    target = nothing
    theta_com = nothing
    secondary = nothing
    if kappa > 0 && collision_u < min(ds * kappa, 1.0)
        target_u = rand(rng)
        target = choose_collision_target(
            targets, local_state.density_m3, energy_eV, target_u,
        )
        scattering_u = rand(rng)
        azimuth_u = rand(rng)
        theta_com = sample_scattering_angle_from_uniform(scattering_u)
        projectile_after, target_after = elastic_collision(
            velocity1, (0.0, 0.0, 0.0), O_MASS_KG, target.mass_kg,
            theta_com, 2pi * azimuth_u,
        )
        particle.velocity_m_s = projectile_after
        particle.collisions += 1
        if target.species == :O &&
           kinetic_energy_eV(target_after) > minimum_secondary_energy_eV
            secondary = HotOParticle(
                position1, target_after, particle.weight_s1, true, 0,
            )
        end
    end
    HotOTransportStep(
        position0, velocity0, position1, velocity1, particle.velocity_m_s, dt, ds,
        kappa, collision_u, target_u, scattering_u, azimuth_u, target,
        theta_com, secondary,
    )
end

"""
Transport one hot O particle in a spherically symmetric MGITM atmosphere.
Collisions are sampled by optical depth. Gravity is updated at every spatial
step. Returns a named tuple with the final state and stopping reason.
"""
function transport_particle!(rng, particle::HotOParticle,
                             atmosphere::AtmosphereProfile, targets;
                             step_m=1000.0, lower_altitude_m=80e3,
                             upper_altitude_m=1000e3, max_steps=1_000_000,
                             record_steps=false)
    secondaries = HotOParticle[]
    step_history = HotOTransportStep[]
    for step in 1:max_steps
        position, velocity = particle.position_m, particle.velocity_m_s
        radius = _norm(position)
        altitude = radius - MARS_RADIUS_M
        if altitude <= lower_altitude_m
            particle.alive = false
            return (particle=particle, reason=:lower_boundary, steps=step,
                    secondaries=secondaries, step_history=step_history)
        elseif altitude >= upper_altitude_m
            particle.alive = false
            specific_energy = 0.5 * _dot(velocity, velocity) - MARS_MU_M3_S2 / radius
            reason = specific_energy >= 0 ? :escaped : :upper_boundary_bound
            return (particle=particle, reason=reason, steps=step,
                    secondaries=secondaries, step_history=step_history)
        end
        step_result = advance_hot_o_step!(
            rng, particle, atmosphere, targets; step_m=step_m,
        )
        record_steps && push!(step_history, step_result)
        !isnothing(step_result.secondary) &&
            push!(secondaries, step_result.secondary)
    end
    particle.alive = false
    (particle=particle, reason=:maximum_steps, steps=max_steps,
     secondaries=secondaries, step_history=step_history)
end
