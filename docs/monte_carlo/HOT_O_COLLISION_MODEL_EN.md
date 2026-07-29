# MarsHotO Monte Carlo Physics Model Overview

MarsHotO consists of two coupled components, a hot O source model and a collisional transport model.

## 1. Hot O source model

The source model reads $n_e$, $n_{\mathrm{O_2^+}}$, $T_e$, and $T_i$ from MGITM. It calculates the O2+ dissociative recombination production rate, samples the four reaction branches and the O2+ vibrational states, and generates the initial position, energy, and isotropic velocity direction of each hot O atom in the Mars stationary frame.

See [Hot O altitude and nascent energy distributions](HOT_O_SOURCE_MODEL_EN.md) for details.

## 2. Collision frequency and target species

The neutral densities and energy dependent total cross sections jointly determine the mean free path, the collision probability during one numerical step, and the target species when a collision occurs.

See [Hot O collision cross sections with the neutral atmosphere](HOT_O_CROSS_SECTIONS_EN.md) for details.

## 3. Scattering angle and two body kinematics

The current model uses the analytical fit of Rahmati to the Kharchenko O and O differential cross section. The scattering angle is sampled in the center of mass frame. The default minimum angle is zero, so scattering below 10 degrees is retained. Postcollision velocities are calculated from momentum and energy conservation for an elastic collision with a stationary target, followed by transformation from the center of mass frame to the Mars stationary frame.

See [COM scattering angle and two body collisions](HOT_O_SCATTERING_TWO_BODY_EN.md) for details.

## 4. Complete Monte Carlo sequence

```text
Read the MGITM atmosphere and plasma profiles
    -> calculate the O2+ dissociative recombination source Q_hotO(z)
    -> sample source altitude, reaction branch, vibration, and nascent velocity
    -> calculate the mean free path from neutral densities and total cross sections
    -> advance the particle under Martian gravity
    -> sample whether a collision occurs during the step
    -> select the target species using n_s sigma_s
    -> sample the Rahmati COM scattering angle from the analytical inverse CDF
    -> sample a uniform azimuth
    -> scatter in COM and transform the postcollision velocities to the Mars frame
    -> continue tracking the primary hot O and any recoil secondary O
    -> estimate altitude and energy distributions from residence time
```

The source is represented through the spherical production rate

```math
4\pi r^2Q(z)\,dz.
```

In the stratified implementation, each source altitude contains a fixed number of primary macroparticles. The rate represented by one macroparticle is the source rate of that spherical shell divided by the number of primary particles at that altitude. A secondary O atom inherits the macroparticle rate of its parent. Extending the nearest subsolar profile as a spherically symmetric atmosphere is a current model approximation.

## 5. Main inputs and code

| Component | File |
|---|---|
| Dissociative recombination and vibration | `data/chemistry/o2plus_dissociative_recombination.toml` |
| Total collision cross sections | `data/cross_sections/rahmati_total_cross_sections.toml` |
| Rahmati analytical COM scattering model | `src/scattering.jl` |
| MGITM atmospheres | `MGITM/` |
| Initial particles | `src/source_particles.jl` |
| Scattering angle sampling | `src/scattering.jl` |
| Two body collisions | `src/collision_kinematics.jl` |
| Single particle transport | `src/transport.jl` |
| Monte Carlo ensembles | `src/ensembles.jl` |
| Complete example | `examples/run_hot_o_corona.jl` |

The overall procedure follows Rahmati for the adaptive step length, collision probability, COM scattering angle, energy loss relation, and stopping conditions. The file `data/cross_sections/scattering_angle_distribution.txt` is retained only as MarsASPEN reference data and is not used for runtime sampling.
