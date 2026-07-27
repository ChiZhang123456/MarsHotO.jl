"""Reproduce hot O nascent energy maps from an MGITM vertical profile.

This script creates:
1. A Lillis Figure 1 style conditional energy probability map.
2. An Ali Rahmati Figure 2.4 style spectral production rate map.

The nascent energy broadening is calculated from Maxwellian electron and O2+
velocities followed by two-body reaction kinematics. Collisions with the
neutral atmosphere are not part of these source maps.
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from plot_mgitm_hot_o_profiles import (
    INPUT_FILE,
    calculate_derived_profiles,
    load_nearest_subsolar_profile,
)


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "figures"
OUTPUT_PREFIX = (
    OUTPUT_DIR / "mgitm_ls000_f070_hot_o_nascent_energy_with_vibration"
)
SOURCE_DATA_FILE = (
    OUTPUT_DIR / "mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_map.csv"
)

RANDOM_SEED = 20260727
REACTIONS_PER_ALTITUDE = 200_000
ENERGY_EDGES_EV = np.linspace(0.0, 7.0, 281)

ELECTRON_MASS_KG = 9.1093837139e-31
ATOMIC_MASS_UNIT_KG = 1.66053906892e-27
O_MASS_KG = 15.999 * ATOMIC_MASS_UNIT_KG
O2_ION_MASS_KG = 2.0 * O_MASS_KG
BOLTZMANN_JK = 1.380649e-23
ELECTRON_VOLT_J = 1.602176634e-19

# User-specified reaction branches. Do not alter without explicit approval.
BRANCH_RELEASE_ENERGY_EV = np.array([6.99, 5.02, 3.06, 0.83])
BRANCH_PROBABILITY = np.array([0.265, 0.473, 0.204, 0.058])

# Mars exobase O2+ vibrational distribution from Fox and Hac (1997), Table 2.
# The published rounded fractions are normalized before sampling.
VIBRATIONAL_QUANTUM_NUMBER = np.arange(9)
VIBRATIONAL_FRACTION_RAW = np.array(
    [0.800, 0.074, 0.043, 0.035, 0.025, 0.015, 0.0047, 0.00027, 0.00021]
)
VIBRATIONAL_PROBABILITY = (
    VIBRATIONAL_FRACTION_RAW / VIBRATIONAL_FRACTION_RAW.sum()
)
VIBRATIONAL_QUANTUM_ENERGY_EV = 0.23


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def sample_maxwellian_velocity(
    rng: np.random.Generator,
    temperature_k: float,
    mass_kg: float,
    sample_count: int,
) -> np.ndarray:
    component_sigma_ms = np.sqrt(BOLTZMANN_JK * temperature_k / mass_kg)
    return rng.normal(0.0, component_sigma_ms, size=(sample_count, 3))


def sample_nascent_o_energies(
    rng: np.random.Generator,
    electron_temperature_k: float,
    ion_temperature_k: float,
    reaction_count: int,
) -> np.ndarray:
    electron_velocity = sample_maxwellian_velocity(
        rng,
        electron_temperature_k,
        ELECTRON_MASS_KG,
        reaction_count,
    )
    ion_velocity = sample_maxwellian_velocity(
        rng,
        ion_temperature_k,
        O2_ION_MASS_KG,
        reaction_count,
    )

    total_reactant_mass = ELECTRON_MASS_KG + O2_ION_MASS_KG
    center_of_mass_velocity = (
        ELECTRON_MASS_KG * electron_velocity + O2_ION_MASS_KG * ion_velocity
    ) / total_reactant_mass

    relative_velocity = electron_velocity - ion_velocity
    reduced_mass = (
        ELECTRON_MASS_KG * O2_ION_MASS_KG / total_reactant_mass
    )
    relative_energy_ev = (
        0.5
        * reduced_mass
        * np.sum(relative_velocity * relative_velocity, axis=1)
        / ELECTRON_VOLT_J
    )

    branch_index = rng.choice(
        len(BRANCH_PROBABILITY),
        size=reaction_count,
        p=BRANCH_PROBABILITY,
    )
    vibrational_index = rng.choice(
        len(VIBRATIONAL_PROBABILITY),
        size=reaction_count,
        p=VIBRATIONAL_PROBABILITY,
    )
    available_energy_ev = (
        BRANCH_RELEASE_ENERGY_EV[branch_index]
        + VIBRATIONAL_QUANTUM_ENERGY_EV
        * VIBRATIONAL_QUANTUM_NUMBER[vibrational_index]
        + relative_energy_ev
    )

    product_energy_com_ev = 0.5 * available_energy_ev
    product_speed_com_ms = np.sqrt(
        2.0 * product_energy_com_ev * ELECTRON_VOLT_J / O_MASS_KG
    )

    direction = rng.normal(size=(reaction_count, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    product_velocity_com = product_speed_com_ms[:, None] * direction

    product_velocity_1 = center_of_mass_velocity + product_velocity_com
    product_velocity_2 = center_of_mass_velocity - product_velocity_com

    product_energy_1_ev = (
        0.5
        * O_MASS_KG
        * np.sum(product_velocity_1 * product_velocity_1, axis=1)
        / ELECTRON_VOLT_J
    )
    product_energy_2_ev = (
        0.5
        * O_MASS_KG
        * np.sum(product_velocity_2 * product_velocity_2, axis=1)
        / ELECTRON_VOLT_J
    )
    return np.concatenate((product_energy_1_ev, product_energy_2_ev))


def calculate_energy_maps(
    profile: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(RANDOM_SEED)
    energy_centers_ev = 0.5 * (ENERGY_EDGES_EV[:-1] + ENERGY_EDGES_EV[1:])
    energy_bin_width_ev = np.diff(ENERGY_EDGES_EV)
    probability_density_ev1 = np.empty(
        (len(profile), len(energy_centers_ev)),
        dtype=float,
    )

    for altitude_index, row in profile.iterrows():
        energies_ev = sample_nascent_o_energies(
            rng,
            electron_temperature_k=float(row["Te_K"]),
            ion_temperature_k=float(row["Ti_K"]),
            reaction_count=REACTIONS_PER_ALTITUDE,
        )
        counts, _ = np.histogram(energies_ev, bins=ENERGY_EDGES_EV)
        probability_density_ev1[altitude_index, :] = (
            counts / counts.sum() / energy_bin_width_ev
        )

    spectral_production_cm3s_ev1 = (
        profile["Q_hotO_cm3s"].to_numpy()[:, None] * probability_density_ev1
    )
    return (
        energy_centers_ev,
        probability_density_ev1,
        spectral_production_cm3s_ev1,
    )


def save_source_data(
    profile: pd.DataFrame,
    energy_centers_ev: np.ndarray,
    probability_density_ev1: np.ndarray,
    spectral_production_cm3s_ev1: np.ndarray,
) -> None:
    altitude_grid_km, energy_grid_ev = np.meshgrid(
        profile["altitude_km"].to_numpy(),
        energy_centers_ev,
        indexing="ij",
    )
    source_data = pd.DataFrame(
        {
            "altitude_km": altitude_grid_km.ravel(),
            "energy_eV": energy_grid_ev.ravel(),
            "probability_density_eV-1": probability_density_ev1.ravel(),
            "spectral_production_cm-3_s-1_eV-1": (
                spectral_production_cm3s_ev1.ravel()
            ),
        }
    )
    source_data.to_csv(SOURCE_DATA_FILE, index=False)


def escape_energy_ev(altitude_km: np.ndarray) -> np.ndarray:
    gravitational_constant = 6.67430e-11
    mars_mass_kg = 6.4171e23
    mars_radius_m = 3388.25e3
    radial_distance_m = mars_radius_m + altitude_km * 1.0e3
    return (
        gravitational_constant * mars_mass_kg * O_MASS_KG
        / radial_distance_m
        / ELECTRON_VOLT_J
    )


def make_probability_figure(
    altitude_km: np.ndarray,
    energy_centers_ev: np.ndarray,
    probability_density_ev1: np.ndarray,
) -> plt.Figure:
    configure_matplotlib()
    fig, axis = plt.subplots(figsize=(4.5, 4.0), constrained_layout=True)
    positive_values = probability_density_ev1[probability_density_ev1 > 0.0]
    image = axis.pcolormesh(
        energy_centers_ev,
        altitude_km,
        probability_density_ev1,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=float(np.percentile(positive_values, 99.5)),
        rasterized=True,
    )
    axis.plot(
        escape_energy_ev(altitude_km),
        altitude_km,
        color="white",
        lw=1.2,
        ls="--",
        label="Escape energy",
    )
    axis.set(
        xlabel="Nascent O energy (eV)",
        ylabel="Altitude (km)",
        xlim=(0.0, 7.0),
        ylim=(float(altitude_km.min()), float(altitude_km.max())),
        title=(
            r"Nascent hot O probability with O$_2^+$ vibration, "
            r"$L_s=0^\circ$, F070"
        ),
    )
    axis.tick_params(axis="y", labelleft=True)
    legend = axis.legend(loc="upper right", frameon=False)
    for text in legend.get_texts():
        text.set_color("white")
    colorbar = fig.colorbar(image, ax=axis, pad=0.03)
    colorbar.set_label(r"Probability density (eV$^{-1}$)")
    return fig


def make_production_figure(
    altitude_km: np.ndarray,
    energy_centers_ev: np.ndarray,
    spectral_production_cm3s_ev1: np.ndarray,
) -> plt.Figure:
    configure_matplotlib()
    fig, axis = plt.subplots(figsize=(4.5, 4.0), constrained_layout=True)
    log_production = np.log10(
        np.maximum(spectral_production_cm3s_ev1, 1.0e-6)
    )
    image = axis.pcolormesh(
        energy_centers_ev,
        altitude_km,
        log_production,
        shading="auto",
        cmap="turbo",
        vmin=-4.0,
        vmax=4.5,
        rasterized=True,
    )
    axis.set(
        xlabel="Nascent O energy (eV)",
        ylabel="Altitude (km)",
        xlim=(0.0, 7.0),
        ylim=(float(altitude_km.min()), float(altitude_km.max())),
        title=(
            r"Hot O production with O$_2^+$ vibration, "
            r"$L_s=0^\circ$, F070"
        ),
    )
    colorbar = fig.colorbar(image, ax=axis, pad=0.03)
    colorbar.set_label(
        r"$\log_{10}\,[Q(E,z)\;(\mathrm{cm^{-3}\,s^{-1}\,eV^{-1}})]$"
    )
    return fig


def save_figure(fig: plt.Figure, suffix: str) -> None:
    stem = Path(f"{OUTPUT_PREFIX}_{suffix}")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    profile, metadata = load_nearest_subsolar_profile(INPUT_FILE)
    profile = calculate_derived_profiles(profile)
    (
        energy_centers_ev,
        probability_density_ev1,
        spectral_production_cm3s_ev1,
    ) = calculate_energy_maps(profile)
    save_source_data(
        profile,
        energy_centers_ev,
        probability_density_ev1,
        spectral_production_cm3s_ev1,
    )

    probability_figure = make_probability_figure(
        profile["altitude_km"].to_numpy(),
        energy_centers_ev,
        probability_density_ev1,
    )
    save_figure(probability_figure, "probability")

    production_figure = make_production_figure(
        profile["altitude_km"].to_numpy(),
        energy_centers_ev,
        spectral_production_cm3s_ev1,
    )
    save_figure(production_figure, "production")

    integrated_production = np.sum(
        spectral_production_cm3s_ev1 * np.diff(ENERGY_EDGES_EV)[None, :],
        axis=1,
    )
    relative_error = np.max(
        np.abs(integrated_production - profile["Q_hotO_cm3s"].to_numpy())
        / profile["Q_hotO_cm3s"].to_numpy()
    )
    peak_index = np.unravel_index(
        np.argmax(spectral_production_cm3s_ev1),
        spectral_production_cm3s_ev1.shape,
    )

    print(f"Input: {INPUT_FILE}")
    print(
        "Selected columns: "
        f"lon={metadata['selected_longitudes_deg']}, "
        f"lat={metadata['selected_latitudes_deg']}, "
        f"SZA={metadata['minimum_grid_sza_deg']:.4f} deg"
    )
    print(f"Random seed: {RANDOM_SEED}")
    print(f"Reactions per altitude: {REACTIONS_PER_ALTITUDE}")
    print(
        "Normalized vibrational probabilities: "
        f"{VIBRATIONAL_PROBABILITY.tolist()}"
    )
    print(
        "Mean vibrational energy added per reaction: "
        f"{np.sum(VIBRATIONAL_PROBABILITY * VIBRATIONAL_QUANTUM_NUMBER) * VIBRATIONAL_QUANTUM_ENERGY_EV:.6f} eV"
    )
    print(f"Maximum spectral integration relative error: {relative_error:.3e}")
    print(
        "Maximum spectral production: "
        f"{spectral_production_cm3s_ev1[peak_index]:.6e} cm^-3 s^-1 eV^-1 "
        f"at altitude={profile.loc[peak_index[0], 'altitude_km']:.2f} km, "
        f"energy={energy_centers_ev[peak_index[1]]:.3f} eV"
    )
    print(f"Source data: {SOURCE_DATA_FILE}")


if __name__ == "__main__":
    main()
