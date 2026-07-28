"""Reproduce hot O nascent energy maps from an MGITM vertical profile.

This script creates:
1. A Lillis Figure 1 style conditional energy probability map.
2. A Rahmati Figure 2.4 style spectral production rate map.

The electron and O2+ bulk velocities are zero. Their velocities are sampled
from normalized three-dimensional Maxwellian distributions. Isotropy follows
from the equal Gaussian variance of all three Cartesian components. Collisions
with the neutral atmosphere are not part of these source maps.
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
OUTPUT_DIR = ROOT / "examples" / "figures"
OUTPUT_PREFIX = (
    OUTPUT_DIR / "mgitm_ls000_f070_hot_o_nascent_energy_with_vibration"
)
OUTPUT_FIGURE = Path(f"{OUTPUT_PREFIX}_energy_maps.png")
SOURCE_DATA_FILE = (
    OUTPUT_DIR / "mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_map.csv"
)

RANDOM_SEED = 20260727
REACTIONS_PER_ALTITUDE = 1_000_000
REACTIONS_PER_BATCH = 100_000
ENERGY_EDGES_EV = np.linspace(0.0, 7.0, 281)
PLOT_INTERPOLATION = "bilinear"

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
    bulk_velocity_ms: np.ndarray | None = None,
) -> np.ndarray:
    """Sample a normalized three-dimensional Maxwellian velocity."""
    if temperature_k <= 0.0:
        raise ValueError("Temperature must be positive")
    if mass_kg <= 0.0:
        raise ValueError("Mass must be positive")
    if bulk_velocity_ms is None:
        bulk_velocity_ms = np.zeros(3)
    bulk_velocity_ms = np.asarray(bulk_velocity_ms, dtype=float)
    if bulk_velocity_ms.shape != (3,):
        raise ValueError("Bulk velocity must contain three components")

    component_sigma_ms = np.sqrt(BOLTZMANN_JK * temperature_k / mass_kg)
    return rng.normal(
        loc=bulk_velocity_ms,
        scale=component_sigma_ms,
        size=(sample_count, 3),
    )


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

    cosine_polar = rng.uniform(-1.0, 1.0, size=reaction_count)
    azimuth = rng.uniform(0.0, 2.0 * np.pi, size=reaction_count)
    sine_polar = np.sqrt(np.maximum(1.0 - cosine_polar**2, 0.0))
    direction = np.column_stack(
        (
            sine_polar * np.cos(azimuth),
            sine_polar * np.sin(azimuth),
            cosine_polar,
        )
    )
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
        counts = np.zeros(len(energy_centers_ev), dtype=np.int64)
        reactions_remaining = REACTIONS_PER_ALTITUDE
        while reactions_remaining > 0:
            batch_size = min(REACTIONS_PER_BATCH, reactions_remaining)
            energies_ev = sample_nascent_o_energies(
                rng,
                electron_temperature_k=float(row["Te_K"]),
                ion_temperature_k=float(row["Ti_K"]),
                reaction_count=batch_size,
            )
            batch_counts, _ = np.histogram(energies_ev, bins=ENERGY_EDGES_EV)
            counts += batch_counts
            reactions_remaining -= batch_size
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


def make_combined_figure(
    altitude_km: np.ndarray,
    energy_centers_ev: np.ndarray,
    probability_density_ev1: np.ndarray,
    spectral_production_cm3s_ev1: np.ndarray,
) -> plt.Figure:
    configure_matplotlib()
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.0, 4.0),
        sharey=True,
        constrained_layout=True,
    )
    positive_values = probability_density_ev1[probability_density_ev1 > 0.0]
    probability_image = axes[0].imshow(
        probability_density_ev1,
        origin="lower",
        aspect="auto",
        extent=(
            float(energy_centers_ev[0]),
            float(energy_centers_ev[-1]),
            float(altitude_km[0]),
            float(altitude_km[-1]),
        ),
        interpolation=PLOT_INTERPOLATION,
        cmap="viridis",
        vmin=0.0,
        vmax=float(np.percentile(positive_values, 99.5)),
        rasterized=True,
    )
    axes[0].plot(
        escape_energy_ev(altitude_km),
        altitude_km,
        color="white",
        lw=1.2,
        ls="--",
        label="Escape energy",
    )
    axes[0].set(
        xlabel="Nascent O energy (eV)",
        ylabel="Altitude (km)",
        xlim=(0.0, 5.0),
        ylim=(float(altitude_km.min()), float(altitude_km.max())),
        title="Conditional energy probability",
    )
    axes[0].tick_params(axis="y", labelleft=True)
    legend = axes[0].legend(loc="upper right", frameon=False)
    for text in legend.get_texts():
        text.set_color("white")
    probability_colorbar = fig.colorbar(
        probability_image,
        ax=axes[0],
        pad=0.03,
    )
    probability_colorbar.set_label(r"Probability density (eV$^{-1}$)")

    log_production = np.log10(
        np.maximum(spectral_production_cm3s_ev1, 1.0e-6)
    )
    production_image = axes[1].imshow(
        log_production,
        origin="lower",
        aspect="auto",
        extent=(
            float(energy_centers_ev[0]),
            float(energy_centers_ev[-1]),
            float(altitude_km[0]),
            float(altitude_km[-1]),
        ),
        interpolation=PLOT_INTERPOLATION,
        cmap="turbo",
        vmin=-4.0,
        vmax=4.5,
        rasterized=True,
    )
    axes[1].set(
        xlabel="Nascent O energy (eV)",
        xlim=(0.0, 5.0),
        ylim=(float(altitude_km.min()), float(altitude_km.max())),
        title="Spectral production rate",
    )
    production_colorbar = fig.colorbar(
        production_image,
        ax=axes[1],
        pad=0.03,
    )
    production_colorbar.set_label(
        r"$\log_{10}\,[Q(E,z)\;(\mathrm{cm^{-3}\,s^{-1}\,eV^{-1}})]$"
    )

    for label, axis in zip(("a", "b"), axes):
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
        )
    fig.suptitle(
        r"Nascent hot O with O$_2^+$ vibration, $L_s=0^\circ$, F070",
        fontsize=10,
    )
    return fig


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

    combined_figure = make_combined_figure(
        profile["altitude_km"].to_numpy(),
        energy_centers_ev,
        probability_density_ev1,
        spectral_production_cm3s_ev1,
    )
    combined_figure.savefig(OUTPUT_FIGURE, dpi=400, bbox_inches="tight")
    plt.close(combined_figure)

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
    print(f"Reactions per batch: {REACTIONS_PER_BATCH}")
    print(f"Plot interpolation: {PLOT_INTERPOLATION}")
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
