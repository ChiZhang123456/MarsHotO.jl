"""Plot the scattering and energy-loss physics used by MarsHotO."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "examples" / "figures" /
    "hot_o_collision_cross_sections_and_scattering.png"
)
SCATTERING_FILE = (
    ROOT / "data" / "cross_sections" /
    "scattering_angle_distribution.txt"
)

PROJECTILE_MASS_AMU = 16.0
RANDOM_SEED = 20260727
SAMPLE_COUNT = 1_000_000

TARGETS = {
    "O": (16.0, 6.4e-15),
    "CO": (28.0, 1.8e-14),
    r"N$_2$": (28.0, 1.8e-14),
    r"O$_2$": (32.0, 1.8e-14),
    r"CO$_2$": (44.0, 2.0e-14),
}

COLORS = {
    "O": "#2878B5",
    "CO": "#9C755F",
    r"N$_2$": "#59A14F",
    r"O$_2$": "#F28E2B",
    r"CO$_2$": "#E15759",
}


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "legend.frameon": False,
        }
    )


def load_scattering_distribution() -> tuple[np.ndarray, np.ndarray]:
    table = np.loadtxt(
        SCATTERING_FILE, comments="#", skiprows=8, dtype=float
    )
    order = np.argsort(table[:, 0])
    return table[order, 0], table[order, 1]


def fractional_energy_loss(
    theta_rad: np.ndarray,
    target_mass_amu: float,
) -> np.ndarray:
    mass_ratio = PROJECTILE_MASS_AMU / target_mass_amu
    discriminant = np.maximum(
        1.0 - (mass_ratio * np.sin(theta_rad)) ** 2, 0.0
    )
    speed_ratio = np.maximum(
        (
            mass_ratio * np.cos(theta_rad)
            + np.sqrt(discriminant)
        )
        / (1.0 + mass_ratio),
        0.0,
    )
    return 1.0 - speed_ratio**2


def total_cross_section(
    energy_ev: np.ndarray,
    sigma_3ev_cm2: float,
) -> np.ndarray:
    return sigma_3ev_cm2 * (energy_ev / 3.0) ** -0.2


def make_figure() -> plt.Figure:
    configure_matplotlib()
    random_grid, angle_grid_deg = load_scattering_distribution()
    theta_deg = np.linspace(0.12, 180.0, 2000)
    theta_rad = np.deg2rad(theta_deg)

    rng = np.random.default_rng(RANDOM_SEED)
    samples_deg = np.interp(
        rng.random(SAMPLE_COUNT), random_grid, angle_grid_deg
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(9.0, 3.2),
        constrained_layout=True,
    )

    energy_ev = np.geomspace(0.01, 100.0, 600)
    cross_section_groups = {
        "O": TARGETS["O"],
        r"CO, N$_2$, O$_2$": TARGETS["CO"],
        r"CO$_2$": TARGETS[r"CO$_2$"],
    }
    cross_section_colors = {
        "O": COLORS["O"],
        r"CO, N$_2$, O$_2$": COLORS[r"O$_2$"],
        r"CO$_2$": COLORS[r"CO$_2$"],
    }
    for species, (_, sigma_3ev_cm2) in cross_section_groups.items():
        axes[0].plot(
            energy_ev,
            total_cross_section(energy_ev, sigma_3ev_cm2),
            color=cross_section_colors[species],
            linewidth=1.4,
            label=species,
        )
    axes[0].axvline(
        3.0, color="0.35", linestyle="--", linewidth=0.9, label="3 eV"
    )
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel="Relative collision energy (eV)",
        ylabel=r"Total cross section (cm$^2$)",
        title="Energy-dependent cross sections",
    )
    axes[0].legend(loc="best")

    bin_edges_deg = np.linspace(0.0, 180.0, 91)
    sampled_count, _ = np.histogram(samples_deg, bins=bin_edges_deg)
    bin_width_deg = np.diff(bin_edges_deg)
    sampled_density = sampled_count / (SAMPLE_COUNT * bin_width_deg)
    expected_fraction = np.diff(
        np.interp(bin_edges_deg, angle_grid_deg, random_grid)
    )
    expected_density = expected_fraction / bin_width_deg
    bin_center = 0.5 * (bin_edges_deg[:-1] + bin_edges_deg[1:])
    axes[1].stairs(
        sampled_density,
        bin_edges_deg,
        fill=True,
        color="#E07B39",
        alpha=0.55,
        label=f"Monte Carlo, N = {SAMPLE_COUNT:,}",
    )
    axes[1].plot(
        bin_center,
        expected_density,
        color="#222222",
        linewidth=1.3,
        label="Lookup-table probability",
    )
    axes[1].set(
        xlim=(0.0, 180.0),
        yscale="log",
        ylim=(1.0e-6, 1.0),
        xlabel=r"LAB scattering angle, $\theta$ (deg)",
        ylabel=r"Probability density (deg$^{-1}$)",
        title="Kallio and Barabash angle distribution",
    )
    axes[1].legend(loc="best")

    energy_loss_groups = {
        "O": TARGETS["O"][0],
        r"CO, N$_2$": TARGETS["CO"][0],
        r"O$_2$": TARGETS[r"O$_2$"][0],
        r"CO$_2$": TARGETS[r"CO$_2$"][0],
    }
    energy_loss_colors = {
        "O": COLORS["O"],
        r"CO, N$_2$": COLORS["CO"],
        r"O$_2$": COLORS[r"O$_2$"],
        r"CO$_2$": COLORS[r"CO$_2$"],
    }
    for species, mass_amu in energy_loss_groups.items():
        theta_max_deg = (
            180.0 if PROJECTILE_MASS_AMU < mass_amu else
            np.rad2deg(np.arcsin(mass_amu / PROJECTILE_MASS_AMU))
        )
        species_theta_deg = np.linspace(0.12, theta_max_deg, 2000)
        axes[2].plot(
            species_theta_deg,
            fractional_energy_loss(
                np.deg2rad(species_theta_deg), mass_amu
            ),
            color=energy_loss_colors[species],
            linewidth=1.3,
            label=species,
        )
    axes[2].set(
        xlim=(0.0, 180.0),
        ylim=(0.0, 1.02),
        xlabel=r"LAB scattering angle, $\theta$ (deg)",
        ylabel=r"Fractional energy loss, $\Delta E/E$",
        title="Elastic energy transfer",
    )
    axes[2].legend(ncol=2, loc="upper left")
    axes[2].axvline(
        90.0, color=COLORS["O"], linestyle=":", linewidth=0.9
    )
    axes[2].text(
        92.0, 0.04, "O + O LAB limit",
        color=COLORS["O"], fontsize=7, va="bottom", rotation=90,
    )

    for label, axis in zip(("a", "b", "c"), axes):
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
        axis.grid(True, which="both", color="0.91", linewidth=0.5)

    figure.suptitle(
        "Hot O collision physics with Kallio and Barabash angle sampling",
        fontsize=10,
    )
    return figure


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure = make_figure()
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
