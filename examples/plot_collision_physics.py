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
PROJECTILE_MASS_AMU = 16.0
RAHMATI_BETA = -1.85

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


def inverse_scattering_angle(random_number: np.ndarray) -> np.ndarray:
    """Rahmati inverse CDF with theta_min = 0."""
    exponent = RAHMATI_BETA + 2.0
    theta_rad = 2.0 * np.arcsin(random_number ** (1.0 / exponent))
    return np.rad2deg(theta_rad)


def fractional_energy_loss(
    theta_com_rad: np.ndarray,
    target_mass_amu: float,
) -> np.ndarray:
    return (
        2.0
        * PROJECTILE_MASS_AMU
        * target_mass_amu
        / (PROJECTILE_MASS_AMU + target_mass_amu) ** 2
        * (1.0 - np.cos(theta_com_rad))
    )


def total_cross_section(
    energy_ev: np.ndarray,
    sigma_3ev_cm2: float,
) -> np.ndarray:
    return sigma_3ev_cm2 * (energy_ev / 3.0) ** -0.2


def make_figure() -> plt.Figure:
    configure_matplotlib()
    theta_deg = np.linspace(0.0, 180.0, 2000)
    theta_rad = np.deg2rad(theta_deg)

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

    random_number = np.linspace(0.0, 1.0, 4000)
    sampled_angle_deg = inverse_scattering_angle(random_number)
    axes[1].plot(
        random_number,
        sampled_angle_deg,
        color="#2878B5",
        linewidth=1.5,
        label=r"$\beta=-1.85,\ \theta_{\min}=0^\circ$",
    )
    axes[1].set(
        xlim=(0.0, 1.0),
        ylim=(0.0, 180.0),
        xlabel="Uniform random number, R",
        ylabel=r"COM scattering angle, $\theta$ (deg)",
        title="Rahmati inverse-CDF sampling",
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
        species_theta_deg = np.linspace(0.0, 180.0, 2000)
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
        xlabel=r"COM scattering angle, $\theta$ (deg)",
        ylabel=r"Fractional energy loss, $\Delta E/E$",
        title="Elastic energy transfer",
    )
    axes[2].legend(ncol=2, loc="upper left")

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
        "Hot O collision physics with Rahmati COM angle sampling",
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
