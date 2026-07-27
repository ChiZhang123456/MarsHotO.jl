"""Plot the scattering and energy-loss physics used by MarsHotO."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "examples" / "figures" / "hot_o_collision_physics.png"

ALPHA_CM2_SR = 0.36e-16
BETA = -1.85
THETA_MIN_DEG = 10.0
PROJECTILE_MASS_AMU = 16.0
RANDOM_SEED = 20260727
SAMPLE_COUNT = 1_000_000

TARGETS = {
    "O": (16.0, 6.4e-15),
    "CO": (28.0, 1.8e-14),
    r"N$_2$": (28.0, 1.8e-14),
    r"O$_2$": (32.0, 1.8e-14),
    "Ar": (40.0, 1.2e-14),
    r"CO$_2$": (44.0, 2.0e-14),
}

COLORS = {
    "O": "#2878B5",
    "CO": "#9C755F",
    r"N$_2$": "#59A14F",
    r"O$_2$": "#F28E2B",
    "Ar": "#B07AA1",
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


def differential_cross_section(theta_rad: np.ndarray) -> np.ndarray:
    return ALPHA_CM2_SR * np.sin(theta_rad / 2.0) ** BETA


def angle_pdf(theta_rad: np.ndarray, theta_min_rad: float) -> np.ndarray:
    exponent = BETA + 2.0
    denominator = 1.0 - np.sin(theta_min_rad / 2.0) ** exponent
    probability = (
        exponent
        * np.sin(theta_rad / 2.0) ** (BETA + 1.0)
        * np.cos(theta_rad / 2.0)
        / (2.0 * denominator)
    )
    return np.where(theta_rad >= theta_min_rad, probability, 0.0)


def sample_scattering_angle(
    rng: np.random.Generator,
    count: int,
    theta_min_rad: float,
) -> np.ndarray:
    exponent = BETA + 2.0
    lower = np.sin(theta_min_rad / 2.0) ** exponent
    random_number = rng.random(count)
    return 2.0 * np.arcsin(
        (lower + random_number * (1.0 - lower)) ** (1.0 / exponent)
    )


def fractional_energy_loss(
    theta_rad: np.ndarray,
    target_mass_amu: float,
) -> np.ndarray:
    return (
        2.0
        * PROJECTILE_MASS_AMU
        * target_mass_amu
        / (PROJECTILE_MASS_AMU + target_mass_amu) ** 2
        * (1.0 - np.cos(theta_rad))
    )


def total_cross_section(
    energy_ev: np.ndarray,
    sigma_3ev_cm2: float,
) -> np.ndarray:
    return sigma_3ev_cm2 * (energy_ev / 3.0) ** -0.2


def make_figure() -> plt.Figure:
    configure_matplotlib()
    theta_min_rad = np.deg2rad(THETA_MIN_DEG)
    theta_deg = np.geomspace(0.1, 180.0, 2000)
    theta_rad = np.deg2rad(theta_deg)
    theta_cut_deg = np.geomspace(THETA_MIN_DEG, 179.9, 1200)
    theta_cut_rad = np.deg2rad(theta_cut_deg)

    rng = np.random.default_rng(RANDOM_SEED)
    samples_deg = np.rad2deg(
        sample_scattering_angle(rng, SAMPLE_COUNT, theta_min_rad)
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(7.2, 6.0),
        constrained_layout=True,
    )

    axes[0, 0].plot(
        theta_deg,
        differential_cross_section(theta_rad),
        color="#2878B5",
        linewidth=1.6,
    )
    axes[0, 0].axvline(
        THETA_MIN_DEG,
        color="0.25",
        linestyle="--",
        linewidth=1.0,
        label=rf"$\theta_{{\min}}={THETA_MIN_DEG:g}^\circ$",
    )
    axes[0, 0].set(
        xscale="log",
        yscale="log",
        xlabel=r"COM scattering angle, $\theta$ (deg)",
        ylabel=r"$d\sigma/d\Omega$ (cm$^2$ sr$^{-1}$)",
        title="Kharchenko O-O differential cross section",
    )
    axes[0, 0].legend(loc="best")

    bin_edges_deg = np.geomspace(THETA_MIN_DEG, 180.0, 61)
    axes[0, 1].hist(
        samples_deg,
        bins=bin_edges_deg,
        density=True,
        color="#E07B39",
        alpha=0.55,
        label=f"Monte Carlo, N = {SAMPLE_COUNT:,}",
    )
    axes[0, 1].plot(
        theta_cut_deg,
        angle_pdf(theta_cut_rad, theta_min_rad) * np.pi / 180.0,
        color="#222222",
        linewidth=1.3,
        label="Analytical PDF",
    )
    axes[0, 1].set(
        xscale="log",
        yscale="log",
        ylim=(1.0e-5, 1.0e-1),
        xlabel=r"COM scattering angle, $\theta$ (deg)",
        ylabel=r"Probability density (deg$^{-1}$)",
        title="Inverse-CDF scattering-angle sampling",
    )
    axes[0, 1].legend(loc="best")

    for species, (mass_amu, _) in TARGETS.items():
        axes[1, 0].plot(
            theta_deg,
            fractional_energy_loss(theta_rad, mass_amu),
            color=COLORS[species],
            linewidth=1.3,
            label=species,
        )
    axes[1, 0].axvline(
        THETA_MIN_DEG, color="0.35", linestyle="--", linewidth=0.9
    )
    axes[1, 0].set(
        xlim=(0.0, 180.0),
        ylim=(0.0, 1.02),
        xlabel=r"COM scattering angle, $\theta$ (deg)",
        ylabel=r"Fractional energy loss, $\Delta E/E$",
        title="Elastic energy transfer",
    )
    axes[1, 0].legend(ncol=2, loc="upper left")

    energy_ev = np.geomspace(0.01, 100.0, 600)
    for species, (_, sigma_3ev_cm2) in TARGETS.items():
        axes[1, 1].plot(
            energy_ev,
            total_cross_section(energy_ev, sigma_3ev_cm2),
            color=COLORS[species],
            linewidth=1.3,
            label=species,
        )
    axes[1, 1].axvline(
        3.0, color="0.35", linestyle="--", linewidth=0.9, label="3 eV"
    )
    axes[1, 1].set(
        xscale="log",
        yscale="log",
        xlabel="Relative collision energy (eV)",
        ylabel=r"Total cross section (cm$^2$)",
        title=r"$\sigma(E)=\sigma(3\,\mathrm{eV})(E/3\,\mathrm{eV})^{-0.2}$",
    )
    axes[1, 1].legend(ncol=2, loc="best")

    for label, axis in zip(("a", "b", "c", "d"), axes.flat):
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
        "Hot O collision physics used by MarsHotO",
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
