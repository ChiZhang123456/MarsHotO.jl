"""Plot normalized three-dimensional Maxwellian velocity sampling."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "examples" / "figures" / "thermal_energy_sampling_300K.png"
TEMPERATURE_K = 300.0
BOLTZMANN_J_K = 1.380649e-23
BOLTZMANN_EV_K = 8.617333262145e-5
ATOMIC_MASS_UNIT_KG = 1.66053906660e-27
O2_ION_MASS_KG = 32.0 * ATOMIC_MASS_UNIT_KG
RANDOM_SEED = 73
SAMPLE_COUNT = 500_000


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def main() -> None:
    configure_matplotlib()
    thermal_energy_ev = BOLTZMANN_EV_K * TEMPERATURE_K
    component_sigma_ms = np.sqrt(
        BOLTZMANN_J_K * TEMPERATURE_K / O2_ION_MASS_KG
    )
    rng = np.random.default_rng(RANDOM_SEED)
    sampled_velocity_ms = rng.normal(
        loc=0.0,
        scale=component_sigma_ms,
        size=(SAMPLE_COUNT, 3),
    )
    sampled_speed_ms = np.linalg.norm(sampled_velocity_ms, axis=1)
    sampled_energy_ev = (
        0.5
        * O2_ION_MASS_KG
        * sampled_speed_ms**2
        / 1.602176634e-19
    )
    sampled_direction_cosine = (
        sampled_velocity_ms[:, 2] / sampled_speed_ms
    )

    velocity_component_ms = np.linspace(
        -4.5 * component_sigma_ms,
        4.5 * component_sigma_ms,
        900,
    )
    component_pdf = np.exp(
        -0.5 * (velocity_component_ms / component_sigma_ms) ** 2
    ) / (np.sqrt(2.0 * np.pi) * component_sigma_ms)
    energy_ev = np.linspace(0.0, 8.0 * thermal_energy_ev, 900)
    energy_pdf = (
        2.0
        / np.sqrt(np.pi)
        * np.sqrt(energy_ev)
        / thermal_energy_ev**1.5
        * np.exp(-energy_ev / thermal_energy_ev)
    )

    fig, axes = plt.subplots(
        1, 3, figsize=(9.2, 3.15), constrained_layout=True,
    )
    axes[0].hist(
        sampled_velocity_ms[:, 0],
        bins=120,
        density=True,
        color="#9ecae1",
        edgecolor="none",
        alpha=0.65,
        label="Monte Carlo samples",
    )
    axes[0].plot(
        velocity_component_ms,
        component_pdf,
        color="#08519c",
        lw=1.6,
        label="Maxwellian theory",
    )
    axes[0].set(
        xlabel=r"$v_x$ (m s$^{-1}$)",
        ylabel=r"Probability density (s m$^{-1}$)",
        title=r"Velocity component, $T=300$ K",
        xlim=(velocity_component_ms[0], velocity_component_ms[-1]),
    )
    axes[0].legend()

    axes[1].hist(
        sampled_energy_ev,
        bins=120,
        density=True,
        color="#fdd0a2",
        edgecolor="none",
        alpha=0.65,
        label="Monte Carlo samples",
    )
    axes[1].plot(
        energy_ev,
        energy_pdf,
        color="#d94801",
        lw=1.6,
        label="Maxwellian theory",
    )
    axes[1].set(
        xlabel="Kinetic energy (eV)",
        ylabel=r"Probability density (eV$^{-1}$)",
        title="Total kinetic energy",
        xlim=(0.0, energy_ev[-1]),
    )
    axes[1].legend()

    axes[2].hist(
        sampled_direction_cosine,
        bins=100,
        density=True,
        color="#a1d99b",
        edgecolor="none",
        alpha=0.7,
        label="Monte Carlo samples",
    )
    axes[2].axhline(
        0.5, color="#238b45", lw=1.6, label="Isotropic theory",
    )
    axes[2].set(
        xlabel=r"Direction cosine, $\mu=v_z/|\mathbf{v}|$",
        ylabel="Probability density",
        title="Isotropic direction",
        xlim=(-1.0, 1.0),
        ylim=(0.0, 0.65),
    )
    axes[2].legend()

    for label, axis in zip(("a", "b", "c"), axes):
        axis.text(
            -0.12, 1.04, label, transform=axis.transAxes,
            fontsize=10, fontweight="bold", va="bottom",
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    component_integral = np.trapz(component_pdf, velocity_component_ms)
    energy_integral = np.trapz(energy_pdf, energy_ev)
    print(f"1D velocity PDF integral: {component_integral:.10f}")
    print(f"Energy PDF integral over plotted range: {energy_integral:.10f}")
    print(f"Mean sampled energy: {sampled_energy_ev.mean():.6f} eV")
    print(
        "Expected mean energy: "
        f"{1.5 * thermal_energy_ev:.6f} eV"
    )
    print(f"Wrote: {OUTPUT}")


if __name__ == "__main__":
    main()
