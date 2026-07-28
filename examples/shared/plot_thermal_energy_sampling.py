"""Plot the zero-mode half-normal energy model and inverse-CDF mapping."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import halfnorm


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "examples" / "figures" / "thermal_energy_sampling_300K.png"
TEMPERATURE_K = 300.0
BOLTZMANN_EV_K = 8.617333262145e-5
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
    thermal_scale_ev = BOLTZMANN_EV_K * TEMPERATURE_K
    distribution = halfnorm(loc=0.0, scale=thermal_scale_ev)
    rng = np.random.default_rng(RANDOM_SEED)
    sampled_energy_ev = np.abs(
        rng.normal(
            loc=0.0,
            scale=thermal_scale_ev,
            size=SAMPLE_COUNT,
        )
    )

    energy_ev = np.linspace(0.0, distribution.ppf(0.9999), 700)
    quantile = np.linspace(1.0e-6, 0.9999, 700)
    quantile_energy_ev = distribution.ppf(quantile)

    fig, axes = plt.subplots(
        1, 2, figsize=(7.2, 3.15), constrained_layout=True,
    )
    axes[0].hist(
        sampled_energy_ev,
        bins=120,
        density=True,
        color="#9ecae1",
        edgecolor="none",
        alpha=0.65,
        label="Monte Carlo samples",
    )
    axes[0].plot(
        energy_ev,
        distribution.pdf(energy_ev),
        color="#08519c",
        lw=1.6,
        label="Half-normal PDF",
    )
    axes[0].set(
        xlabel="Kinetic energy (eV)",
        ylabel=r"Probability density (eV$^{-1}$)",
        title=r"Zero-mode energy distribution, $T=300$ K",
        xlim=(0.0, energy_ev[-1]),
    )
    axes[0].legend()

    axes[1].plot(quantile, quantile_energy_ev, color="#d94801", lw=1.7)
    axes[1].axvline(0.5, color="0.45", lw=0.8, ls="--")
    axes[1].scatter(
        [0.5], [distribution.ppf(0.5)],
        color="#d94801", edgecolor="white", linewidth=0.5, s=24, zorder=3,
    )
    axes[1].annotate(
        f"median = {distribution.ppf(0.5):.3f} eV",
        xy=(0.5, distribution.ppf(0.5)),
        xytext=(0.57, 0.38),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "-", "color": "0.35", "lw": 0.7},
    )
    axes[1].set(
        xlabel=r"Uniform random quantile $u$",
        ylabel="Sampled kinetic energy (eV)",
        title=r"Inverse CDF, $E=F_E^{-1}(u)$",
        xlim=(0.0, 1.0),
        ylim=(0.0, quantile_energy_ev[-1]),
    )

    for label, axis in zip(("a", "b"), axes):
        axis.text(
            -0.12, 1.04, label, transform=axis.transAxes,
            fontsize=10, fontweight="bold", va="bottom",
        )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Mean sampled energy: {sampled_energy_ev.mean():.6f} eV")
    print(
        "Expected mean energy: "
        f"{np.sqrt(2.0 / np.pi) * thermal_scale_ev:.6f} eV"
    )
    print(f"Wrote: {OUTPUT}")


if __name__ == "__main__":
    main()
