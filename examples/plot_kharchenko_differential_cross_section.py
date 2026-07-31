"""Explain the Kharchenko O + O differential scattering cross section."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ALPHA = 0.36e-16  # cm^2 sr^-1
BETA = -1.85


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def main() -> None:
    theta_deg = np.geomspace(0.01, 180.0, 3000)
    theta = np.deg2rad(theta_deg)

    dcs = ALPHA * np.sin(theta / 2.0) ** BETA
    # Axisymmetric solid-angle ring: dOmega = 2*pi*sin(theta)*dtheta.
    sigma_total = 8.0 * np.pi * ALPHA / (BETA + 2.0)
    pdf_per_rad = 2.0 * np.pi * np.sin(theta) * dcs / sigma_total
    pdf_per_deg = pdf_per_rad * np.pi / 180.0
    cdf = np.sin(theta / 2.0) ** (BETA + 2.0)

    fig = plt.figure(figsize=(10.5, 3.45), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.9])

    ax1 = fig.add_subplot(gs[0, 0])
    ax1.loglog(theta_deg, dcs, color="#176B87", lw=2.2)
    ax1.set_xlabel(r"COM scattering angle, $\theta$ (deg)")
    ax1.set_ylabel(r"$d\sigma/d\Omega$ (cm$^2$ sr$^{-1}$)")
    ax1.set_title("Differential cross section")
    ax1.grid(which="both", color="0.9", lw=0.6)
    ax1.text(
        0.05,
        0.06,
        r"$\frac{d\sigma}{d\Omega}=\alpha\sin^{\beta}(\theta/2)$"
        "\n"
        r"$\alpha=0.36\times10^{-16}$, $\beta=-1.85$",
        transform=ax1.transAxes,
        va="bottom",
    )

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.loglog(theta_deg, pdf_per_deg, color="#D55E00", lw=2.2)
    ax2.set_xlabel(r"COM scattering angle, $\theta$ (deg)")
    ax2.set_ylabel(r"Probability density (deg$^{-1}$)", color="#D55E00")
    ax2.tick_params(axis="y", colors="#D55E00")
    ax2.grid(which="both", color="0.9", lw=0.6)
    ax2b = ax2.twinx()
    ax2b.semilogx(theta_deg, cdf, color="#4C956C", lw=2.0)
    ax2b.set_ylim(0.0, 1.02)
    ax2b.set_ylabel(r"Cumulative probability, $P(\Theta\leq\theta)$", color="#4C956C")
    ax2b.tick_params(axis="y", colors="#4C956C")
    ax2b.spines["top"].set_visible(False)
    ax2.set_title("Angle probability includes $d\Omega$")

    ax3 = fig.add_subplot(gs[0, 2], projection="polar")
    ax3.set_theta_zero_location("N")
    ax3.set_theta_direction(-1)
    ax3.set_thetamin(-90)
    ax3.set_thetamax(90)
    ax3.set_rlim(0, 1.12)
    ax3.set_rticks([])
    ax3.set_xticks(np.deg2rad([-60, -30, 0, 30, 60, 90]))
    ax3.set_xticklabels([r"$60^\circ$", r"$30^\circ$", r"$0^\circ$", r"$30^\circ$", r"$60^\circ$", r"$90^\circ$"])
    th0, dth = np.deg2rad(42), np.deg2rad(9)
    ax3.bar(th0, 1.0, width=dth, bottom=0, color="#7DB7D5", alpha=0.65, edgecolor="#176B87")
    ax3.plot([0, 0], [0, 1.07], color="0.2", lw=1.5)
    ax3.plot([th0, th0], [0, 1.03], color="#176B87", lw=1.5)
    ax3.annotate(r"$\theta$", xy=(th0 / 2, 0.48), ha="center", color="#176B87", fontsize=12)
    ax3.text(th0, 0.78, r"ring $d\Omega$", ha="center", va="center", rotation=-42)
    ax3.set_title(r"$d\Omega$: a tiny direction patch", pad=15)

    fig.suptitle("Kharchenko O + O angular scattering model", fontsize=12, fontweight="bold")

    output = Path(__file__).resolve().parent / "figures"
    output.mkdir(parents=True, exist_ok=True)
    stem = output / "kharchenko_differential_cross_section_explained"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"sigma_total = {sigma_total:.6e} cm^2")
    print(f"CDF(10 deg) = {np.sin(np.deg2rad(10) / 2) ** (BETA + 2):.6f}")
    print(stem.with_suffix('.png'))


if __name__ == "__main__":
    main()
