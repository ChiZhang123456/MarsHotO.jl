"""Plot weighted hot O energy distributions at selected model times."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples" / "output" / "hot_o_time_snapshots.dat"
OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "hot_o_energy_altitude_time_snapshots.png"
)
TIMES_S = np.array([0.0, 10.0, 50.0, 100.0])
PANEL_LABELS = ("a", "b", "c", "d")


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        }
    )


def main() -> None:
    table = np.loadtxt(INPUT)
    time_s = np.unique(table[:, 0])
    altitude_km = np.unique(table[:, 1])
    energy_eV = np.unique(table[:, 2])
    np.testing.assert_allclose(time_s, TIMES_S)

    probability = table[:, 3].reshape(
        time_s.size, altitude_km.size, energy_eV.size
    )
    particle_count = table[:, 4].reshape(probability.shape)
    altitude_count = particle_count.sum(axis=2)

    energy_step = float(np.median(np.diff(energy_eV)))
    altitude_step = float(np.median(np.diff(altitude_km)))
    energy_edges = np.r_[
        energy_eV - energy_step / 2, energy_eV[-1] + energy_step / 2
    ]
    altitude_edges = np.r_[
        altitude_km - altitude_step / 2,
        altitude_km[-1] + altitude_step / 2,
    ]

    configure_matplotlib()
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(7.3, 6.2),
        constrained_layout=True,
        sharex=True,
        sharey=True,
    )
    colormap = mpl.colormaps["turbo"].copy()
    colormap.set_bad(colormap(0.0))
    image = None
    for index, axis in enumerate(axes.flat):
        display = probability[index].copy()
        display[altitude_count[index] < 20, :] = np.nan
        image = axis.pcolormesh(
            energy_edges,
            altitude_edges,
            display,
            cmap=colormap,
            vmin=0.0,
            vmax=0.06,
            shading="auto",
            rasterized=True,
        )
        axis.set(
            xlim=(0.0, 7.0),
            ylim=(100.0, 1000.0),
            title=rf"$t={TIMES_S[index]:.0f}$ s",
        )
        axis.text(
            0.025,
            0.965,
            PANEL_LABELS[index],
            transform=axis.transAxes,
            va="top",
            ha="left",
            fontsize=10,
            fontweight="bold",
            color="white",
        )

    for axis in axes[:, 0]:
        axis.set_ylabel("Altitude (km)")
    for axis in axes[-1, :]:
        axis.set_xlabel("Hot O energy (eV)")

    colorbar = figure.colorbar(
        image,
        ax=axes,
        location="right",
        pad=0.02,
        shrink=0.98,
    )
    colorbar.set_label("Weighted probability per 0.05 eV bin")
    figure.suptitle(
        "Evolution of an initially released hot O ensemble",
        fontsize=11,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)

    print(f"output={OUTPUT.resolve()}")
    for index, snapshot_time_s in enumerate(TIMES_S):
        valid_altitudes = altitude_count[index] >= 20
        row_sums = probability[index, valid_altitudes].sum(axis=1)
        print(
            f"time={snapshot_time_s:.0f}s, "
            f"particles={int(particle_count[index].sum())}, "
            f"resolved_altitude_bins={int(valid_altitudes.sum())}, "
            f"row_sum_min_max={row_sums.min():.12g},"
            f"{row_sums.max():.12g}"
        )


if __name__ == "__main__":
    main()
