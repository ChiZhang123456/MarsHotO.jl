"""Plot weighted hot O energy probability density near 200 km."""

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
    / "hot_o_energy_probability_density_200km.png"
)
SELECTED_TIMES_S = (0.0, 50.0, 100.0)
ALTITUDE_MIN_KM = 195.0
ALTITUDE_MAX_KM = 205.0


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "axes.linewidth": 0.8,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "legend.frameon": False,
        }
    )


def main() -> None:
    table = np.loadtxt(INPUT)
    time_s = np.unique(table[:, 0])
    altitude_km = np.unique(table[:, 1])
    energy_eV = np.unique(table[:, 2])
    shape = (time_s.size, altitude_km.size, energy_eV.size)
    particle_count = table[:, 4].reshape(shape)
    weighted_rate_s1 = table[:, 5].reshape(shape)

    energy_step_eV = float(np.median(np.diff(energy_eV)))
    altitude_indices = np.flatnonzero(
        (altitude_km > ALTITUDE_MIN_KM)
        & (altitude_km < ALTITUDE_MAX_KM)
    )
    if altitude_indices.size == 0:
        raise RuntimeError("No altitude bins found near 200 km")

    configure_matplotlib()
    figure, axis = plt.subplots(
        figsize=(5.6, 4.0),
        constrained_layout=True,
    )
    colors = ("#202020", "#2878B5", "#D95319")

    for selected_time_s, color in zip(SELECTED_TIMES_S, colors):
        time_matches = np.flatnonzero(
            np.isclose(time_s, selected_time_s)
        )
        if time_matches.size != 1:
            raise RuntimeError(
                f"Snapshot t={selected_time_s:g} s is unavailable"
            )
        time_index = int(time_matches[0])
        weighted_spectrum = weighted_rate_s1[
            time_index, altitude_indices, :
        ].sum(axis=0)
        total_weight_s1 = float(weighted_spectrum.sum())
        if total_weight_s1 <= 0:
            raise RuntimeError(
                f"No weighted particles at t={selected_time_s:g} s"
            )

        probability_density_eV1 = (
            weighted_spectrum / total_weight_s1 / energy_step_eV
        )
        display_density = np.where(
            probability_density_eV1 > 0,
            probability_density_eV1,
            np.nan,
        )
        axis.plot(
            energy_eV,
            display_density,
            color=color,
            linewidth=1.8,
            drawstyle="steps-mid",
            label=rf"$t={selected_time_s:.0f}$ s",
        )

        integral = float(
            np.nansum(probability_density_eV1) * energy_step_eV
        )
        particles = int(
            particle_count[time_index, altitude_indices, :].sum()
        )
        mean_energy_eV = float(
            np.sum(
                probability_density_eV1 * energy_eV
            )
            * energy_step_eV
        )
        print(
            f"time={selected_time_s:.0f}s, "
            f"particles={particles}, "
            f"integral={integral:.12g}, "
            f"mean_energy_eV={mean_energy_eV:.6g}"
        )

    axis.set(
        xlim=(0.0, 7.0),
        ylim=(5.0e-5, 30.0),
        yscale="log",
        xlabel="Hot O energy (eV)",
        ylabel=r"Probability density ($\mathrm{eV^{-1}}$)",
        title="Hot O energy distribution near 200 km",
    )
    axis.grid(True, which="major", color="0.88", linewidth=0.6)
    axis.legend(loc="upper right")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(f"altitude_bin_centers_km={altitude_km[altitude_indices]}")
    print(f"energy_bin_width_eV={energy_step_eV:.6g}")
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
