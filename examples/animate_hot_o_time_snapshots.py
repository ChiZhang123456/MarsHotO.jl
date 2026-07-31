"""Animate the altitude and energy evolution of hot O snapshot flux."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT
    / "examples"
    / "output"
    / "hot_o_time_animation_snapshots.dat"
)
DEFAULT_OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "hot_o_energy_altitude_time_evolution_legacy_flux.gif"
)
MARS_RADIUS_M = 3389.5e3


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fps", type=float, default=10.0)
    parser.add_argument("--minimum-particles", type=int, default=20)
    arguments = parser.parse_args()

    table = np.loadtxt(arguments.input)
    time_s = np.unique(table[:, 0])
    altitude_km = np.unique(table[:, 1])
    energy_eV = np.unique(table[:, 2])
    shape = (time_s.size, altitude_km.size, energy_eV.size)
    particle_count = table[:, 4].reshape(shape)
    weighted_rate_s1 = table[:, 5].reshape(shape)
    altitude_count = particle_count.sum(axis=2)

    area_m2 = 4.0 * np.pi * (
        MARS_RADIUS_M + 1000.0 * altitude_km
    ) ** 2
    flux_cm2_s1 = weighted_rate_s1 / area_m2[None, :, None] / 1.0e4
    resolved = (
        altitude_count[:, :, None] >= arguments.minimum_particles
    ) & (flux_cm2_s1 > 0)
    positive_flux = flux_cm2_s1[resolved]
    if positive_flux.size == 0:
        raise RuntimeError("No resolved positive snapshot flux")

    log_flux = np.full(flux_cm2_s1.shape, np.nan)
    log_flux[resolved] = np.log10(flux_cm2_s1[resolved])
    frame_color_min = []
    frame_color_max = []
    for frame_index in range(time_s.size):
        frame_values = log_flux[frame_index][
            np.isfinite(log_flux[frame_index])
        ]
        if frame_values.size > 0:
            frame_color_min.append(np.percentile(frame_values, 1.0))
            frame_color_max.append(np.percentile(frame_values, 99.5))
    color_min = float(np.min(frame_color_min))
    color_max = float(np.max(frame_color_max))

    energy_step = float(np.median(np.diff(energy_eV)))
    altitude_step = float(np.median(np.diff(altitude_km)))
    energy_edges = np.r_[
        energy_eV - energy_step / 2,
        energy_eV[-1] + energy_step / 2,
    ]
    altitude_edges = np.r_[
        altitude_km - altitude_step / 2,
        altitude_km[-1] + altitude_step / 2,
    ]

    configure_matplotlib()
    colormap = mpl.colormaps["turbo"].copy()
    colormap.set_bad(colormap(0.0))
    figure, axis = plt.subplots(
        figsize=(5.8, 4.5), constrained_layout=True
    )
    image = axis.pcolormesh(
        energy_edges,
        altitude_edges,
        log_flux[0],
        cmap=colormap,
        vmin=color_min,
        vmax=color_max,
        shading="auto",
        rasterized=True,
    )
    axis.set(
        xlim=(0.0, 7.0),
        ylim=(100.0, 1000.0),
        xlabel="Hot O energy (eV)",
        ylabel="Altitude (km)",
    )
    title = axis.set_title(rf"$t={time_s[0]:.0f}$ s")
    colorbar = figure.colorbar(image, ax=axis, pad=0.025)
    colorbar.set_label(
        r"$\log_{10}[\Phi_k\;(\mathrm{cm^{-2}\,s^{-1}\ per\ bin})]$"
    )
    figure.suptitle(
        "Area normalized hot O snapshot flux",
        fontsize=11,
    )

    def update(frame_index: int):
        image.set_array(log_flux[frame_index].ravel())
        title.set_text(rf"$t={time_s[frame_index]:.0f}$ s")
        return image, title

    animation = FuncAnimation(
        figure,
        update,
        frames=time_s.size,
        interval=1000.0 / arguments.fps,
        blit=False,
        repeat=True,
    )
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    animation.save(
        arguments.output,
        writer=PillowWriter(fps=arguments.fps),
        dpi=120,
    )
    plt.close(figure)

    print(f"output={arguments.output.resolve()}")
    print(f"frames={time_s.size}")
    print(f"time_range_s={time_s[0]:.6g},{time_s[-1]:.6g}")
    print(f"log10_flux_color_range={color_min:.6g},{color_max:.6g}")


if __name__ == "__main__":
    main()
