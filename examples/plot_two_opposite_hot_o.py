"""Plot trajectories and collision histories for two opposite hot O atoms."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
TRAJECTORY_INPUT = (
    ROOT / "examples" / "output" / "two_opposite_hot_o_trajectories.dat"
)
COLLISION_INPUT = (
    ROOT / "examples" / "output" / "two_opposite_hot_o_collisions.dat"
)
OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "two_opposite_hot_o_collision_trajectories.png"
)

SPECIES_NAMES = {1: "O", 2: "CO", 3: r"N$_2$", 4: r"O$_2$", 5: r"CO$_2$"}
SPECIES_COLORS = {
    1: "#6A3D9A",
    2: "#1F78B4",
    3: "#33A02C",
    4: "#FF7F00",
    5: "#E31A1C",
}
PARTICLE_COLORS = {1: "#1F77B4", 2: "#D62728"}
PARTICLE_MARKERS = {1: "o", 2: "^"}


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
            "legend.fontsize": 7,
            "legend.frameon": False,
        }
    )


def load_table(path: Path, columns: int) -> np.ndarray:
    table = np.loadtxt(path)
    if table.size == 0:
        return np.empty((0, columns))
    return np.atleast_2d(table)


def main() -> None:
    trajectory = load_table(TRAJECTORY_INPUT, 9)
    collisions = load_table(COLLISION_INPUT, 12)
    configure_matplotlib()

    figure = plt.figure(figsize=(8.2, 6.6), constrained_layout=True)
    grid = figure.add_gridspec(2, 2)
    axis_trajectory = figure.add_subplot(grid[0, 0], projection="3d")
    axis_energy = figure.add_subplot(grid[0, 1])
    axis_angle = figure.add_subplot(grid[1, 0])
    axis_loss = figure.add_subplot(grid[1, 1])

    for particle_id in (1, 2):
        rows = trajectory[trajectory[:, 0] == particle_id]
        axis_trajectory.plot(
            rows[:, 4],
            rows[:, 5],
            rows[:, 6],
            color=PARTICLE_COLORS[particle_id],
            linewidth=1.4,
            label=f"O{particle_id}",
        )
        axis_trajectory.scatter(
            rows[0, 4],
            rows[0, 5],
            rows[0, 6],
            color=PARTICLE_COLORS[particle_id],
            marker="*",
            s=45,
        )
        axis_energy.plot(
            rows[:, 2],
            rows[:, 7],
            color=PARTICLE_COLORS[particle_id],
            linewidth=1.4,
            label=f"O{particle_id}",
        )

    for species_code, species_name in SPECIES_NAMES.items():
        rows = collisions[collisions[:, 7] == species_code]
        if rows.size == 0:
            continue
        axis_trajectory.scatter(
            rows[:, 4],
            rows[:, 5],
            rows[:, 6],
            color=SPECIES_COLORS[species_code],
            edgecolor="black",
            linewidth=0.3,
            s=18,
            label=f"collision with {species_name}",
        )
        axis_energy.scatter(
            rows[:, 2],
            rows[:, 10],
            color=SPECIES_COLORS[species_code],
            edgecolor="black",
            linewidth=0.3,
            s=18,
            zorder=3,
        )

    for particle_id in (1, 2):
        rows = collisions[collisions[:, 0] == particle_id]
        if rows.size == 0:
            continue
        for species_code in SPECIES_NAMES:
            selected = rows[rows[:, 7] == species_code]
            if selected.size == 0:
                continue
            axis_angle.scatter(
                selected[:, 1],
                selected[:, 8],
                color=SPECIES_COLORS[species_code],
                marker=PARTICLE_MARKERS[particle_id],
                edgecolor="black",
                linewidth=0.3,
                s=28,
            )
            axis_loss.scatter(
                selected[:, 8],
                100 * selected[:, 11],
                color=SPECIES_COLORS[species_code],
                marker=PARTICLE_MARKERS[particle_id],
                edgecolor="black",
                linewidth=0.3,
                s=28,
            )
        axis_angle.plot(
            rows[:, 1],
            rows[:, 8],
            color=PARTICLE_COLORS[particle_id],
            linewidth=0.7,
            alpha=0.6,
        )

    axis_trajectory.set(
        xlabel="East displacement (km)",
        ylabel="North displacement (km)",
        zlabel="Altitude (km)",
        title="Particle trajectories",
    )
    axis_trajectory.view_init(elev=24, azim=-55)
    axis_trajectory.legend(loc="upper left", bbox_to_anchor=(-0.05, 1.0))

    axis_energy.set(
        xlabel="Elapsed time (s)",
        ylabel="Kinetic energy (eV)",
        title="Energy history",
    )
    axis_energy.grid(True, color="0.88", linewidth=0.6)
    axis_energy.legend(loc="best")

    axis_angle.set(
        xlabel="Collision number",
        ylabel=r"COM scattering angle, $\theta$ (deg)",
        title="Sampled scattering angles",
    )
    axis_angle.grid(True, color="0.88", linewidth=0.6)

    axis_loss.set(
        xlabel=r"COM scattering angle, $\theta$ (deg)",
        ylabel="Projectile energy loss (%)",
        title="Energy loss at each collision",
    )
    axis_loss.grid(True, color="0.88", linewidth=0.6)

    for label, axis in zip(
        ("a", "b", "c", "d"),
        (axis_trajectory, axis_energy, axis_angle, axis_loss),
    ):
        text_method = axis.text2D if axis is axis_trajectory else axis.text
        text_method(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
        )

    particle_handles = [
        mpl.lines.Line2D(
            [], [], color="black", marker=PARTICLE_MARKERS[particle_id],
            linestyle="none", markersize=5, label=f"O{particle_id}",
        )
        for particle_id in (1, 2)
    ]
    species_handles = [
        mpl.lines.Line2D(
            [], [], color=SPECIES_COLORS[code], marker="o",
            linestyle="none", markersize=5, label=name,
        )
        for code, name in SPECIES_NAMES.items()
        if np.any(collisions[:, 7] == code)
    ]
    axis_loss.legend(
        handles=particle_handles + species_handles,
        loc="best",
        ncol=2,
    )

    figure.suptitle(
        "Two opposite hot O atoms: gravity and neutral collisions\n"
        r"$z_0=180$ km, $E_0=3.495$ eV per O",
        fontsize=10,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
