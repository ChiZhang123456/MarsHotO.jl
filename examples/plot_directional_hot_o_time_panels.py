"""Plot upward and downward hot O distributions at 0, 50, and 100 s."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples/output/directional_hot_o_time_snapshots.dat"
OUTPUT = ROOT / "examples/figures/directional_hot_o_time_panels.png"
SELECTED_TIMES_S = (0.0, 50.0, 100.0)

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "mathtext.fontset": "dejavusans",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    }
)

table = np.loadtxt(INPUT)
times = np.unique(table[:, 0])
altitudes = np.unique(table[:, 1])
energies = np.unique(table[:, 2])
shape = (times.size, altitudes.size, energies.size)
upward = table[:, 3].reshape(shape)
downward = table[:, 4].reshape(shape)
indices = [int(np.argmin(np.abs(times - time))) for time in SELECTED_TIMES_S]

selected_positive = []
for index in indices:
    selected_positive.extend(upward[index][upward[index] > 0])
    selected_positive.extend(downward[index][downward[index] > 0])
log_positive = np.log10(np.asarray(selected_positive))
vmin, vmax = np.percentile(log_positive, [1.0, 99.8])
floor = 10**vmin

energy_step = float(np.median(np.diff(energies)))
altitude_step = float(np.median(np.diff(altitudes)))
energy_edges = np.r_[energies - energy_step / 2, energies[-1] + energy_step / 2]
altitude_edges = np.r_[
    altitudes - altitude_step / 2,
    altitudes[-1] + altitude_step / 2,
]

figure, axes = plt.subplots(
    3, 2, figsize=(9.0, 10.0), sharex=True, sharey=True,
    constrained_layout=True,
)
panel_labels = iter("abcdef")
image = None
for row, (time, index) in enumerate(zip(SELECTED_TIMES_S, indices)):
    for column, (values, direction) in enumerate(
        ((upward[index], "Upward"), (downward[index], "Downward"))
    ):
        axis = axes[row, column]
        image = axis.pcolormesh(
            energy_edges,
            altitude_edges,
            np.log10(np.maximum(values, floor)),
            shading="flat",
            cmap="turbo",
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        axis.set_xlim(0, 6)
        axis.set_ylim(100, 400)
        axis.set_title(f"{direction},  t = {time:g} s")
        axis.text(
            0.02, 0.96, next(panel_labels), transform=axis.transAxes,
            ha="left", va="top", fontsize=11, fontweight="bold",
            color="white",
        )
        if column == 0:
            axis.set_ylabel("Altitude (km)")
        if row == 2:
            axis.set_xlabel("Hot O energy (eV)")

colorbar = figure.colorbar(image, ax=axes, pad=0.02, shrink=0.96)
colorbar.set_label(r"$\log_{10}[n(E,z)]$ (cm$^{-3}$ eV$^{-1}$)")
figure.suptitle(
    "Directional hot O time evolution\n"
    "10,000 DR events per 1 km source altitude",
    fontsize=13,
)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
plt.close(figure)
print(f"output={OUTPUT.resolve()}")
