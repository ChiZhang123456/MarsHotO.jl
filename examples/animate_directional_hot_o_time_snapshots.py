"""Animate upward and downward hot O energy-altitude distributions."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples/output/directional_hot_o_time_snapshots.dat"
OUTPUT = ROOT / "examples/figures/hot_o_energy_altitude_time_evolution.gif"
PREVIEW = ROOT / "examples/figures/directional_hot_o_time_snapshots_50s.png"

mpl.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "dejavusans",
    "font.size": 10,
})

table = np.loadtxt(INPUT)
times = np.unique(table[:, 0])
altitudes = np.unique(table[:, 1])
energies = np.unique(table[:, 2])
shape = (times.size, altitudes.size, energies.size)
upward = table[:, 3].reshape(shape)
downward = table[:, 4].reshape(shape)
positive = np.concatenate((upward[upward > 0], downward[downward > 0]))
vmin, vmax = np.percentile(np.log10(positive), [1, 99.8])
energy_edges = np.r_[energies - 0.05, energies[-1] + 0.05]
altitude_edges = np.r_[altitudes - 2.5, altitudes[-1] + 2.5]

figure, axes = plt.subplots(1, 2, figsize=(10.5, 5.0), constrained_layout=True)
images = []
for axis, values, title in zip(axes, (upward, downward), ("Upward", "Downward")):
    image = axis.pcolormesh(
        energy_edges, altitude_edges,
        np.log10(np.maximum(values[0], 10**vmin)),
        shading="flat", cmap="turbo", vmin=vmin, vmax=vmax,
    )
    axis.set(xlabel="Hot O energy (eV)", ylabel="Altitude (km)", title=title)
    axis.set_xlim(0, 6)
    axis.set_ylim(100, 400)
    images.append(image)
colorbar = figure.colorbar(images[0], ax=axes, pad=0.02)
colorbar.set_label(r"$\log_{10}[n(E,z)]$ (cm$^{-3}$ eV$^{-1}$)")
suptitle = figure.suptitle("")

def update(frame: int):
    for image, values in zip(images, (upward, downward)):
        image.set_array(np.log10(np.maximum(values[frame], 10**vmin)).ravel())
    suptitle.set_text(
        f"Directional hot O distribution, t = {times[frame]:g} s\n"
        "10,000 DR events per 1 km source altitude"
    )
    return (*images, suptitle)

animation = FuncAnimation(figure, update, frames=times.size, interval=350, blit=False)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
animation.save(OUTPUT, writer=PillowWriter(fps=3), dpi=150)
update(int(np.argmin(np.abs(times - 50))))
figure.savefig(PREVIEW, dpi=300)
plt.close(figure)
print(f"animation={OUTPUT.resolve()}")
print(f"preview={PREVIEW.resolve()}")
