"""Plot the MarsHotO density in each altitude and energy bin."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "examples" / "output" / "hot_o_altitude_energy_distribution.dat"
OUTPUT = ROOT / "examples" / "figures" / "hot_o_altitude_energy_distribution.png"


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
        }
    )


def main() -> None:
    table = np.loadtxt(INPUT)
    altitude = np.unique(table[:, 0])
    energy = np.unique(table[:, 1])
    density = table[:, 2].reshape(altitude.size, energy.size)
    positive = density[density > 0]
    if positive.size == 0:
        raise RuntimeError("The Monte Carlo density grid contains no positive values")

    log_density = np.log10(np.maximum(density, positive.min()))
    configure_matplotlib()
    figure, axis = plt.subplots(figsize=(5.2, 4.0), constrained_layout=True)
    image = axis.imshow(
        log_density,
        origin="lower",
        aspect="auto",
        extent=(energy[0], energy[-1], altitude[0], altitude[-1]),
        interpolation="bilinear",
        cmap="turbo",
        vmin=float(np.percentile(np.log10(positive), 2)),
        vmax=float(np.percentile(np.log10(positive), 99)),
        rasterized=True,
    )
    axis.set(
        xlabel="Hot O energy (eV)",
        ylabel="Altitude (km)",
        title="Monte Carlo hot O altitude-energy distribution\n"
        "Spherically extended subsolar MGITM profile",
    )
    colorbar = figure.colorbar(image, ax=axis, pad=0.03)
    colorbar.set_label(
        r"$\log_{10}\,[n_k(z)\;(\mathrm{m^{-3}})]$"
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
