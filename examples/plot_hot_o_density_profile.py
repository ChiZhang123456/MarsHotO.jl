"""Plot the energy-integrated hot O number density versus altitude."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIRECTORY = (
    ROOT / "examples" / "output" / "run_1p51m_current"
)


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
            "legend.fontsize": 8,
        }
    )


def load_density_grid(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = np.loadtxt(path)
    altitude_km = np.unique(table[:, 0])
    energy_eV = np.unique(table[:, 1])
    expected_rows = altitude_km.size * energy_eV.size
    if table.shape[0] != expected_rows:
        raise RuntimeError(
            f"Unexpected grid size in {path}: "
            f"{table.shape[0]} != {expected_rows}"
        )
    density_m3_per_bin = table[:, 2].reshape(
        altitude_km.size, energy_eV.size
    )
    if np.any(density_m3_per_bin < 0):
        raise RuntimeError(f"Negative density in {path}")
    return altitude_km, energy_eV, density_m3_per_bin


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "run_directory",
        nargs="?",
        type=Path,
        default=DEFAULT_RUN_DIRECTORY,
    )
    arguments = parser.parse_args()
    run_directory = arguments.run_directory.resolve()
    batch_paths = sorted(run_directory.glob("batch_??.dat"))
    if not batch_paths:
        raise RuntimeError(f"No batch density files in {run_directory}")

    altitude_reference = None
    energy_reference = None
    total_density_batches_cm3 = []
    for path in batch_paths:
        altitude_km, energy_eV, density_m3_per_bin = (
            load_density_grid(path)
        )
        if altitude_reference is None:
            altitude_reference = altitude_km
            energy_reference = energy_eV
        else:
            np.testing.assert_allclose(
                altitude_km, altitude_reference
            )
            np.testing.assert_allclose(energy_eV, energy_reference)
        total_density_batches_cm3.append(
            np.sum(density_m3_per_bin, axis=1) / 1.0e6
        )

    density_stack_cm3 = np.stack(total_density_batches_cm3)
    density_mean_cm3 = np.mean(density_stack_cm3, axis=0)
    density_standard_error_cm3 = np.std(
        density_stack_cm3, axis=0, ddof=1
    ) / np.sqrt(density_stack_cm3.shape[0])
    if np.any(density_mean_cm3 <= 0):
        raise RuntimeError("The total density profile must be positive")

    np.savez_compressed(
        run_directory / "hot_o_density_profile_1p51m.npz",
        altitude_km=altitude_reference,
        density_cm3=density_mean_cm3,
        standard_error_cm3=density_standard_error_cm3,
        batch_count=density_stack_cm3.shape[0],
    )

    configure_matplotlib()
    figure, axis = plt.subplots(
        figsize=(4.3, 4.8), constrained_layout=True
    )
    density_lower = np.maximum(
        density_mean_cm3 - density_standard_error_cm3,
        np.finfo(float).tiny,
    )
    density_upper = (
        density_mean_cm3 + density_standard_error_cm3
    )
    axis.fill_betweenx(
        altitude_reference,
        density_lower,
        density_upper,
        color="#2166ac",
        alpha=0.25,
        linewidth=0,
        label="Standard error",
    )
    axis.plot(
        density_mean_cm3,
        altitude_reference,
        color="#2166ac",
        linewidth=1.6,
        label="Batch mean",
    )
    axis.set(
        xscale="log",
        xlabel=r"Hot O number density (cm$^{-3}$)",
        ylabel="Altitude (km)",
        title="Hot O number density profile",
        ylim=(
            float(altitude_reference[0]),
            float(altitude_reference[-1]),
        ),
    )
    axis.grid(True, color="0.88", linewidth=0.6)
    axis.legend(frameon=False, loc="upper right")
    axis.text(
        0.97,
        0.83,
        (
            f"{density_stack_cm3.shape[0]} independent batches\n"
            "1,510,000 primary particles"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        color="0.25",
    )

    output = run_directory / "hot_o_density_profile_1p51m.png"
    figure.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(figure)
    print(f"output={output}")
    print(
        "peak_altitude_km="
        f"{altitude_reference[np.argmax(density_mean_cm3)]:.1f}"
    )
    print(f"peak_density_cm3={np.max(density_mean_cm3):.8e}")


if __name__ == "__main__":
    main()
