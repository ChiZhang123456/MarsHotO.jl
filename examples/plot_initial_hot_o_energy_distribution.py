"""Plot the sampled hot O energy distribution at Monte Carlo initialization."""

from __future__ import annotations

import argparse
import struct
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIRECTORY = (
    ROOT / "examples" / "output" / "run_paired_dr_crossings"
)
DEFAULT_OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "hot_o_initial_energy_altitude_distribution.png"
)

EV_J = 1.602176634e-19
AMU_KG = 1.66053906892e-27
O_MASS_KG = 15.999 * AMU_KG
HEADER = struct.Struct("<8sIIQQQQ16x")
EVENT_BIRTH = 1
EVENT_DTYPE = np.dtype(
    {
        "names": [
            "particle_id",
            "parent_id",
            "weight_s1",
            "time_s",
            "altitude_km",
            "velocity_x_m_s",
            "velocity_y_m_s",
            "velocity_z_m_s",
            "radial_velocity_m_s",
            "event_index",
            "collisions",
            "surface_index",
            "event_code",
            "direction",
        ],
        "formats": [
            "<i8",
            "<i8",
            "<f8",
            "<f8",
            "<f8",
            "<f8",
            "<f8",
            "<f8",
            "<f8",
            "<i4",
            "<i4",
            "<i2",
            "i1",
            "i1",
        ],
        "offsets": [
            0,
            8,
            16,
            24,
            32,
            40,
            48,
            56,
            64,
            72,
            76,
            80,
            82,
            83,
        ],
        "itemsize": 88,
    }
)


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


def read_header(path: Path) -> dict[str, int]:
    with path.open("rb") as stream:
        values = HEADER.unpack(stream.read(HEADER.size))
    magic, version, record_size, events, primary, secondary, tracked = values
    if magic != b"MHOTE001":
        raise RuntimeError(f"Unexpected event-file magic in {path}")
    if version != 1 or record_size != EVENT_DTYPE.itemsize:
        raise RuntimeError(
            f"Unsupported version or record size in {path}: "
            f"{version}, {record_size}"
        )
    expected_size = HEADER.size + events * record_size
    if path.stat().st_size != expected_size:
        raise RuntimeError(
            f"Incomplete event file {path}: "
            f"{path.stat().st_size} != {expected_size}"
        )
    return {
        "events": events,
        "primary": primary,
        "secondary": secondary,
        "tracked": tracked,
    }


def accumulate_primary_births(
    path: Path,
    altitude_edges_km: np.ndarray,
    energy_edges_eV: np.ndarray,
    counts: np.ndarray,
) -> int:
    header = read_header(path)
    records = np.memmap(
        path,
        mode="r",
        dtype=EVENT_DTYPE,
        offset=HEADER.size,
        shape=(header["events"],),
    )
    selected = (
        (records["event_code"] == EVENT_BIRTH)
        & (records["parent_id"] == 0)
        & (records["time_s"] == 0.0)
    )
    birth = records[selected]
    speed_squared = (
        birth["velocity_x_m_s"] ** 2
        + birth["velocity_y_m_s"] ** 2
        + birth["velocity_z_m_s"] ** 2
    )
    energy_eV = 0.5 * O_MASS_KG * speed_squared / EV_J
    altitude_bin = (
        np.searchsorted(
            altitude_edges_km, birth["altitude_km"], side="right"
        )
        - 1
    )
    energy_bin = (
        np.searchsorted(energy_edges_eV, energy_eV, side="right") - 1
    )
    valid = (
        (altitude_bin >= 0)
        & (altitude_bin < counts.shape[0])
        & (energy_bin >= 0)
        & (energy_bin < counts.shape[1])
    )
    np.add.at(counts, (altitude_bin[valid], energy_bin[valid]), 1)
    return int(np.count_nonzero(valid))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "run_directory",
        nargs="?",
        type=Path,
        default=DEFAULT_RUN_DIRECTORY,
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    run_directory = arguments.run_directory.resolve()
    paths = sorted(run_directory.glob("batch_*.bin"))
    if not paths:
        raise RuntimeError(f"No crossing-event files in {run_directory}")

    altitude_centers_km = np.arange(100.0, 251.0, 1.0)
    altitude_edges_km = np.arange(99.5, 250.5 + 1.0e-9, 1.0)
    energy_edges_eV = np.arange(0.0, 7.0 + 0.05, 0.05)
    counts = np.zeros(
        (altitude_centers_km.size, energy_edges_eV.size - 1),
        dtype=np.int64,
    )

    total_primary_births = 0
    for path in paths:
        total_primary_births += accumulate_primary_births(
            path, altitude_edges_km, energy_edges_eV, counts
        )

    particles_per_altitude = counts.sum(axis=1)
    if np.any(particles_per_altitude == 0):
        missing = altitude_centers_km[particles_per_altitude == 0]
        raise RuntimeError(f"No primary births at altitudes {missing}")
    probability = counts / particles_per_altitude[:, None]
    row_sums = probability.sum(axis=1)
    if not np.allclose(row_sums, 1.0, rtol=0.0, atol=1.0e-12):
        raise RuntimeError("Initial energy probabilities do not sum to one")

    configure_matplotlib()
    figure, axis = plt.subplots(
        figsize=(5.8, 4.3), constrained_layout=True
    )
    image = axis.pcolormesh(
        energy_edges_eV,
        altitude_edges_km,
        probability,
        cmap="turbo",
        vmin=0.0,
        vmax=0.06,
        shading="auto",
        rasterized=True,
    )
    axis.set(
        xlim=(0.0, 7.0),
        ylim=(100.0, 251.0),
        xlabel="Initial hot O energy (eV)",
        ylabel="Source altitude (km)",
        title=r"Monte Carlo initial hot O distribution, $t=0$ s",
    )
    colorbar = figure.colorbar(image, ax=axis, pad=0.025)
    colorbar.set_label("Probability per 0.05 eV bin")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.output, dpi=400, bbox_inches="tight")
    plt.close(figure)

    print(f"output={arguments.output.resolve()}")
    print(f"batch_files={len(paths)}")
    print(f"primary_births={total_primary_births}")
    print(
        "particles_per_altitude_min_max="
        f"{particles_per_altitude.min()},{particles_per_altitude.max()}"
    )
    print(f"probability_row_sum_min_max={row_sums.min():.16g},"
          f"{row_sums.max():.16g}")
    print(f"maximum_probability_per_bin={probability.max():.6f}")
    peak_indices = np.argmax(probability, axis=1)
    energy_centers_eV = (
        energy_edges_eV[:-1] + energy_edges_eV[1:]
    ) / 2
    print(
        "peak_energy_eV_at_100_180_250km="
        + ",".join(
            f"{energy_centers_eV[peak_indices[index]]:.3f}"
            for index in (0, 80, 150)
        )
    )


if __name__ == "__main__":
    main()
