"""Calculate directional hot O flux from particle crossing-event files."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIRECTORY = (
    ROOT / "examples" / "output" / "run_1p51m_crossings"
)
MARS_RADIUS_M = 3389.5e3
EV_J = 1.602176634e-19
AMU_KG = 1.66053906892e-27
O_MASS_KG = 15.999 * AMU_KG
HEADER = struct.Struct("<8sIIQQQQ16x")
EVENT_CROSSING = 2
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


def crossing_flux_for_file(
    path: Path,
    altitude_km: np.ndarray,
    energy_edges_eV: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    header = read_header(path)
    records = np.memmap(
        path,
        mode="r",
        dtype=EVENT_DTYPE,
        offset=HEADER.size,
        shape=(header["events"],),
    )
    crossing_indices = np.flatnonzero(
        records["event_code"] == EVENT_CROSSING
    )
    upward_rate_s1 = np.zeros(
        (altitude_km.size, energy_edges_eV.size - 1)
    )
    downward_rate_s1 = np.zeros_like(upward_rate_s1)

    chunk_size = 2_000_000
    for start in range(0, crossing_indices.size, chunk_size):
        indices = crossing_indices[start : start + chunk_size]
        event = records[indices]
        surface = event["surface_index"].astype(np.int64) - 1
        speed_squared = (
            event["velocity_x_m_s"] ** 2
            + event["velocity_y_m_s"] ** 2
            + event["velocity_z_m_s"] ** 2
        )
        energy_eV = 0.5 * O_MASS_KG * speed_squared / EV_J
        energy_bin = np.searchsorted(
            energy_edges_eV, energy_eV, side="right"
        ) - 1
        valid = (
            (surface >= 0)
            & (surface < altitude_km.size)
            & (energy_bin >= 0)
            & (energy_bin < energy_edges_eV.size - 1)
        )
        surface = surface[valid]
        energy_bin = energy_bin[valid]
        weight = event["weight_s1"][valid]
        direction = event["direction"][valid]
        upward = direction > 0
        np.add.at(
            upward_rate_s1,
            (surface[upward], energy_bin[upward]),
            weight[upward],
        )
        np.add.at(
            downward_rate_s1,
            (surface[~upward], energy_bin[~upward]),
            weight[~upward],
        )

    area_m2 = 4 * np.pi * (
        MARS_RADIUS_M + 1000 * altitude_km
    ) ** 2
    upward_flux_cm2_s1 = upward_rate_s1 / area_m2[:, None] / 1.0e4
    downward_flux_cm2_s1 = downward_rate_s1 / area_m2[:, None] / 1.0e4
    return upward_flux_cm2_s1, downward_flux_cm2_s1, header


def draw_map(
    flux: np.ndarray,
    altitude_km: np.ndarray,
    energy_eV: np.ndarray,
    title: str,
    output: Path,
    vmin: float,
    vmax: float,
) -> None:
    positive = flux[flux > 0]
    log_flux = np.log10(np.maximum(flux, positive.min()))
    figure, axis = plt.subplots(figsize=(5.2, 4.0), constrained_layout=True)
    image = axis.imshow(
        log_flux,
        origin="lower",
        aspect="auto",
        extent=(
            energy_eV[0],
            energy_eV[-1],
            altitude_km[0],
            altitude_km[-1],
        ),
        interpolation="bilinear",
        cmap="turbo",
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )
    axis.set(
        xlabel="Hot O energy (eV)",
        ylabel="Altitude (km)",
        title=title,
    )
    colorbar = figure.colorbar(image, ax=axis, pad=0.03)
    colorbar.set_label(
        r"$\log_{10}[\Phi_k\;(\mathrm{cm^{-2}\,s^{-1}\ per\ bin})]$"
    )
    figure.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(figure)


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
    paths = sorted(run_directory.glob("batch_*.bin"))
    if not paths:
        raise RuntimeError(f"No crossing-event files in {run_directory}")

    altitude_km = np.arange(100.0, 2000.0 + 10.0, 10.0)
    energy_edges_eV = np.linspace(0.01, 7.0, 141)
    energy_eV = (energy_edges_eV[:-1] + energy_edges_eV[1:]) / 2
    upward_batches = []
    downward_batches = []
    headers = []
    for path in paths:
        upward, downward, header = crossing_flux_for_file(
            path, altitude_km, energy_edges_eV,
        )
        upward_batches.append(upward)
        downward_batches.append(downward)
        headers.append(header)

    upward_stack = np.stack(upward_batches)
    downward_stack = np.stack(downward_batches)
    upward_flux = np.mean(upward_stack, axis=0)
    downward_flux = np.mean(downward_stack, axis=0)
    upward_standard_error = np.std(
        upward_stack, axis=0, ddof=1,
    ) / np.sqrt(len(paths))
    downward_standard_error = np.std(
        downward_stack, axis=0, ddof=1,
    ) / np.sqrt(len(paths))

    np.savez_compressed(
        run_directory / "hot_o_directional_crossing_flux_1p51m.npz",
        altitude_km=altitude_km,
        energy_eV=energy_eV,
        upward_flux_cm2_s1_per_bin=upward_flux,
        downward_flux_cm2_s1_per_bin=downward_flux,
        net_upward_flux_cm2_s1_per_bin=upward_flux - downward_flux,
        upward_standard_error_cm2_s1_per_bin=upward_standard_error,
        downward_standard_error_cm2_s1_per_bin=downward_standard_error,
    )
    output_table = np.column_stack(
        [
            np.repeat(altitude_km, energy_eV.size),
            np.tile(energy_eV, altitude_km.size),
            upward_flux.ravel(),
            downward_flux.ravel(),
            (upward_flux - downward_flux).ravel(),
            upward_standard_error.ravel(),
            downward_standard_error.ravel(),
        ]
    )
    np.savetxt(
        run_directory / "hot_o_directional_crossing_flux_1p51m.dat",
        output_table,
        header=(
            "altitude_km energy_eV upward_flux_cm-2_s-1_per_bin "
            "downward_flux_cm-2_s-1_per_bin "
            "net_upward_flux_cm-2_s-1_per_bin "
            "upward_standard_error_cm-2_s-1_per_bin "
            "downward_standard_error_cm-2_s-1_per_bin"
        ),
    )

    configure_matplotlib()
    positive = np.concatenate(
        [upward_flux[upward_flux > 0], downward_flux[downward_flux > 0]]
    )
    vmin = float(np.percentile(np.log10(positive), 2))
    vmax = float(np.percentile(np.log10(positive), 99))
    draw_map(
        upward_flux,
        altitude_km,
        energy_eV,
        "Upward hot O radial flux",
        run_directory / "upward_hot_o_crossing_flux_1p51m.png",
        vmin,
        vmax,
    )
    draw_map(
        downward_flux,
        altitude_km,
        energy_eV,
        "Downward hot O radial flux",
        run_directory / "downward_hot_o_crossing_flux_1p51m.png",
        vmin,
        vmax,
    )

    figure, axes = plt.subplots(
        1, 2, figsize=(8.2, 4.0), constrained_layout=True, sharey=True,
    )
    image = None
    for axis, flux, title, label in zip(
        axes,
        (upward_flux, downward_flux),
        ("Upward, $v_r>0$", "Downward, $v_r<0$"),
        ("a", "b"),
    ):
        positive_flux = flux[flux > 0]
        log_flux = np.log10(np.maximum(flux, positive_flux.min()))
        image = axis.imshow(
            log_flux,
            origin="lower",
            aspect="auto",
            extent=(
                energy_eV[0],
                energy_eV[-1],
                altitude_km[0],
                altitude_km[-1],
            ),
            interpolation="bilinear",
            cmap="turbo",
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        axis.set(xlabel="Hot O energy (eV)", title=title)
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
        )
    axes[0].set_ylabel("Altitude (km)")
    colorbar = figure.colorbar(image, ax=axes, pad=0.02)
    colorbar.set_label(
        r"$\log_{10}[\Phi_k\;(\mathrm{cm^{-2}\,s^{-1}\ per\ bin})]$"
    )
    figure.suptitle(
        "Directional hot O flux from altitude-surface crossings",
        fontsize=10,
    )
    figure.savefig(
        run_directory / "upward_downward_hot_o_crossing_flux_1p51m.png",
        dpi=400,
        bbox_inches="tight",
    )
    plt.close(figure)

    metadata = {
        "batch_files": len(paths),
        "primary_particles": sum(item["primary"] for item in headers),
        "secondary_particles": sum(item["secondary"] for item in headers),
        "tracked_particles": sum(item["tracked"] for item in headers),
        "event_records": sum(item["events"] for item in headers),
        "domain_km": [100.0, 2000.0],
        "crossing_spacing_km": 10.0,
        "energy_bins": 140,
        "flux_definition": "sum(weight_s1) / (4 pi r^2)",
        "flux_unit": "cm^-2 s^-1 per energy bin",
    }
    (run_directory / "directional_flux_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
