"""Calculate directional hot O flux from particle crossing-event files."""

from __future__ import annotations

import argparse
import json
import shutil
import struct
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIRECTORY = (
    ROOT / "examples" / "output" / "run_paired_dr_crossings"
)
MARS_RADIUS_M = 3389.5e3
MARS_MU_M3_S2 = 4.282837e13
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
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
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
    total_stack = upward_stack + downward_stack
    total_flux = np.mean(total_stack, axis=0)
    upward_standard_error = np.std(
        upward_stack, axis=0, ddof=1,
    ) / np.sqrt(len(paths))
    downward_standard_error = np.std(
        downward_stack, axis=0, ddof=1,
    ) / np.sqrt(len(paths))
    total_standard_error = np.std(
        total_stack, axis=0, ddof=1,
    ) / np.sqrt(len(paths))

    np.savez_compressed(
        run_directory / "hot_o_directional_crossing_flux.npz",
        altitude_km=altitude_km,
        energy_eV=energy_eV,
        upward_flux_cm2_s1_per_bin=upward_flux,
        downward_flux_cm2_s1_per_bin=downward_flux,
        total_flux_cm2_s1_per_bin=total_flux,
        net_upward_flux_cm2_s1_per_bin=upward_flux - downward_flux,
        upward_standard_error_cm2_s1_per_bin=upward_standard_error,
        downward_standard_error_cm2_s1_per_bin=downward_standard_error,
        total_standard_error_cm2_s1_per_bin=total_standard_error,
    )
    output_table = np.column_stack(
        [
            np.repeat(altitude_km, energy_eV.size),
            np.tile(energy_eV, altitude_km.size),
            upward_flux.ravel(),
            downward_flux.ravel(),
            total_flux.ravel(),
            (upward_flux - downward_flux).ravel(),
            upward_standard_error.ravel(),
            downward_standard_error.ravel(),
            total_standard_error.ravel(),
        ]
    )
    np.savetxt(
        run_directory / "hot_o_directional_crossing_flux.dat",
        output_table,
        header=(
            "altitude_km energy_eV upward_flux_cm-2_s-1_per_bin "
            "downward_flux_cm-2_s-1_per_bin "
            "total_flux_cm-2_s-1_per_bin "
            "net_upward_flux_cm-2_s-1_per_bin "
            "upward_standard_error_cm-2_s-1_per_bin "
            "downward_standard_error_cm-2_s-1_per_bin "
            "total_standard_error_cm-2_s-1_per_bin"
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

    low_altitude = (altitude_km >= 100.0) & (altitude_km <= 300.0)
    low_altitude_km = altitude_km[low_altitude]
    low_upward_flux = upward_flux[low_altitude]
    low_downward_flux = downward_flux[low_altitude]
    low_positive = np.concatenate(
        [
            low_upward_flux[low_upward_flux > 0],
            low_downward_flux[low_downward_flux > 0],
        ]
    )
    low_vmin = float(np.percentile(np.log10(low_positive), 2))
    low_vmax = float(np.percentile(np.log10(low_positive), 99))
    low_figure, low_axes = plt.subplots(
        1, 2, figsize=(8.2, 3.8), constrained_layout=True, sharey=True,
    )
    low_image = None
    for axis, flux, title, label in zip(
        low_axes,
        (low_upward_flux, low_downward_flux),
        ("Upward, $v_r>0$", "Downward, $v_r<0$"),
        ("a", "b"),
    ):
        positive_flux = flux[flux > 0]
        log_flux = np.log10(np.maximum(flux, positive_flux.min()))
        low_image = axis.imshow(
            log_flux,
            origin="lower",
            aspect="auto",
            extent=(
                energy_eV[0],
                energy_eV[-1],
                low_altitude_km[0],
                low_altitude_km[-1],
            ),
            interpolation="bilinear",
            cmap="turbo",
            vmin=low_vmin,
            vmax=low_vmax,
            rasterized=True,
        )
        axis.set(
            xlabel="Hot O energy (eV)",
            title=title,
            ylim=(100, 300),
            xlim=(0, 6),
        )
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
        )
    low_axes[0].set_ylabel("Altitude (km)")
    low_colorbar = low_figure.colorbar(low_image, ax=low_axes, pad=0.02)
    low_colorbar.set_label(
        r"$\log_{10}[\Phi_k\;(\mathrm{cm^{-2}\,s^{-1}\ per\ bin})]$"
    )
    low_figure.suptitle(
        "Directional hot O flux, 100 to 300 km",
        fontsize=10,
    )
    low_output = (
        run_directory / "upward_downward_hot_o_crossing_flux_100_300km.png"
    )
    low_figure.savefig(
        low_output,
        dpi=400,
        bbox_inches="tight",
    )
    shutil.copy2(
        low_output,
        ROOT / "examples" / "figures" /
        "hot_o_directional_flux_100_300km.png",
    )
    plt.close(low_figure)

    altitude_200_index = int(np.argmin(np.abs(altitude_km - 200.0)))
    if not np.isclose(altitude_km[altitude_200_index], 200.0):
        raise RuntimeError("The directional grid does not contain 200 km")
    upward_200 = upward_flux[altitude_200_index]
    downward_200 = downward_flux[altitude_200_index]
    upward_200_error = upward_standard_error[altitude_200_index]
    downward_200_error = downward_standard_error[altitude_200_index]
    spectrum_figure, spectrum_axis = plt.subplots(
        figsize=(5.2, 4.0), constrained_layout=True,
    )
    upward_valid = upward_200 > 0
    downward_valid = downward_200 > 0
    spectrum_axis.plot(
        energy_eV[upward_valid],
        upward_200[upward_valid],
        color="#2166ac",
        linewidth=1.5,
        label="Upward",
    )
    upward_error_valid = (
        upward_valid & (upward_200_error < 0.8 * upward_200)
    )
    spectrum_axis.fill_between(
        energy_eV[upward_error_valid],
        upward_200[upward_error_valid]
        - upward_200_error[upward_error_valid],
        upward_200[upward_error_valid]
        + upward_200_error[upward_error_valid],
        color="#2166ac",
        alpha=0.20,
        linewidth=0,
    )
    spectrum_axis.plot(
        energy_eV[downward_valid],
        downward_200[downward_valid],
        color="#d95f02",
        linewidth=1.5,
        label="Downward",
    )
    downward_error_valid = (
        downward_valid & (downward_200_error < 0.8 * downward_200)
    )
    spectrum_axis.fill_between(
        energy_eV[downward_error_valid],
        downward_200[downward_error_valid]
        - downward_200_error[downward_error_valid],
        downward_200[downward_error_valid]
        + downward_200_error[downward_error_valid],
        color="#d95f02",
        alpha=0.20,
        linewidth=0,
    )
    spectrum_axis.set(
        yscale="log",
        xlabel="Hot O energy (eV)",
        ylabel=r"Flux (cm$^{-2}$ s$^{-1}$ per bin)",
        title="Directional hot O energy spectra at 200 km",
        xlim=(energy_eV[0], energy_eV[-1]),
    )
    positive_spectrum = np.concatenate(
        [upward_200[upward_valid], downward_200[downward_valid]]
    )
    spectrum_axis.set_ylim(
        positive_spectrum.min() / 2,
        positive_spectrum.max() * 2,
    )
    spectrum_axis.grid(True, color="0.88", linewidth=0.6)
    spectrum_axis.legend(frameon=False)
    spectrum_figure.savefig(
        run_directory / "upward_downward_hot_o_flux_spectrum_200km.png",
        dpi=400,
        bbox_inches="tight",
    )
    plt.close(spectrum_figure)

    altitude_300_index = int(np.argmin(np.abs(altitude_km - 300.0)))
    if not np.isclose(altitude_km[altitude_300_index], 300.0):
        raise RuntimeError("The directional grid does not contain 300 km")
    upward_300 = upward_flux[altitude_300_index]
    downward_300 = downward_flux[altitude_300_index]
    total_300 = total_flux[altitude_300_index]
    upward_300_error = upward_standard_error[altitude_300_index]
    downward_300_error = downward_standard_error[altitude_300_index]
    total_300_error = total_standard_error[altitude_300_index]
    radius_300_m = MARS_RADIUS_M + 300.0e3
    escape_energy_300_eV = (
        MARS_MU_M3_S2 * O_MASS_KG / radius_300_m / EV_J
    )
    escape_energy_bins = energy_eV >= escape_energy_300_eV

    spectrum_300_figure, spectrum_300_axis = plt.subplots(
        figsize=(5.2, 4.0), constrained_layout=True,
    )
    upward_300_valid = upward_300 > 0
    downward_300_valid = downward_300 > 0
    total_300_valid = total_300 > 0
    spectrum_300_axis.plot(
        energy_eV[total_300_valid],
        total_300[total_300_valid],
        color="#222222",
        linewidth=2.0,
        label="Total (upward + downward)",
        zorder=3,
    )
    total_300_error_valid = (
        total_300_valid & (total_300_error < 0.8 * total_300)
    )
    spectrum_300_axis.fill_between(
        energy_eV[total_300_error_valid],
        total_300[total_300_error_valid]
        - total_300_error[total_300_error_valid],
        total_300[total_300_error_valid]
        + total_300_error[total_300_error_valid],
        color="#555555",
        alpha=0.14,
        linewidth=0,
        zorder=1,
    )
    spectrum_300_axis.plot(
        energy_eV[upward_300_valid],
        upward_300[upward_300_valid],
        color="#2166ac",
        linewidth=1.5,
        linestyle="--",
        label="Upward",
        zorder=4,
    )
    upward_300_error_valid = (
        upward_300_valid & (upward_300_error < 0.8 * upward_300)
    )
    spectrum_300_axis.fill_between(
        energy_eV[upward_300_error_valid],
        upward_300[upward_300_error_valid]
        - upward_300_error[upward_300_error_valid],
        upward_300[upward_300_error_valid]
        + upward_300_error[upward_300_error_valid],
        color="#2166ac",
        alpha=0.20,
        linewidth=0,
    )
    spectrum_300_axis.plot(
        energy_eV[downward_300_valid],
        downward_300[downward_300_valid],
        color="#d95f02",
        linewidth=1.5,
        label="Downward",
    )
    downward_300_error_valid = (
        downward_300_valid
        & (downward_300_error < 0.8 * downward_300)
    )
    spectrum_300_axis.fill_between(
        energy_eV[downward_300_error_valid],
        downward_300[downward_300_error_valid]
        - downward_300_error[downward_300_error_valid],
        downward_300[downward_300_error_valid]
        + downward_300_error[downward_300_error_valid],
        color="#d95f02",
        alpha=0.20,
        linewidth=0,
    )
    spectrum_300_axis.axvline(
        escape_energy_300_eV,
        color="0.25",
        linestyle="--",
        linewidth=1.0,
        label=rf"$E_{{esc}}={escape_energy_300_eV:.2f}$ eV",
    )
    spectrum_300_axis.set(
        yscale="log",
        xlabel="Hot O energy (eV)",
        ylabel=r"Flux (cm$^{-2}$ s$^{-1}$ per bin)",
        title="Directional hot O energy spectra at 300 km",
        xlim=(0, 6),
    )
    positive_300_spectrum = np.concatenate(
        [
            upward_300[upward_300_valid],
            downward_300[downward_300_valid],
            total_300[total_300_valid],
        ]
    )
    spectrum_300_axis.set_ylim(
        positive_300_spectrum.min() / 2,
        positive_300_spectrum.max() * 2,
    )
    spectrum_300_axis.grid(True, color="0.88", linewidth=0.6)
    spectrum_300_axis.legend(frameon=False)
    spectrum_300_output = (
        run_directory / "upward_downward_total_hot_o_flux_spectrum_300km"
    )
    for suffix in (".png", ".pdf", ".svg"):
        spectrum_300_figure.savefig(
            spectrum_300_output.with_suffix(suffix),
            dpi=600 if suffix == ".png" else None,
            bbox_inches="tight",
        )
    shutil.copy2(
        spectrum_300_output.with_suffix(".png"),
        ROOT / "examples" / "figures" /
        "hot_o_directional_flux_spectrum_300km.png",
    )
    for suffix in (".pdf", ".svg"):
        shutil.copy2(
            spectrum_300_output.with_suffix(suffix),
            ROOT / "examples" / "figures" /
            f"hot_o_directional_flux_spectrum_300km{suffix}",
        )
    plt.close(spectrum_300_figure)

    upward_300_batch_total = np.sum(
        upward_stack[:, altitude_300_index, :], axis=1
    )
    downward_300_batch_total = np.sum(
        downward_stack[:, altitude_300_index, :], axis=1
    )
    upward_300_batch_escape = np.sum(
        upward_stack[:, altitude_300_index, escape_energy_bins], axis=1
    )
    upward_300_total = float(np.mean(upward_300_batch_total))
    downward_300_total = float(np.mean(downward_300_batch_total))
    upward_300_escape = float(np.mean(upward_300_batch_escape))
    upward_300_total_error = float(
        np.std(upward_300_batch_total, ddof=1) / np.sqrt(len(paths))
    )
    downward_300_total_error = float(
        np.std(downward_300_batch_total, ddof=1) / np.sqrt(len(paths))
    )
    upward_300_escape_error = float(
        np.std(upward_300_batch_escape, ddof=1) / np.sqrt(len(paths))
    )
    projected_area_300_cm2 = np.pi * (radius_300_m * 100.0) ** 2
    spherical_area_300_cm2 = 4.0 * projected_area_300_cm2
    escape_summary = {
        "altitude_km": 300.0,
        "mars_radius_km": MARS_RADIUS_M / 1000.0,
        "local_escape_energy_eV": escape_energy_300_eV,
        "energy_bin_selection": (
            "Energy-bin centers greater than or equal to local escape energy"
        ),
        "upward_flux_all_energy_cm2_s1": upward_300_total,
        "upward_flux_all_energy_standard_error_cm2_s1": (
            upward_300_total_error
        ),
        "downward_flux_all_energy_cm2_s1": downward_300_total,
        "downward_flux_all_energy_standard_error_cm2_s1": (
            downward_300_total_error
        ),
        "upward_escape_capable_flux_cm2_s1": upward_300_escape,
        "upward_escape_capable_flux_standard_error_cm2_s1": (
            upward_300_escape_error
        ),
        "projected_area_pi_r2_cm2": projected_area_300_cm2,
        "all_upward_rate_projected_area_s1": (
            upward_300_total * projected_area_300_cm2
        ),
        "all_upward_rate_projected_area_standard_error_s1": (
            upward_300_total_error * projected_area_300_cm2
        ),
        "escape_rate_projected_area_s1": (
            upward_300_escape * projected_area_300_cm2
        ),
        "escape_rate_projected_area_standard_error_s1": (
            upward_300_escape_error * projected_area_300_cm2
        ),
        "spherical_area_4pi_r2_cm2": spherical_area_300_cm2,
        "all_upward_rate_spherical_area_s1": (
            upward_300_total * spherical_area_300_cm2
        ),
        "all_upward_rate_spherical_area_standard_error_s1": (
            upward_300_total_error * spherical_area_300_cm2
        ),
        "escape_rate_spherical_area_s1": (
            upward_300_escape * spherical_area_300_cm2
        ),
        "escape_rate_spherical_area_standard_error_s1": (
            upward_300_escape_error * spherical_area_300_cm2
        ),
    }
    (run_directory / "escape_flux_300km.json").write_text(
        json.dumps(escape_summary, indent=2), encoding="utf-8"
    )

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
    print(json.dumps(escape_summary, indent=2))


if __name__ == "__main__":
    main()
