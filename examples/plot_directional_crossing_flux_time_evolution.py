"""Cumulative steady-state crossing-flux contribution versus particle age."""

from pathlib import Path
import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import numpy as np

from plot_directional_hot_o_flux import (
    EVENT_CROSSING,
    EVENT_DTYPE,
    EV_J,
    HEADER,
    MARS_RADIUS_M,
    O_MASS_KG,
    read_header,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN = ROOT / "examples/output/run_paired_dr_10000_crossings"
TIMES_S = np.arange(0.0, 100.0 + 5.0, 5.0)
ALTITUDES_KM = np.arange(100.0, 300.0 + 10.0, 10.0)
ENERGY_EDGES_EV = np.arange(0.0, 6.0 + 0.05, 0.05)
ENERGY_EV = (ENERGY_EDGES_EV[:-1] + ENERGY_EDGES_EV[1:]) / 2


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
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


def cumulative_flux_for_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    header = read_header(path)
    records = np.memmap(
        path, mode="r", dtype=EVENT_DTYPE, offset=HEADER.size,
        shape=(header["events"],),
    )
    shape = (TIMES_S.size, ALTITUDES_KM.size, ENERGY_EV.size)
    upward_increment = np.zeros(shape)
    downward_increment = np.zeros(shape)
    chunk_size = 2_000_000
    for start in range(0, header["events"], chunk_size):
        event = records[start : start + chunk_size]
        selected = (
            (event["event_code"] == EVENT_CROSSING)
            & (event["surface_index"] >= 1)
            & (event["surface_index"] <= ALTITUDES_KM.size)
            & (event["time_s"] <= TIMES_S[-1])
        )
        event = event[selected]
        if event.size == 0:
            continue
        time_bin = np.searchsorted(TIMES_S, event["time_s"], side="left")
        surface = event["surface_index"].astype(np.int64) - 1
        speed_squared = (
            event["velocity_x_m_s"] ** 2
            + event["velocity_y_m_s"] ** 2
            + event["velocity_z_m_s"] ** 2
        )
        energy = 0.5 * O_MASS_KG * speed_squared / EV_J
        energy_bin = np.searchsorted(
            ENERGY_EDGES_EV, energy, side="right"
        ) - 1
        valid = (energy_bin >= 0) & (energy_bin < ENERGY_EV.size)
        time_bin = time_bin[valid]
        surface = surface[valid]
        energy_bin = energy_bin[valid]
        weight = event["weight_s1"][valid]
        upward = event["direction"][valid] > 0
        np.add.at(
            upward_increment,
            (time_bin[upward], surface[upward], energy_bin[upward]),
            weight[upward],
        )
        np.add.at(
            downward_increment,
            (time_bin[~upward], surface[~upward], energy_bin[~upward]),
            weight[~upward],
        )
    area_m2 = 4 * np.pi * (
        MARS_RADIUS_M + 1000 * ALTITUDES_KM
    ) ** 2
    normalization = area_m2[None, :, None] * 1e4
    upward = np.cumsum(upward_increment, axis=0) / normalization
    downward = np.cumsum(downward_increment, axis=0) / normalization
    return upward, downward


def load_batches(paths: list[Path]) -> tuple[np.ndarray, ...]:
    upward_batches = []
    downward_batches = []
    for index, path in enumerate(paths, start=1):
        upward, downward = cumulative_flux_for_file(path)
        upward_batches.append(upward)
        downward_batches.append(downward)
        print(f"processed {index}/{len(paths)}: {path.name}", flush=True)
    upward_stack = np.stack(upward_batches)
    downward_stack = np.stack(downward_batches)
    upward_mean = upward_stack.mean(axis=0)
    downward_mean = downward_stack.mean(axis=0)
    upward_se = upward_stack.std(axis=0, ddof=1) / np.sqrt(len(paths))
    downward_se = downward_stack.std(axis=0, ddof=1) / np.sqrt(len(paths))
    return upward_mean, downward_mean, upward_se, downward_se


def color_limits(upward: np.ndarray, downward: np.ndarray) -> tuple[float, float]:
    positive = np.concatenate((upward[upward > 0], downward[downward > 0]))
    return tuple(np.percentile(np.log10(positive), [1.0, 99.5]))


def log_values(values: np.ndarray, vmin: float) -> np.ndarray:
    return np.log10(np.maximum(values, 10**vmin))


def draw_static(
    upward: np.ndarray, downward: np.ndarray, output_base: Path,
) -> None:
    selected_times = (0.0, 50.0, 100.0)
    indices = [int(np.where(TIMES_S == value)[0][0]) for value in selected_times]
    vmin, vmax = color_limits(upward, downward)
    figure, axes = plt.subplots(
        3, 2, figsize=(7.2, 8.4), sharex=True, sharey=True,
        constrained_layout=True,
    )
    labels = iter("abcdef")
    image = None
    extent = (0, 6, 100, 300)
    for row, (time_s, time_index) in enumerate(zip(selected_times, indices)):
        for column, (values, direction) in enumerate(
            ((upward[time_index], "Upward"),
             (downward[time_index], "Downward"))
        ):
            axis = axes[row, column]
            image = axis.imshow(
                log_values(values, vmin), origin="lower", aspect="auto",
                extent=extent, interpolation="bilinear", cmap="magma",
                vmin=vmin, vmax=vmax, rasterized=True,
            )
            axis.set_title(
                rf"{direction}, $\tau\leq {time_s:g}$ s"
            )
            axis.text(
                0.02, 0.97, next(labels), transform=axis.transAxes,
                ha="left", va="top", color="white", fontweight="bold",
                fontsize=9,
            )
            if column == 0:
                axis.set_ylabel("Altitude (km)")
            if row == 2:
                axis.set_xlabel("Hot O energy (eV)")
    colorbar = figure.colorbar(image, ax=axes, pad=0.02, shrink=0.97)
    colorbar.set_label(
        r"$\log_{10}[\Phi_k(E,z;\tau\leq t)]$ "
        r"(cm$^{-2}$ s$^{-1}$ per bin)"
    )
    figure.suptitle(
        "Cumulative contribution to steady-state directional hot O flux",
        fontsize=11,
    )
    figure.savefig(output_base.with_suffix(".png"), dpi=400)
    figure.savefig(output_base.with_suffix(".pdf"))
    figure.savefig(output_base.with_suffix(".svg"))
    plt.close(figure)


def draw_animation(
    upward: np.ndarray, downward: np.ndarray, output: Path,
) -> None:
    vmin, vmax = color_limits(upward, downward)
    figure, axes = plt.subplots(
        1, 2, figsize=(9.5, 4.4), sharex=True, sharey=True,
        constrained_layout=True,
    )
    images = []
    for axis, values, direction in zip(
        axes, (upward, downward), ("Upward", "Downward")
    ):
        image = axis.imshow(
            log_values(values[0], vmin), origin="lower", aspect="auto",
            extent=(0, 6, 100, 300), interpolation="bilinear",
            cmap="magma", vmin=vmin, vmax=vmax, animated=True,
        )
        axis.set(
            xlabel="Hot O energy (eV)", ylabel="Altitude (km)",
            title=direction,
        )
        images.append(image)
    colorbar = figure.colorbar(images[0], ax=axes, pad=0.02)
    colorbar.set_label(
        r"$\log_{10}[\Phi_k(E,z;\tau\leq t)]$ "
        r"(cm$^{-2}$ s$^{-1}$ per bin)"
    )
    title = figure.suptitle("")

    def update(frame: int):
        images[0].set_data(log_values(upward[frame], vmin))
        images[1].set_data(log_values(downward[frame], vmin))
        title.set_text(
            "Cumulative contribution to steady-state directional flux\n"
            rf"particle flight age $\tau\leq {TIMES_S[frame]:g}$ s"
        )
        return (*images, title)

    animation = FuncAnimation(
        figure, update, frames=TIMES_S.size, interval=350, blit=False,
    )
    animation.save(output, writer=PillowWriter(fps=3), dpi=150)
    plt.close(figure)


def write_source_data(
    output_dir: Path, upward: np.ndarray, downward: np.ndarray,
    upward_se: np.ndarray, downward_se: np.ndarray,
) -> None:
    np.savez_compressed(
        output_dir / "hot_o_directional_crossing_flux_time_evolution.npz",
        time_s=TIMES_S, altitude_km=ALTITUDES_KM, energy_eV=ENERGY_EV,
        upward_flux_cm2_s1_per_bin=upward,
        downward_flux_cm2_s1_per_bin=downward,
        upward_standard_error_cm2_s1_per_bin=upward_se,
        downward_standard_error_cm2_s1_per_bin=downward_se,
    )
    rows = []
    for it, time_s in enumerate(TIMES_S):
        for ia, altitude in enumerate(ALTITUDES_KM):
            for ie, energy in enumerate(ENERGY_EV):
                rows.append((
                    time_s, altitude, energy,
                    upward[it, ia, ie], downward[it, ia, ie],
                    upward_se[it, ia, ie], downward_se[it, ia, ie],
                ))
    np.savetxt(
        output_dir / "hot_o_directional_crossing_flux_time_evolution.dat",
        np.asarray(rows),
        header=(
            "maximum_flight_age_s altitude_km energy_eV "
            "upward_flux_cm-2_s-1_per_bin "
            "downward_flux_cm-2_s-1_per_bin upward_standard_error "
            "downward_standard_error"
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", nargs="?", type=Path, default=DEFAULT_RUN)
    args = parser.parse_args()
    run_directory = args.run_directory.resolve()
    paths = sorted(run_directory.glob("batch_*.bin"))
    if not paths:
        raise RuntimeError(f"No crossing-event files in {run_directory}")
    configure_matplotlib()
    upward, downward, upward_se, downward_se = load_batches(paths)
    write_source_data(run_directory, upward, downward, upward_se, downward_se)
    figure_dir = ROOT / "examples/figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    draw_static(
        upward, downward,
        figure_dir / "hot_o_directional_crossing_flux_time_panels",
    )
    draw_animation(
        upward, downward,
        figure_dir / "hot_o_directional_crossing_flux_time_evolution.gif",
    )
    print("completed directional crossing-flux time evolution")


if __name__ == "__main__":
    main()
