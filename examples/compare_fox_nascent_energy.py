"""Compare the Fox hot O energy grid with the MarsHotO MGITM result."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FOX_GRID = Path(
    r"D:\Work_Work\Mars\batsrus_mhd\SWMF\PT\AMPS\srcMars\EnergyGridFox.h"
)
DEFAULT_MHOT_CSV = (
    ROOT
    / "examples"
    / "figures"
    / "mgitm_ls000_f070_hot_o_nascent_energy_with_vibration_map.csv"
)
DEFAULT_OUTPUT = (
    ROOT / "examples" / "figures" / "fox_mhot_nascent_energy_comparison"
)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_fox_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the 65 by 300 probability-mass table from EnergyGridFox.h."""
    text = path.read_text(encoding="utf-8")
    number = r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?"
    pattern = re.compile(
        rf"\{{\s*({number})\s*,\s*({number})\s*,\s*({number})\s*,"
        rf"\s*\{{(.*?)\}}\s*,\s*({number})\s*,\s*({number})\s*,"
        rf"\s*({number})\s*,\s*({number})\s*\}}",
        re.DOTALL,
    )
    altitude_km = []
    probability_mass = []
    energy_eV = None
    for match in pattern.finditer(text):
        values = np.array(
            [float(value) for value in re.findall(number, match.group(4))]
        )
        if values.size != 300:
            raise ValueError(f"Expected 300 Fox energy values, found {values.size}")
        emin_eV = float(match.group(5))
        step_eV = float(match.group(7))
        current_energy = emin_eV + step_eV * np.arange(values.size)
        if energy_eV is None:
            energy_eV = current_energy
        elif not np.allclose(energy_eV, current_energy):
            raise ValueError("Fox energy grids are not identical between altitudes")
        altitude_km.append(float(match.group(1)))
        probability_mass.append(values / values.sum())
    if len(altitude_km) != 65:
        raise ValueError(f"Expected 65 Fox altitude levels, found {len(altitude_km)}")
    return (
        np.asarray(altitude_km),
        np.asarray(energy_eV),
        np.asarray(probability_mass),
    )


def load_mhot_grid(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = pd.read_csv(path)
    altitude_km = np.sort(table["altitude_km"].unique())
    energy_eV = np.sort(table["energy_eV"].unique())
    probability_mass = (
        table.pivot(
            index="altitude_km",
            columns="energy_eV",
            values="energy_bin_probability",
        )
        .loc[altitude_km, energy_eV]
        .to_numpy()
    )
    probability_mass /= probability_mass.sum(axis=1, keepdims=True)
    return altitude_km, energy_eV, probability_mass


def interpolate_fox_density(
    fox_altitude_km: np.ndarray,
    fox_energy_eV: np.ndarray,
    fox_density_eV1: np.ndarray,
    target_altitude_km: np.ndarray,
    target_energy_eV: np.ndarray,
) -> np.ndarray:
    energy_interpolated = np.vstack(
        [
            np.interp(target_energy_eV, fox_energy_eV, row, left=0.0, right=0.0)
            for row in fox_density_eV1
        ]
    )
    return np.vstack(
        [
            np.interp(
                target_altitude_km,
                fox_altitude_km,
                energy_interpolated[:, index],
            )
            for index in range(target_energy_eV.size)
        ]
    ).T


def interpolate_altitude_profile(
    altitude_km: np.ndarray,
    density_eV1: np.ndarray,
    target_altitude_km: float,
) -> np.ndarray:
    """Linearly interpolate every energy channel to one target altitude."""
    return np.array(
        [
            np.interp(target_altitude_km, altitude_km, density_eV1[:, index])
            for index in range(density_eV1.shape[1])
        ]
    )


def probability_mean(energy_eV: np.ndarray, density_eV1: np.ndarray) -> float:
    return float(np.trapz(energy_eV * density_eV1, energy_eV))


def make_figure(
    fox_altitude_km: np.ndarray,
    fox_energy_eV: np.ndarray,
    fox_density_eV1: np.ndarray,
    mhot_altitude_km: np.ndarray,
    mhot_energy_eV: np.ndarray,
    mhot_density_eV1: np.ndarray,
) -> plt.Figure:
    altitude_mask_fox = (fox_altitude_km >= 100.0) & (fox_altitude_km <= 250.0)
    energy_mask_fox = fox_energy_eV <= 5.0
    altitude_mask_mhot = (mhot_altitude_km >= 100.0) & (mhot_altitude_km <= 250.0)
    energy_mask_mhot = mhot_energy_eV <= 5.0

    fa = fox_altitude_km[altitude_mask_fox]
    fe = fox_energy_eV[energy_mask_fox]
    fd = fox_density_eV1[np.ix_(altitude_mask_fox, energy_mask_fox)]
    ma = mhot_altitude_km[altitude_mask_mhot]
    me = mhot_energy_eV[energy_mask_mhot]
    md = mhot_density_eV1[np.ix_(altitude_mask_mhot, energy_mask_mhot)]
    fox_on_mhot = interpolate_fox_density(fa, fe, fd, ma, me)
    difference = md - fox_on_mhot

    common_max = float(np.nanpercentile(np.concatenate((fd.ravel(), md.ravel())), 99.5))
    difference_max = float(np.nanpercentile(np.abs(difference), 99.5))

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    probability_cmap = "viridis"
    difference_cmap = "RdBu_r"

    mesh_a = axes[0, 0].pcolormesh(
        fe, fa, fd, shading="nearest", cmap=probability_cmap, vmin=0.0, vmax=common_max
    )
    mesh_b = axes[0, 1].pcolormesh(
        me, ma, md, shading="nearest", cmap=probability_cmap, vmin=0.0, vmax=common_max
    )
    mesh_c = axes[1, 0].pcolormesh(
        me,
        ma,
        difference,
        shading="nearest",
        cmap=difference_cmap,
        vmin=-difference_max,
        vmax=difference_max,
    )

    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xlim(0.0, 5.0)
        axis.set_ylim(100.0, 250.0)
        axis.set_xlabel("Nascent O energy (eV)")
        axis.set_ylabel("Altitude (km)")

    axes[0, 0].set_title("Fox energy grid")
    axes[0, 1].set_title("MarsHotO, MGITM $L_s=0^\circ$, F070")
    axes[1, 0].set_title("MarsHotO minus Fox")
    shared_colorbar = fig.colorbar(mesh_b, ax=axes[0, :], shrink=0.91, pad=0.02)
    shared_colorbar.set_label("Conditional probability density (eV$^{-1}$)")
    difference_colorbar = fig.colorbar(mesh_c, ax=axes[1, 0], shrink=0.91, pad=0.02)
    difference_colorbar.set_label(r"$\Delta p(E\mid z)$ (eV$^{-1}$)")

    profile_axis = axes[1, 1]
    target_altitude = 140.0
    fox_profile = interpolate_altitude_profile(fa, fd, target_altitude)
    mhot_profile = interpolate_altitude_profile(ma, md, target_altitude)
    profile_axis.plot(
        fe, fox_profile, color="#0072B2", linestyle="--", linewidth=1.5,
        label="Fox",
    )
    profile_axis.plot(
        me, mhot_profile, color="#D55E00", linewidth=1.7, label="MarsHotO",
    )
    mean_fox = probability_mean(fe, fox_profile)
    mean_mhot = probability_mean(me, mhot_profile)
    print(
        f"{target_altitude:.0f} km interpolated mean energies: "
        f"Fox={mean_fox:.4f} eV, MarsHotO={mean_mhot:.4f} eV"
    )
    profile_axis.set_xlim(0.0, 5.0)
    profile_axis.set_ylim(bottom=0.0)
    profile_axis.set_xlabel("Nascent O energy (eV)")
    profile_axis.set_ylabel("Probability density (eV$^{-1}$)")
    profile_axis.set_title("Energy distribution at 140 km")
    profile_axis.legend(loc="upper right")

    for label, axis in zip("abcd", axes.ravel()):
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            fontsize=10,
            fontweight="bold",
            ha="left",
            va="top",
            color="black",
        )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fox-grid", type=Path, default=DEFAULT_FOX_GRID)
    parser.add_argument("--mhot-csv", type=Path, default=DEFAULT_MHOT_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    configure_matplotlib()
    fox_altitude, fox_energy, fox_mass = load_fox_grid(args.fox_grid)
    mhot_altitude, mhot_energy, mhot_mass = load_mhot_grid(args.mhot_csv)
    fox_step_eV = float(np.median(np.diff(fox_energy)))
    mhot_step_eV = float(np.median(np.diff(mhot_energy)))
    fox_density = fox_mass / fox_step_eV
    mhot_density = mhot_mass / mhot_step_eV

    figure = make_figure(
        fox_altitude,
        fox_energy,
        fox_density,
        mhot_altitude,
        mhot_energy,
        mhot_density,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output.with_suffix(".png"), dpi=400, bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {args.output.with_suffix('.png')}")


if __name__ == "__main__":
    main()
