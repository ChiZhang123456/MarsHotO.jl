"""Plot effective differential cross sections used by the AMPS Mars setup.

The active Mars configuration loads one O-O angular profile at 3 eV. The
same scattering-angle sampler is called for every background species. O-CO2
and O-N2 use constant total cross sections, so their effective differential
cross sections below are the normalized O-O angular profile scaled to those
totals. Consequently, all three maps are energy independent.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_AMPS_TABLE = Path(
    r"D:\Work_Work\Mars\batsrus_mhd\SWMF\PT\AMPS\srcMars"
    r"\Kharchenko-2000-fig-3.h"
)
DEFAULT_OUTPUT = (
    ROOT / "examples" / "figures" / "amps_active_differential_cross_sections"
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
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def load_active_fox_profile(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read the active _FOX_DATA_ angle and d-sigma/d-Omega arrays."""
    text = path.read_text(encoding="utf-8")
    block_match = re.search(
        r"_FOX_DATA_\[_FOX_LENGTH_\]\s*=\s*\{(.*?)\};", text, re.DOTALL
    )
    if block_match is None:
        raise ValueError("Could not find _FOX_DATA_ in the AMPS table")
    pairs = re.findall(
        r"\{\s*([0-9.eE+-]+)\s*,\s*([0-9.eE+-]+)\s*\}",
        block_match.group(1),
    )
    if len(pairs) != 137:
        raise ValueError(f"Expected 137 active Fox points, found {len(pairs)}")
    values = np.asarray([(float(a), float(b)) for a, b in pairs])
    return values[:, 0], values[:, 1]


def log_interpolate_profile(
    angle_data_deg: np.ndarray,
    dcs_data_m2_sr: np.ndarray,
    angle_grid_deg: np.ndarray,
) -> np.ndarray:
    return np.exp(
        np.interp(angle_grid_deg, angle_data_deg, np.log(dcs_data_m2_sr))
    )


def integrated_total_cross_section(
    angle_deg: np.ndarray, dcs_m2_sr: np.ndarray
) -> float:
    angle_rad = np.deg2rad(angle_deg)
    return float(2.0 * np.pi * np.trapz(dcs_m2_sr * np.sin(angle_rad), angle_rad))


def make_figure(
    energy_eV: np.ndarray,
    angle_deg: np.ndarray,
    angular_shape_m2_sr: np.ndarray,
    source_total_m2: float,
) -> plt.Figure:
    species = (
        (r"O-CO$_2$", 2.0e-18),
        (r"O-O", source_total_m2),
        (r"O-N$_2$", 1.8e-18),
    )
    maps = []
    for _, total_m2 in species:
        profile = angular_shape_m2_sr * total_m2 / source_total_m2
        maps.append(np.repeat(profile[:, None], energy_eV.size, axis=1))

    fig, axes = plt.subplots(
        1, 3, figsize=(9.0, 2.75), sharex=True, sharey=True,
        constrained_layout=True,
    )
    for label, total_m2, axis, values in zip(
        (item[0] for item in species),
        (item[1] for item in species),
        axes,
        maps,
    ):
        norm = LogNorm(vmin=float(values.min()), vmax=float(values.max()))
        mesh = axis.pcolormesh(
            energy_eV,
            angle_deg,
            values,
            shading="nearest",
            cmap="viridis",
            norm=norm,
        )
        axis.set_xscale("log")
        axis.set_xlim(0.01, 10.0)
        axis.set_ylim(angle_deg.min(), angle_deg.max())
        axis.set_xlabel("Collision energy (eV)")
        axis.set_title(
            f"{label}\n$\\sigma_{{\\mathrm{{tot}}}}={total_m2:.2g}$ m$^2$"
        )
        axis.text(
            0.96,
            0.96,
            "energy independent",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=7,
            color="white",
        )
        colorbar = fig.colorbar(mesh, ax=axis, shrink=0.92, pad=0.02)
        colorbar.set_label(r"$d\sigma/d\Omega$ (m$^2$ sr$^{-1}$)")
    axes[0].set_ylabel(r"COM scattering angle, $\theta$ (deg)")
    for label, axis in zip("abc", axes):
        axis.text(
            0.03,
            0.96,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            fontweight="bold",
            color="white",
        )
    return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps-table", type=Path, default=DEFAULT_AMPS_TABLE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    configure_matplotlib()
    angle_data_deg, dcs_data_m2_sr = load_active_fox_profile(args.amps_table)
    integration_angle_deg = np.linspace(
        angle_data_deg.min(), angle_data_deg.max(), 200_000
    )
    integration_dcs_m2_sr = log_interpolate_profile(
        angle_data_deg, dcs_data_m2_sr, integration_angle_deg
    )
    source_total_m2 = integrated_total_cross_section(
        integration_angle_deg, integration_dcs_m2_sr
    )
    angle_deg = np.linspace(angle_data_deg.min(), angle_data_deg.max(), 600)
    dcs_m2_sr = log_interpolate_profile(
        angle_data_deg, dcs_data_m2_sr, angle_deg
    )
    energy_eV = np.geomspace(0.01, 10.0, 400)
    figure = make_figure(energy_eV, angle_deg, dcs_m2_sr, source_total_m2)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output.with_suffix(".png"), dpi=400, bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(args.output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(figure)
    print(f"Integrated active O-O profile: {source_total_m2:.8e} m^2")
    print(f"Saved {args.output.with_suffix('.png')}")


if __name__ == "__main__":
    main()
