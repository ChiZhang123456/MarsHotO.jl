"""Plot discrete hot O reaction branches without thermal broadening."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

from plot_mgitm_hot_o_profiles import (
    INPUT_FILE,
    calculate_derived_profiles,
    load_nearest_subsolar_profile,
)


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "examples" / "figures"
OUTPUT_STEM = OUTPUT_DIR / "mgitm_ls000_f070_discrete_hot_o_energy_lines"
SOURCE_DATA_FILE = OUTPUT_DIR / "mgitm_ls000_f070_discrete_hot_o_energy_lines.csv"

# Single O energies are one half of the total reaction exothermicities.
ENERGY_EV = np.array([0.415, 1.530, 2.510, 3.495])
BRANCH_PROBABILITY = np.array([0.058, 0.204, 0.473, 0.265])


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def altitude_edges(altitude_centers_km: np.ndarray) -> np.ndarray:
    midpoint = 0.5 * (altitude_centers_km[:-1] + altitude_centers_km[1:])
    first = altitude_centers_km[0] - (
        midpoint[0] - altitude_centers_km[0]
    )
    last = altitude_centers_km[-1] + (
        altitude_centers_km[-1] - midpoint[-1]
    )
    return np.concatenate(([first], midpoint, [last]))


def colored_vertical_line(
    axis: plt.Axes,
    energy_ev: float,
    altitude_edges_km: np.ndarray,
    values: np.ndarray,
    cmap: str,
    norm: Normalize,
    linewidth: float = 7.0,
) -> None:
    segments = np.stack(
        (
            np.column_stack(
                (
                    np.full(len(values), energy_ev),
                    altitude_edges_km[:-1],
                )
            ),
            np.column_stack(
                (
                    np.full(len(values), energy_ev),
                    altitude_edges_km[1:],
                )
            ),
        ),
        axis=1,
    )
    collection = LineCollection(
        segments,
        cmap=cmap,
        norm=norm,
        linewidth=linewidth,
        capstyle="butt",
    )
    collection.set_array(values)
    axis.add_collection(collection)


def prepare_source_data(profile: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for energy_ev, probability in zip(ENERGY_EV, BRANCH_PROBABILITY):
        branch_production = probability * profile["Q_hotO_cm3s"].to_numpy()
        for altitude_km, total_production, production in zip(
            profile["altitude_km"],
            profile["Q_hotO_cm3s"],
            branch_production,
        ):
            rows.append(
                {
                    "altitude_km": altitude_km,
                    "energy_eV": energy_ev,
                    "branch_probability": probability,
                    "total_hot_O_production_cm-3_s-1": total_production,
                    "branch_hot_O_production_cm-3_s-1": production,
                }
            )
    return pd.DataFrame(rows)


def make_figure(profile: pd.DataFrame) -> plt.Figure:
    configure_matplotlib()
    altitude_km = profile["altitude_km"].to_numpy()
    edges_km = altitude_edges(altitude_km)

    probability_norm = Normalize(vmin=0.0, vmax=0.50)
    branch_production = (
        BRANCH_PROBABILITY[:, None]
        * profile["Q_hotO_cm3s"].to_numpy()[None, :]
    )
    log_branch_production = np.log10(branch_production)
    production_norm = Normalize(
        vmin=float(np.floor(np.min(log_branch_production))),
        vmax=float(np.ceil(np.max(log_branch_production))),
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(7.2, 4.1),
        sharey=True,
        constrained_layout=True,
    )

    for energy_ev, probability in zip(ENERGY_EV, BRANCH_PROBABILITY):
        colored_vertical_line(
            axes[0],
            energy_ev,
            edges_km,
            np.full(len(altitude_km), probability),
            cmap="viridis",
            norm=probability_norm,
        )
        axes[0].text(
            energy_ev,
            edges_km[-1] + 1.5,
            f"{100.0 * probability:.1f}%",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    for branch_index, energy_ev in enumerate(ENERGY_EV):
        colored_vertical_line(
            axes[1],
            energy_ev,
            edges_km,
            log_branch_production[branch_index],
            cmap="turbo",
            norm=production_norm,
        )

    axes[0].set_title("Fixed branch probability")
    axes[1].set_title("Branch production rate")
    axes[0].set_ylabel("Altitude (km)")

    for axis in axes:
        axis.set_xlabel("Nascent O energy (eV)")
        axis.set_xlim(0.0, 4.0)
        axis.set_ylim(edges_km[0], edges_km[-1] + 8.0)
        axis.set_xticks([0, 1, 2, 3, 4])
        axis.grid(axis="y", color="0.90", lw=0.6)
        axis.set_facecolor("#F5F5F5")

    probability_mappable = mpl.cm.ScalarMappable(
        norm=probability_norm,
        cmap="viridis",
    )
    probability_colorbar = fig.colorbar(
        probability_mappable,
        ax=axes[0],
        pad=0.03,
        fraction=0.05,
    )
    probability_colorbar.set_label("Branch probability")

    production_mappable = mpl.cm.ScalarMappable(
        norm=production_norm,
        cmap="turbo",
    )
    production_colorbar = fig.colorbar(
        production_mappable,
        ax=axes[1],
        pad=0.03,
        fraction=0.05,
    )
    production_colorbar.set_label(
        r"$\log_{10}[Q_i(z)\;(\mathrm{cm^{-3}\,s^{-1}})]$"
    )

    axes[0].text(
        0.02,
        0.98,
        "a",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontweight="bold",
        fontsize=9,
    )
    axes[1].text(
        0.02,
        0.98,
        "b",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontweight="bold",
        fontsize=9,
    )
    fig.suptitle(
        r"Discrete hot O source channels without thermal broadening, "
        r"$L_s=0^\circ$, F070",
        fontsize=10,
    )
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    profile, metadata = load_nearest_subsolar_profile(INPUT_FILE)
    profile = calculate_derived_profiles(profile)
    source_data = prepare_source_data(profile)
    source_data.to_csv(SOURCE_DATA_FILE, index=False)

    probability_sum = float(BRANCH_PROBABILITY.sum())
    grouped_sum = (
        source_data.groupby("altitude_km")[
            "branch_hot_O_production_cm-3_s-1"
        ].sum()
    )
    total_production = (
        source_data.groupby("altitude_km")[
            "total_hot_O_production_cm-3_s-1"
        ].first()
    )
    maximum_relative_error = float(
        np.max(np.abs(grouped_sum - total_production) / total_production)
    )

    fig = make_figure(profile)
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    print(f"Input: {INPUT_FILE}")
    print(
        "Selected columns: "
        f"lon={metadata['selected_longitudes_deg']}, "
        f"lat={metadata['selected_latitudes_deg']}, "
        f"SZA={metadata['minimum_grid_sza_deg']:.4f} deg"
    )
    print(f"Branch probability sum: {probability_sum:.12f}")
    print(
        "Maximum branch-sum relative error: "
        f"{maximum_relative_error:.3e}"
    )
    print(f"Output: {OUTPUT_STEM.with_suffix('.png')}")
    print(f"Source data: {SOURCE_DATA_FILE}")


if __name__ == "__main__":
    main()
