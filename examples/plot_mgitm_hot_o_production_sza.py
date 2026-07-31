"""Plot MGITM hot O production versus solar zenith angle and altitude."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_mgitm_hot_o_profiles import (  # noqa: E402
    INPUT_FILE,
    MGITM_COLUMNS,
    angular_distance_deg,
    dissociative_recombination_coefficient_cm3s,
    read_subsolar_longitude,
)


OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "mgitm_ls000_f070_hot_o_production_sza_altitude.png"
)
SUBSOLAR_LATITUDE_DEG = 0.0
SZA_BIN_WIDTH_DEG = 3.0
LOG10_Q_MIN = -4.0
LOG10_Q_MAX = 4.5
TARGET_SZA_DEG = np.array([2.0, 10.0, 30.0, 60.0, 90.0, 120.0])
TARGET_SZA_HALF_WIDTH_DEG = 1.5


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


def calculate_sza_deg(
    longitude_deg: np.ndarray,
    latitude_deg: np.ndarray,
    subsolar_longitude_deg: float,
    subsolar_latitude_deg: float = SUBSOLAR_LATITUDE_DEG,
) -> np.ndarray:
    longitude_difference_rad = np.deg2rad(
        angular_distance_deg(
            np.asarray(longitude_deg), subsolar_longitude_deg
        )
    )
    latitude_rad = np.deg2rad(np.asarray(latitude_deg))
    subsolar_latitude_rad = np.deg2rad(subsolar_latitude_deg)
    cos_sza = (
        np.sin(latitude_rad) * np.sin(subsolar_latitude_rad)
        + np.cos(latitude_rad)
        * np.cos(subsolar_latitude_rad)
        * np.cos(longitude_difference_rad)
    )
    return np.rad2deg(np.arccos(np.clip(cos_sza, -1.0, 1.0)))


def altitude_edges(altitude_km: np.ndarray) -> np.ndarray:
    midpoint = 0.5 * (altitude_km[:-1] + altitude_km[1:])
    return np.concatenate(
        (
            [altitude_km[0] - (midpoint[0] - altitude_km[0])],
            midpoint,
            [altitude_km[-1] + (altitude_km[-1] - midpoint[-1])],
        )
    )


def calculate_production_grid(
    path: Path,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    pd.DataFrame,
]:
    subsolar_longitude_deg = read_subsolar_longitude(path)
    data = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=MGITM_COLUMNS,
        engine="c",
    )
    data["sza_deg"] = calculate_sza_deg(
        data["longitude_deg"].to_numpy(),
        data["latitude_deg"].to_numpy(),
        subsolar_longitude_deg,
    )
    electron_density_cm3 = data["ne_m3"].to_numpy() * 1.0e-6
    o2plus_density_cm3 = data["nO2p_m3"].to_numpy() * 1.0e-6
    coefficient_cm3_s1 = dissociative_recombination_coefficient_cm3s(
        data["Te_K"].to_numpy()
    )
    data["Q_hotO_cm3_s1"] = (
        2.0
        * electron_density_cm3
        * o2plus_density_cm3
        * coefficient_cm3_s1
    )

    sza_edges_deg = np.arange(
        0.0, 180.0 + SZA_BIN_WIDTH_DEG, SZA_BIN_WIDTH_DEG
    )
    data["sza_bin"] = pd.cut(
        data["sza_deg"],
        bins=sza_edges_deg,
        labels=False,
        include_lowest=True,
        right=False,
    )
    data["area_weight"] = np.cos(
        np.deg2rad(data["latitude_deg"].to_numpy())
    )
    data["weighted_Q"] = (
        data["Q_hotO_cm3_s1"] * data["area_weight"]
    )

    grouped = data.groupby(
        ["altitude_km", "sza_bin"], observed=True
    ).agg(
        weighted_Q_sum=("weighted_Q", "sum"),
        area_weight_sum=("area_weight", "sum"),
        column_count=("longitude_deg", "size"),
    )
    grouped["Q_hotO_cm3_s1"] = (
        grouped["weighted_Q_sum"] / grouped["area_weight_sum"]
    )

    altitude_km = np.sort(data["altitude_km"].unique())
    sza_bin_index = np.arange(sza_edges_deg.size - 1)
    production = (
        grouped["Q_hotO_cm3_s1"]
        .unstack("sza_bin")
        .reindex(index=altitude_km, columns=sza_bin_index)
        .to_numpy()
    )
    if np.any(~np.isfinite(production)):
        raise RuntimeError("The SZA and altitude production grid is incomplete")
    if np.any(production <= 0):
        raise RuntimeError("Hot O production rates must be positive")

    column_counts = (
        data[["longitude_deg", "latitude_deg", "sza_bin"]]
        .drop_duplicates()
        .groupby("sza_bin", observed=True)
        .size()
        .reindex(sza_bin_index)
        .to_numpy()
    )
    return (
        altitude_km,
        sza_edges_deg,
        production,
        column_counts,
        data,
    )


def calculate_selected_sza_profiles(
    data: pd.DataFrame,
    altitude_km: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    profiles = []
    column_counts = []
    for target_sza_deg in TARGET_SZA_DEG:
        selected = data.loc[
            np.abs(data["sza_deg"] - target_sza_deg)
            <= TARGET_SZA_HALF_WIDTH_DEG
        ].copy()
        selected_columns = selected[
            ["longitude_deg", "latitude_deg"]
        ].drop_duplicates()
        if selected_columns.empty:
            raise RuntimeError(
                f"No MGITM columns near SZA={target_sza_deg:g} deg"
            )
        profile = (
            selected.groupby("altitude_km", observed=True)
            .agg(
                weighted_Q_sum=("weighted_Q", "sum"),
                area_weight_sum=("area_weight", "sum"),
            )
            .reindex(altitude_km)
        )
        production = (
            profile["weighted_Q_sum"] / profile["area_weight_sum"]
        ).to_numpy()
        if np.any(~np.isfinite(production)) or np.any(production <= 0):
            raise RuntimeError(
                f"Invalid profile near SZA={target_sza_deg:g} deg"
            )
        profiles.append(production)
        column_counts.append(len(selected_columns))
    return np.stack(profiles), np.asarray(column_counts)


def main() -> None:
    (
        altitude_km,
        sza_edges_deg,
        production,
        column_counts,
        data,
    ) = calculate_production_grid(INPUT_FILE)
    selected_profiles, selected_column_counts = (
        calculate_selected_sza_profiles(data, altitude_km)
    )
    log10_production = np.log10(production)

    configure_matplotlib()
    figure, (map_axis, profile_axis) = plt.subplots(
        1,
        2,
        figsize=(8.2, 4.0),
        constrained_layout=True,
        sharey=True,
        gridspec_kw={"width_ratios": (1.45, 0.90)},
    )
    image = map_axis.pcolormesh(
        sza_edges_deg,
        altitude_edges(altitude_km),
        log10_production,
        shading="flat",
        cmap="turbo",
        vmin=LOG10_Q_MIN,
        vmax=LOG10_Q_MAX,
        rasterized=True,
    )
    map_axis.axvline(
        90.0,
        color="white",
        linestyle="--",
        linewidth=1.0,
    )
    map_axis.text(
        91.5,
        altitude_km[-1] - 3.0,
        "Terminator",
        color="white",
        ha="left",
        va="top",
    )
    map_axis.set(
        xlabel="Solar zenith angle (deg)",
        ylabel="Altitude (km)",
        title="SZA and altitude distribution",
        xlim=(0.0, 180.0),
        ylim=(altitude_km[0], altitude_km[-1]),
        xticks=np.arange(0.0, 181.0, 30.0),
    )
    map_axis.text(
        0.02,
        0.98,
        "a",
        transform=map_axis.transAxes,
        ha="left",
        va="top",
        color="white",
        fontweight="bold",
        fontsize=9,
    )
    colorbar = figure.colorbar(image, ax=map_axis, pad=0.025)
    colorbar.set_label(
        r"$\log_{10}[Q_{\mathrm{hot\,O}}"
        r"\;(\mathrm{cm^{-3}\,s^{-1}})]$"
    )

    profile_colors = mpl.colormaps["viridis"](
        np.linspace(0.08, 0.92, TARGET_SZA_DEG.size)
    )
    for target_sza_deg, profile, color in zip(
        TARGET_SZA_DEG,
        selected_profiles,
        profile_colors,
    ):
        profile_axis.plot(
            profile,
            altitude_km,
            color=color,
            linewidth=1.5,
            label=f"{target_sza_deg:g}°",
        )
    profile_axis.set(
        xscale="log",
        xlabel=r"$Q_{\mathrm{hot\,O}}$ (cm$^{-3}$ s$^{-1}$)",
        title="Selected SZA profiles",
        xlim=(10.0**LOG10_Q_MIN, 10.0**LOG10_Q_MAX),
    )
    profile_axis.grid(
        True,
        which="both",
        color="0.88",
        linewidth=0.6,
    )
    profile_axis.legend(
        title="SZA",
        frameon=False,
        loc="upper right",
        ncol=1,
        handlelength=1.8,
    )
    profile_axis.text(
        0.02,
        0.98,
        "b",
        transform=profile_axis.transAxes,
        ha="left",
        va="top",
        color="0.15",
        fontweight="bold",
        fontsize=9,
    )
    figure.suptitle(
        "MGITM hot O photochemical production, "
        r"$L_s=0^\circ$, F070"
        "\nArea-weighted means",
        fontsize=10,
    )

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)

    print(f"output={OUTPUT.resolve()}")
    print(
        f"sza_bin_column_count_min={int(np.min(column_counts))}, "
        f"median={float(np.median(column_counts)):.1f}, "
        f"max={int(np.max(column_counts))}"
    )
    print(
        "selected_sza_column_counts="
        + ", ".join(
            f"{target:g}:{count}"
            for target, count in zip(
                TARGET_SZA_DEG, selected_column_counts
            )
        )
    )
    print(
        "production_range_cm3_s1="
        f"{np.min(production):.8e}, {np.max(production):.8e}"
    )


if __name__ == "__main__":
    main()
