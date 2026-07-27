"""Plot baseline MGITM ion, temperature, and hot O production profiles."""

from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
INPUT_FILE = ROOT / "MGITM" / "MGITM_LS000_F070_150901.dat"
OUTPUT_DIR = ROOT / "examples" / "figures"
OUTPUT_STEM = OUTPUT_DIR / "mgitm_ls000_f070_subsolar_hot_o_profiles"
SOURCE_DATA_FILE = OUTPUT_DIR / "mgitm_ls000_f070_subsolar_hot_o_profiles.csv"

MGITM_COLUMNS = [
    "longitude_deg",
    "latitude_deg",
    "altitude_km",
    "Tn_K",
    "Ti_K",
    "Te_K",
    "nCO2_m3",
    "nO_m3",
    "nN2_m3",
    "nCO_m3",
    "nO2_m3",
    "nO2p_m3",
    "nOp_m3",
    "nCO2p_m3",
    "ne_m3",
    "UN_ms",
    "VN_ms",
    "WN_ms",
]


def read_subsolar_longitude(path: Path) -> float:
    with path.open("r", encoding="utf-8") as stream:
        for line in stream:
            match = re.search(r"subsolar_longitude:\s*([+-]?\d+(?:\.\d+)?)", line)
            if match:
                return float(match.group(1))
    raise ValueError(f"Subsolar longitude was not found in {path}")


def angular_distance_deg(longitude_deg: np.ndarray, reference_deg: float) -> np.ndarray:
    return (longitude_deg - reference_deg + 180.0) % 360.0 - 180.0


def load_nearest_subsolar_profile(path: Path) -> tuple[pd.DataFrame, dict]:
    subsolar_lon_deg = read_subsolar_longitude(path)
    data = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=MGITM_COLUMNS,
        engine="c",
    )

    surface_columns = data[["longitude_deg", "latitude_deg"]].drop_duplicates()
    delta_lon_rad = np.deg2rad(
        angular_distance_deg(surface_columns["longitude_deg"].to_numpy(), subsolar_lon_deg)
    )
    latitude_rad = np.deg2rad(surface_columns["latitude_deg"].to_numpy())
    cos_sza = np.cos(latitude_rad) * np.cos(delta_lon_rad)
    sza_deg = np.rad2deg(np.arccos(np.clip(cos_sza, -1.0, 1.0)))
    minimum_sza_deg = float(np.min(sza_deg))

    tolerance = 1.0e-10
    nearest_columns = surface_columns.loc[
        np.isclose(sza_deg, minimum_sza_deg, rtol=0.0, atol=tolerance)
    ].copy()
    selected = data.merge(
        nearest_columns,
        on=["longitude_deg", "latitude_deg"],
        how="inner",
    )

    profile = (
        selected.groupby("altitude_km", as_index=False)
        .mean(numeric_only=True)
        .sort_values("altitude_km")
        .reset_index(drop=True)
    )

    metadata = {
        "subsolar_longitude_deg": subsolar_lon_deg,
        "minimum_grid_sza_deg": minimum_sza_deg,
        "number_of_averaged_columns": len(nearest_columns),
        "selected_longitudes_deg": sorted(
            nearest_columns["longitude_deg"].unique().tolist()
        ),
        "selected_latitudes_deg": sorted(
            nearest_columns["latitude_deg"].unique().tolist()
        ),
    }
    return profile, metadata


def dissociative_recombination_coefficient_cm3s(te_k: np.ndarray) -> np.ndarray:
    te_k = np.asarray(te_k, dtype=float)
    return np.where(
        te_k <= 1200.0,
        1.95e-7 * (300.0 / te_k) ** 0.70,
        7.39e-8 * (1200.0 / te_k) ** 0.56,
    )


def calculate_derived_profiles(profile: pd.DataFrame) -> pd.DataFrame:
    result = profile.copy()
    result["nO2p_cm3"] = result["nO2p_m3"] * 1.0e-6
    result["ne_cm3"] = result["ne_m3"] * 1.0e-6
    result["k_dr_cm3s"] = dissociative_recombination_coefficient_cm3s(
        result["Te_K"].to_numpy()
    )
    result["alpha_dr_cm3s"] = (
        result["ne_cm3"] * result["nO2p_cm3"] * result["k_dr_cm3s"]
    )
    result["Q_hotO_cm3s"] = 2.0 * result["alpha_dr_cm3s"]
    return result


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 7.5,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
        }
    )


def make_figure(profile: pd.DataFrame, metadata: dict) -> plt.Figure:
    configure_matplotlib()
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(7.2, 4.0),
        sharey=True,
        constrained_layout=True,
    )

    altitude = profile["altitude_km"]

    axes[0].plot(profile["nO2p_cm3"], altitude, color="#2F6B9A", lw=2.0)
    axes[0].set_xscale("log")
    axes[0].set_xlabel(r"$n_{\mathrm{O_2^+}}$ (cm$^{-3}$)")
    axes[0].set_ylabel("Altitude (km)")
    axes[0].grid(True, which="both", color="0.90", lw=0.6)

    temperature_colors = {
        "Tn_K": "#3D3D3D",
        "Ti_K": "#D17C2F",
        "Te_K": "#B33B45",
    }
    temperature_labels = {
        "Tn_K": r"$T_n$",
        "Ti_K": r"$T_i$",
        "Te_K": r"$T_e$",
    }
    for column in ["Tn_K", "Ti_K", "Te_K"]:
        axes[1].plot(
            profile[column],
            altitude,
            lw=1.8,
            color=temperature_colors[column],
            label=temperature_labels[column],
        )
    axes[1].set_xlabel("Temperature (K)")
    axes[1].legend(loc="best")
    axes[1].grid(True, color="0.90", lw=0.6)

    axes[2].plot(profile["Q_hotO_cm3s"], altitude, color="#6A4C93", lw=2.0)
    axes[2].set_xscale("log")
    axes[2].set_xlabel(r"$Q_{\mathrm{hot\,O}}$ (cm$^{-3}$ s$^{-1}$)")
    axes[2].grid(True, which="both", color="0.90", lw=0.6)

    for label, axis in zip(["a", "b", "c"], axes):
        axis.text(
            0.02,
            0.98,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontweight="bold",
            fontsize=9,
        )

    selected_lats = ", ".join(f"{value:g}" for value in metadata["selected_latitudes_deg"])
    selected_lons = ", ".join(f"{value:g}" for value in metadata["selected_longitudes_deg"])
    fig.suptitle(
        "MGITM baseline profiles, "
        r"$L_s=0^\circ$, F070"
        f"\nNearest subsolar grid columns: lon = {selected_lons}°, "
        f"lat = {selected_lats}°, SZA = {metadata['minimum_grid_sza_deg']:.2f}°",
        fontsize=9,
    )
    return fig


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    profile, metadata = load_nearest_subsolar_profile(INPUT_FILE)
    profile = calculate_derived_profiles(profile)
    profile.to_csv(SOURCE_DATA_FILE, index=False)

    fig = make_figure(profile, metadata)
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUTPUT_STEM.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    peak_index = profile["Q_hotO_cm3s"].idxmax()
    print(f"Input: {INPUT_FILE}")
    print(f"Output: {OUTPUT_STEM.with_suffix('.png')}")
    print(f"Source data: {SOURCE_DATA_FILE}")
    print(f"Nearest grid SZA: {metadata['minimum_grid_sza_deg']:.4f} deg")
    print(
        "Selected columns: "
        f"lon={metadata['selected_longitudes_deg']}, "
        f"lat={metadata['selected_latitudes_deg']}"
    )
    print(
        "Peak hot O production: "
        f"{profile.loc[peak_index, 'Q_hotO_cm3s']:.6e} cm^-3 s^-1 "
        f"at {profile.loc[peak_index, 'altitude_km']:.2f} km"
    )


if __name__ == "__main__":
    main()
