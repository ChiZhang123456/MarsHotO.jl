"""Plot MGITM ion and electron temperature profiles near the subsolar point."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_mgitm_hot_o_profiles import (  # noqa: E402
    INPUT_FILE,
    load_nearest_subsolar_profile,
)


OUTPUT = (
    ROOT
    / "examples"
    / "figures"
    / "mgitm_ls000_f070_ion_electron_temperature_profiles.png"
)


def configure_matplotlib() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "dejavusans",
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.linewidth": 1.0,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "legend.frameon": False,
        }
    )


def make_figure(profile, metadata: dict) -> plt.Figure:
    altitude_km = profile["altitude_km"].to_numpy(dtype=float)
    ion_temperature_k = profile["Ti_K"].to_numpy(dtype=float)
    electron_temperature_k = profile["Te_K"].to_numpy(dtype=float)

    for name, values in {
        "altitude": altitude_km,
        "ion temperature": ion_temperature_k,
        "electron temperature": electron_temperature_k,
    }.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"MGITM {name} profile contains nonfinite values")

    configure_matplotlib()
    fig, axis = plt.subplots(figsize=(5.2, 5.4), constrained_layout=True)

    axis.plot(
        ion_temperature_k,
        altitude_km,
        color="#D17C2F",
        linewidth=2.2,
        label=r"Ion temperature, $T_i$",
    )
    axis.plot(
        electron_temperature_k,
        altitude_km,
        color="#B33B45",
        linewidth=2.2,
        label=r"Electron temperature, $T_e$",
    )
    axis.set_xlabel("Temperature (K)")
    axis.set_ylabel("Altitude (km)")
    axis.set_xlim(left=0.0)
    axis.grid(True, color="0.88", linewidth=0.7)
    axis.legend(loc="upper left")

    selected_lats = ", ".join(
        f"{value:g}" for value in metadata["selected_latitudes_deg"]
    )
    selected_lons = ", ".join(
        f"{value:g}" for value in metadata["selected_longitudes_deg"]
    )
    axis.set_title(
        "MGITM plasma temperatures, "
        r"$L_s=0^\circ$, F070"
        f"\nlon = {selected_lons}°, lat = {selected_lats}°, "
        f"SZA = {metadata['minimum_grid_sza_deg']:.2f}°"
    )
    return fig


def main() -> None:
    profile, metadata = load_nearest_subsolar_profile(INPUT_FILE)
    figure = make_figure(profile, metadata)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    plt.close(figure)

    print(f"input={INPUT_FILE.resolve()}")
    print(f"output={OUTPUT.resolve()}")
    print(
        "selected_columns="
        f"lon={metadata['selected_longitudes_deg']}, "
        f"lat={metadata['selected_latitudes_deg']}, "
        f"SZA={metadata['minimum_grid_sza_deg']:.4f} deg"
    )


if __name__ == "__main__":
    main()
