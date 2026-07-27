"""Plot the baseline MGITM O2+, temperature, and hot O source profiles."""

from __future__ import annotations

from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plot_mgitm_hot_o_profiles import (  # noqa: E402
    INPUT_FILE,
    calculate_derived_profiles,
    load_nearest_subsolar_profile,
    make_figure,
)


OUTPUT = ROOT / "examples" / "figures" / "mgitm_ls000_f070_profiles.png"


def main() -> None:
    profile, metadata = load_nearest_subsolar_profile(INPUT_FILE)
    profile = calculate_derived_profiles(profile)
    figure = make_figure(profile, metadata)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=400, bbox_inches="tight")
    print(f"output={OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
