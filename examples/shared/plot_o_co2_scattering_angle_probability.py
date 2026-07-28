"""Plot laboratory-frame O + CO2 direct-scattering angular distributions.

Rice data source:
https://www.ruf.rice.edu/~atmol/direct_data.html

The angular probability density is obtained from the differential cross
section using

    p(theta) = [2*pi*sin(theta)*(d sigma/d Omega)] / sigma_range

and is reported per degree.  Each curve is normalized only over the angular
range tabulated for that projectile energy.
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path
from urllib.request import urlopen

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator


BASE_URL = (
    "https://www.ruf.rice.edu/~atmol/data/direct_scattering/O/"
    "O_CO2_{energy}.txt"
)
ENERGIES = (100, 500, 1500)
COLORS = {100: "#4477AA", 500: "#228833", 1500: "#CC6677"}


def read_rice_table(energy: int) -> tuple[np.ndarray, np.ndarray]:
    """Return laboratory angle (degree) and d sigma/d Omega."""
    with urlopen(BASE_URL.format(energy=energy), timeout=30) as response:
        text = response.read().decode("latin-1").replace("\r", "\n")
    lines = [line for line in text.splitlines() if line.strip()]
    data = np.genfromtxt(StringIO("\n".join(lines[1:])), delimiter=",")
    return data[:, 0], data[:, 2] if data.shape[1] >= 3 else data[:, 1]


def angular_pdf(
    angle_deg: np.ndarray, differential_cross_section: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate the cross section and return a normalized PDF per degree."""
    dense_angle = np.linspace(angle_deg.min(), angle_deg.max(), 2000)
    log_interp = PchipInterpolator(
        angle_deg, np.log(differential_cross_section)
    )
    dense_dcs = np.exp(log_interp(dense_angle))
    theta_rad = np.deg2rad(dense_angle)

    # d sigma / d theta in units of cross section per radian.
    dsigma_dtheta_rad = 2.0 * np.pi * np.sin(theta_rad) * dense_dcs
    sigma_range = np.trapz(dsigma_dtheta_rad, theta_rad)

    # Convert the normalized density from radian^-1 to degree^-1.
    pdf_per_degree = dsigma_dtheta_rad / sigma_range * np.pi / 180.0
    cdf = cumulative_trapezoid(
        pdf_per_degree, dense_angle, initial=0.0
    )
    cdf /= cdf[-1]
    return dense_angle, pdf_per_degree, cdf


def main() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial"],
            "font.size": 8,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    distributions = {}
    for energy in ENERGIES:
        angle, dcs = read_rice_table(energy)
        distributions[energy] = (angle, dcs, *angular_pdf(angle, dcs))

    fig, ax = plt.subplots(figsize=(4.4, 3.4))

    for energy in ENERGIES:
        _, _, dense_angle, _, cdf = distributions[energy]
        label = f"{energy} eV"
        random_number = np.linspace(0.0, 1.0, 2000)
        sampled_angle = np.interp(random_number, cdf, dense_angle)
        ax.plot(
            random_number,
            sampled_angle,
            linewidth=1.6,
            color=COLORS[energy],
            label=label,
        )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.2, 180.0)
    ax.set_yscale("log")
    ax.set_xlabel("Random number, $R$")
    ax.set_ylabel(r"Laboratory scattering angle, $\theta$ (deg)")
    ax.set_title(r"Inverse-CDF sampling: $\theta=F^{-1}(R)$", loc="left")
    ax.legend(title="O projectile energy")

    output_dir = Path(__file__).resolve().parent.parent / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / "o_co2_scattering_angle_probability"
    fig.savefig(stem.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    for energy in ENERGIES:
        _, _, dense_angle, pdf, cdf = distributions[energy]
        median = np.interp(0.5, cdf, dense_angle)
        mode = dense_angle[np.argmax(pdf)]
        print(
            f"{energy:4d} eV: range={dense_angle[0]:.3f}-"
            f"{dense_angle[-1]:.3f} deg, mode={mode:.3f} deg, "
            f"median={median:.3f} deg, integral="
            f"{np.trapz(pdf, dense_angle):.8f}"
        )


if __name__ == "__main__":
    main()
