"""Schematic derivation of dR = F N (d sigma / d Omega) d Omega."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle, Wedge
import numpy as np


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 10,
        "axes.spines.left": False,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.spines.bottom": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
    }
)

BLUE = "#2878A5"
ORANGE = "#D96B27"
GREEN = "#4B9560"
GRAY = "#5B6573"


def arrow(ax, xy1, xy2, color=BLUE, lw=1.4, alpha=1.0):
    ax.add_patch(
        FancyArrowPatch(
            xy1,
            xy2,
            arrowstyle="-|>",
            mutation_scale=10,
            lw=lw,
            color=color,
            alpha=alpha,
        )
    )


def clean(ax, xlim=(0, 10), ylim=(0, 7)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def main():
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.8), constrained_layout=True)

    # 1. Flux times effective area.
    ax = axes[0]
    clean(ax)
    ax.set_title("1. One target", fontweight="bold")
    for y in np.linspace(1.0, 6.0, 9):
        arrow(ax, (0.3, y), (8.8, y), alpha=0.45)
    target = Circle((6.0, 3.5), 0.22, fc=ORANGE, ec="none", zorder=4)
    effective = Circle((6.0, 3.5), 1.25, fc=ORANGE, ec=ORANGE, alpha=0.22, lw=1.5)
    ax.add_patch(effective)
    ax.add_patch(target)
    ax.text(1.0, 6.35, r"flux $F$", color=BLUE, fontsize=12)
    ax.annotate(
        r"effective area $\sigma$",
        xy=(6.65, 4.45),
        xytext=(3.5, 6.25),
        arrowprops=dict(arrowstyle="->", color=ORANGE),
        color=ORANGE,
        ha="center",
    )
    ax.text(5.0, 0.2, r"collision rate $=F\sigma$", ha="center", fontsize=12)

    # 2. N independent targets.
    ax = axes[1]
    clean(ax)
    ax.set_title(r"2. $N$ targets", fontweight="bold")
    ax.add_patch(Rectangle((1.0, 0.8), 7.8, 5.4, fc="#EEF3F6", ec="0.75"))
    pts = [(2.0, 1.6), (4.0, 1.8), (6.4, 1.4), (7.8, 2.7), (2.6, 3.5), (5.2, 3.2), (7.0, 4.6), (3.7, 5.2)]
    for x, y in pts:
        ax.add_patch(Circle((x, y), 0.18, fc=ORANGE, ec="none"))
    for y in np.linspace(1.2, 5.8, 7):
        arrow(ax, (0.15, y), (9.4, y), alpha=0.35)
    ax.text(5.0, 0.15, r"total collision rate $=FN\sigma$", ha="center", fontsize=12)
    ax.text(5.0, 6.45, "thin target: targets act independently", ha="center", color=GRAY, fontsize=9)

    # 3. Only a small outgoing solid angle.
    ax = axes[2]
    clean(ax)
    ax.set_title(r"3. Select directions in $d\Omega$", fontweight="bold")
    origin = (3.0, 3.1)
    arrow(ax, (0.3, 3.1), origin, lw=2.0)
    ax.add_patch(Circle(origin, 0.2, fc=ORANGE, ec="none", zorder=5))
    for deg in [-65, -35, -10, 20, 48, 72]:
        th = np.deg2rad(deg)
        arrow(ax, origin, (origin[0] + 3.4 * np.cos(th), origin[1] + 3.4 * np.sin(th)), color=GRAY, alpha=0.55)
    ax.add_patch(Wedge(origin, 3.7, 30, 48, fc=GREEN, ec=GREEN, alpha=0.25, lw=1.8))
    for deg in [33, 39, 45]:
        th = np.deg2rad(deg)
        arrow(ax, origin, (origin[0] + 3.5 * np.cos(th), origin[1] + 3.5 * np.sin(th)), color=GREEN, lw=2)
    ax.text(6.1, 5.85, r"small direction patch $d\Omega$", color=GREEN, ha="center")
    ax.text(5.1, 0.25, r"cross section in patch: $d\sigma$", ha="center", fontsize=12)

    # 4. Definition and dimensional chain.
    ax = axes[3]
    clean(ax)
    ax.set_title("4. Count only that patch", fontweight="bold")
    ax.text(5, 5.75, r"$\dfrac{d\sigma}{d\Omega}$", ha="center", va="center", fontsize=24, color=GREEN)
    ax.text(5, 4.35, "cross section per unit solid angle", ha="center", color=GRAY)
    ax.text(5, 3.35, r"$d\sigma=\dfrac{d\sigma}{d\Omega}\,d\Omega$", ha="center", fontsize=17)
    ax.text(5, 2.25, r"$dR=FN\,d\sigma$", ha="center", fontsize=17)
    ax.text(
        5,
        1.05,
        r"$dR=FN\dfrac{d\sigma}{d\Omega}d\Omega$",
        ha="center",
        fontsize=17,
        color=ORANGE,
        bbox=dict(boxstyle="round,pad=0.35", fc="#FFF5ED", ec=ORANGE),
    )
    ax.text(5, 0.05, r"units: s$^{-1}$", ha="center", color=GRAY)

    fig.suptitle("From particle flux to an angular scattering rate", fontsize=14, fontweight="bold")
    output = Path(__file__).resolve().parent / "figures"
    output.mkdir(parents=True, exist_ok=True)
    stem = output / "cross_section_event_rate_explained"
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
