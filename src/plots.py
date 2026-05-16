"""
plots.py – Produce the three diagnostic figures for the Hyades analysis.

Figure 1 : Colour-magnitude diagram  (BP-RP vs absolute G)
Figure 2 : Proper-motion vector diagram  (field grey, members coloured)
Figure 3 : Parallax histogram  (full quality-cut sample + member peak)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

# Colour scheme
FIELD_COLOUR   = "#aaaaaa"
MEMBER_COLOUR  = "#2166ac"   # blue
MEMBER_CMAP    = "plasma"     # for CMD points coloured by parallax

plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def plot_cmd(members: pd.DataFrame, output: Path = FIGURES_DIR / "fig1_cmd.png") -> None:
    """
    Figure 1: Colour-magnitude diagram.

    x-axis : BP-RP colour index (mag)
    y-axis : Absolute G magnitude M_G (inverted so brighter = up)
    Points are coloured by parallax to show the cluster depth.
    """
    fig, ax = plt.subplots(figsize=(6, 7))

    sc = ax.scatter(
        members["bp_rp"],
        members["abs_g_mag"],
        c=members["parallax"],
        cmap=MEMBER_CMAP,
        s=12,
        alpha=0.85,
        linewidths=0,
        vmin=members["parallax"].quantile(0.05),
        vmax=members["parallax"].quantile(0.95),
    )

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Parallax (mas)", fontsize=10)

    ax.invert_yaxis()  # brighter stars at top (convention)
    ax.set_xlabel("BP – RP  (mag)", fontsize=12)
    ax.set_ylabel(r"$M_G$  (mag)", fontsize=12)
    ax.set_title(f"Hyades — Colour-Magnitude Diagram\n({len(members):,} members)", fontsize=12)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] Saved {output}")


def plot_pm_diagram(
    field: pd.DataFrame,
    members: pd.DataFrame,
    output: Path = FIGURES_DIR / "fig2_proper_motion.png",
) -> None:
    """
    Figure 2: Proper-motion vector diagram.

    Background field sources (quality-cut but not members) are shown in grey.
    Hyades members are overplotted in blue with a filled circle.
    """
    fig, ax = plt.subplots(figsize=(7, 6))

    # Field stars (down-sampled to at most 5000 for clarity)
    n_field = min(5000, len(field))
    field_sample = field.sample(n=n_field, random_state=42) if len(field) > n_field else field

    ax.scatter(
        field_sample["pmra"],
        field_sample["pmdec"],
        s=3,
        c=FIELD_COLOUR,
        alpha=0.4,
        linewidths=0,
        label=f"Field (n={len(field):,}, showing {n_field:,})",
        rasterized=True,
    )

    ax.scatter(
        members["pmra"],
        members["pmdec"],
        s=15,
        c=MEMBER_COLOUR,
        alpha=0.85,
        linewidths=0,
        label=f"Hyades members (n={len(members):,})",
        zorder=3,
    )

    # Annotate cluster centroid
    ax.axvline(members["pmra"].median(), color=MEMBER_COLOUR, lw=0.8, ls="--", alpha=0.6)
    ax.axhline(members["pmdec"].median(), color=MEMBER_COLOUR, lw=0.8, ls="--", alpha=0.6)

    ax.set_xlabel(r"$\mu_{\alpha*}$  (mas yr$^{-1}$)", fontsize=12)
    ax.set_ylabel(r"$\mu_\delta$  (mas yr$^{-1}$)", fontsize=12)
    ax.set_title("Hyades — Proper-Motion Diagram", fontsize=12)
    ax.legend(fontsize=9, markerscale=2)

    # Zoom to a window that shows both the field spread and the Hyades clump
    pmra_med  = members["pmra"].median()
    pmdec_med = members["pmdec"].median()
    ax.set_xlim(pmra_med  - 30, pmra_med  + 30)
    ax.set_ylim(pmdec_med - 25, pmdec_med + 25)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] Saved {output}")


def plot_parallax_histogram(
    quality_sample: pd.DataFrame,
    members: pd.DataFrame,
    output: Path = FIGURES_DIR / "fig3_parallax_histogram.png",
) -> None:
    """
    Figure 3: Parallax histogram.

    The full quality-cut sample (grey) is overlaid with the member sub-sample
    (blue) to show the Hyades peak standing out from the field distribution.
    """
    fig, ax = plt.subplots(figsize=(7, 5))

    bins = np.linspace(
        quality_sample["parallax"].min(),
        quality_sample["parallax"].max(),
        80,
    )

    ax.hist(
        quality_sample["parallax"],
        bins=bins,
        color=FIELD_COLOUR,
        alpha=0.7,
        label=f"Quality-cut sample (n={len(quality_sample):,})",
        density=True,
    )
    ax.hist(
        members["parallax"],
        bins=bins,
        color=MEMBER_COLOUR,
        alpha=0.85,
        label=f"Hyades members (n={len(members):,})",
        density=True,
    )

    # Mark the cluster median
    med_plx = members["parallax"].median()
    ax.axvline(med_plx, color=MEMBER_COLOUR, lw=1.5, ls="--",
               label=f"Member median: {med_plx:.1f} mas")

    ax.set_xlabel("Parallax  (mas)", fontsize=12)
    ax.set_ylabel("Normalised count  (density)", fontsize=12)
    ax.set_title("Hyades — Parallax Distribution", fontsize=12)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] Saved {output}")
