# makes the three figures for the Hyades analysis
# fig 1: colour-magnitude diagram - basically an HR diagram using gaia colours
# fig 2: proper motion plot - shows the cluster moving as a clump against the background stars
# fig 3: parallax histogram - shows the cluster poking up as a spike in the distance distribution

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd

FIGURES_DIR = Path("figures")
FIGURES_DIR.mkdir(exist_ok=True)

# colours used across all three plots
FIELD_COLOUR   = "#aaaaaa"  # grey for background/field stars
MEMBER_COLOUR  = "#2166ac"  # blue for cluster members
MEMBER_CMAP    = "plasma"   # colour scale for CMD points (coloured by parallax)

# clean up the default matplotlib style a bit
plt.rcParams.update({
    "font.family": "sans-serif",
    "axes.spines.top": False,    # remove the top border line
    "axes.spines.right": False,  # remove the right border line
    "figure.dpi": 150,
})


def plot_cmd(members: pd.DataFrame, output: Path = FIGURES_DIR / "fig1_cmd.png") -> None:
    # colour-magnitude diagram
    # x axis is colour (blue-red), y axis is how intrinsically bright the star is
    # the cluster should form a clear main sequence diagonal line
    # we colour each point by its parallax to show the cluster has some physical depth
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

    ax.invert_yaxis()  # brighter stars go at the top by convention
    ax.set_xlabel("BP - RP  (mag)", fontsize=12)
    ax.set_ylabel(r"$M_G$  (mag)", fontsize=12)
    ax.set_title(f"Hyades - Colour-Magnitude Diagram\n({len(members):,} members)", fontsize=12)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] Saved {output}")


def plot_pm_diagram(
    field: pd.DataFrame,
    members: pd.DataFrame,
    output: Path = FIGURES_DIR / "fig2_proper_motion.png",
) -> None:
    # proper motion diagram - each star is a dot at its (pmra, pmdec) coordinates
    # the Hyades all move together so they show up as a tight clump
    # field stars are spread all over the place
    # if the field is huge we only plot 5000 of them so the plot doesn't get too crowded
    fig, ax = plt.subplots(figsize=(7, 6))

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

    # crosshairs showing where the cluster centre sits
    ax.axvline(members["pmra"].median(), color=MEMBER_COLOUR, lw=0.8, ls="--", alpha=0.6)
    ax.axhline(members["pmdec"].median(), color=MEMBER_COLOUR, lw=0.8, ls="--", alpha=0.6)

    ax.set_xlabel(r"$\mu_{\alpha*}$  (mas yr$^{-1}$)", fontsize=12)
    ax.set_ylabel(r"$\mu_\delta$  (mas yr$^{-1}$)", fontsize=12)
    ax.set_title("Hyades - Proper-Motion Diagram", fontsize=12)
    ax.legend(fontsize=9, markerscale=2)

    # zoom in enough to see the cluster clump but also see how spread out the field is
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
    # parallax histogram - shows all the quality-cut stars in grey
    # and the cluster members in blue on top
    # the Hyades should show up as a clear spike around 21.5 mas
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
        density=True,  # normalise so both histograms are on the same scale
    )
    ax.hist(
        members["parallax"],
        bins=bins,
        color=MEMBER_COLOUR,
        alpha=0.85,
        label=f"Hyades members (n={len(members):,})",
        density=True,
    )

    # vertical line marking the median parallax of the cluster
    med_plx = members["parallax"].median()
    ax.axvline(med_plx, color=MEMBER_COLOUR, lw=1.5, ls="--",
               label=f"Member median: {med_plx:.1f} mas")

    ax.set_xlabel("Parallax  (mas)", fontsize=12)
    ax.set_ylabel("Normalised count  (density)", fontsize=12)
    ax.set_title("Hyades - Parallax Distribution", fontsize=12)
    ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"[plots] Saved {output}")
