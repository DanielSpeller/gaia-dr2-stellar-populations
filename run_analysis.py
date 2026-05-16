#!/usr/bin/env python3
"""
run_analysis.py – End-to-end pipeline for Hyades membership from Gaia DR2.

Usage
-----
    python run_analysis.py           # normal run (uses cached raw data if present)
    python run_analysis.py --fresh   # force re-download from ESA TAP

Outputs
-------
    data/raw/hyades_gaia_dr2_raw.csv      (gitignored – large, reproducible)
    data/processed/members.csv            (committed – cleaned member catalogue)
    figures/fig1_cmd.png
    figures/fig2_proper_motion.png
    figures/fig3_parallax_histogram.png
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Add project root to path so `src` is importable when run from any directory
sys.path.insert(0, str(Path(__file__).parent))

from src.query      import download_raw, RAW_OUTPUT
from src.membership import apply_quality_cuts, select_members
from src.plots      import plot_cmd, plot_pm_diagram, plot_parallax_histogram

PROCESSED_DIR = Path("data/processed")
MEMBERS_CSV   = PROCESSED_DIR / "members.csv"


def main(fresh: bool = False) -> None:
    print("=" * 60)
    print("  Gaia DR2 – Hyades Open Cluster Analysis")
    print("=" * 60)

    # ── Step 1: Download (or load cache) ─────────────────────────────────────
    raw_df = download_raw(output_path=RAW_OUTPUT, overwrite=fresh)
    print(f"[main] Raw catalogue: {len(raw_df):,} rows, {len(raw_df.columns)} columns\n")

    # ── Step 2: Quality-cut sample (used as field background in PM plot) ──────
    quality_df = apply_quality_cuts(raw_df)

    # ── Step 3: Full membership selection ────────────────────────────────────
    members_df = select_members(raw_df)

    if len(members_df) == 0:
        print("[main] ERROR: no members found – check membership thresholds.")
        sys.exit(1)

    # ── Step 4: Save processed catalogue ─────────────────────────────────────
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    members_df.to_csv(MEMBERS_CSV, index=False)
    print(f"\n[main] Member catalogue saved to {MEMBERS_CSV} ({len(members_df):,} rows)")

    # Print a brief summary
    print("\n── Membership summary ──────────────────────────────────────────")
    for col, unit in [("parallax", "mas"), ("pmra", "mas/yr"), ("pmdec", "mas/yr")]:
        med = members_df[col].median()
        std = members_df[col].std()
        print(f"  {col:20s}: median={med:8.3f} {unit},  σ={std:.3f}")
    print(f"  {'distance (pc)':20s}: median={1000/members_df['parallax'].median():.1f} pc")
    print()

    # ── Step 5: Produce figures ───────────────────────────────────────────────
    # Field for PM diagram = quality-cut sources that are NOT members
    field_df = quality_df[~quality_df["source_id"].isin(members_df["source_id"])].copy()

    print("[main] Generating figures …")
    plot_cmd(members_df)
    plot_pm_diagram(field_df, members_df)
    plot_parallax_histogram(quality_df, members_df)

    print("\n[main] Done.  All outputs written successfully.")
    print("  Figures : figures/fig1_cmd.png, fig2_proper_motion.png, fig3_parallax_histogram.png")
    print(f"  Members : {MEMBERS_CSV}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Hyades Gaia DR2 analysis.")
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Force re-download from ESA TAP even if local cache exists.",
    )
    args = parser.parse_args()
    main(fresh=args.fresh)
