#!/usr/bin/env python3
# runs the full Hyades analysis from start to finish
# downloads star data from the Gaia archive, finds cluster members, saves a catalogue, makes three plots
#
# normal run (uses the cached download if it exists):
#   python run_analysis.py
#
# force a fresh download from the Gaia website:
#   python run_analysis.py --fresh

import argparse
import sys
from pathlib import Path

import pandas as pd

# make sure python can find the src folder regardless of where you run this from
sys.path.insert(0, str(Path(__file__).parent))

from src.query      import download_raw, RAW_OUTPUT
from src.membership import apply_quality_cuts, select_members
from src.plots      import plot_cmd, plot_pm_diagram, plot_parallax_histogram

PROCESSED_DIR = Path("data/processed")
MEMBERS_CSV   = PROCESSED_DIR / "members.csv"


def main(fresh: bool = False) -> None:
    print("=" * 60)
    print("  Gaia DR2 - Hyades Open Cluster Analysis")
    print("=" * 60)

    # step 1: grab the data (or load it from disk if we already downloaded it)
    raw_df = download_raw(output_path=RAW_OUTPUT, overwrite=fresh)
    print(f"[main] Raw catalogue: {len(raw_df):,} rows, {len(raw_df.columns)} columns\n")

    # step 2: apply quality cuts to get a clean field sample
    # we need this separately so we can use it as the background in the proper motion plot
    quality_df = apply_quality_cuts(raw_df)

    # step 3: run the full membership selection (quality cuts + box + sigma clipping)
    members_df = select_members(raw_df)

    if len(members_df) == 0:
        print("[main] ERROR: no members found - check membership thresholds.")
        sys.exit(1)

    # step 4: save the member catalogue to csv
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    members_df.to_csv(MEMBERS_CSV, index=False)
    print(f"\n[main] Member catalogue saved to {MEMBERS_CSV} ({len(members_df):,} rows)")

    # print a quick summary so we can sanity check the numbers look right
    print("\n-- Membership summary ------------------------------------------")
    for col, unit in [("parallax", "mas"), ("pmra", "mas/yr"), ("pmdec", "mas/yr")]:
        med = members_df[col].median()
        std = members_df[col].std()
        print(f"  {col:20s}: median={med:8.3f} {unit},  sd={std:.3f}")
    print(f"  {'distance (pc)':20s}: median={1000/members_df['parallax'].median():.1f} pc")
    print()

    # step 5: make the three figures
    # field stars for the PM plot = quality-cut stars that didn't make it into the member list
    field_df = quality_df[~quality_df["source_id"].isin(members_df["source_id"])].copy()

    print("[main] Generating figures ...")
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
