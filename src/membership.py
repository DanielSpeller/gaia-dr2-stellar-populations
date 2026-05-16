"""
membership.py – Quality cuts and cluster membership selection for the Hyades.

Strategy
--------
1. Hard quality cuts applied to the raw DR2 catalogue (astrometric & photometric).
2. A coarse PM + parallax box isolates the Hyades locus from the field.
3. Iterative 3-sigma clipping on (pmra, pmdec, parallax) refines membership.

Hyades reference values (van Leeuwen 2009; Gaia Collaboration 2018b)
--------------------------------------------------------------------
  parallax : 21.5 ± 1.5 mas   (spread reflects ~10 pc physical depth)
  pmra     : +104.0 ± 2.0 mas/yr
  pmdec    :  -28.0 ± 2.0 mas/yr
"""

import numpy as np
import pandas as pd

# ── Hyades expected kinematics (used as starting point for sigma-clipping) ──
HYADES_PARALLAX   = 21.5   # mas
HYADES_PMRA       = 104.0  # mas/yr
HYADES_PMDEC      = -28.0  # mas/yr

# ── Quality-cut thresholds ───────────────────────────────────────────────────
MIN_PARALLAX_SNR      = 5.0   # parallax_over_error > 5 → relative error < 20 %
MAX_ASTR_EXCESS_NOISE = 1.0   # mas; larger values indicate bad astrometric fit
MAX_EXCESS_NOISE_SIG  = 2.0   # significance of excess noise
MIN_VISIBILITY_PERIODS = 7    # independent observing epochs

MIN_G_FLUX_SNR        = 50.0  # phot_g_mean_flux_over_error > 50
MIN_BP_FLUX_SNR       = 10.0  # BP & RP SNR can be lower for faint/red stars
MIN_RP_FLUX_SNR       = 10.0
MAX_BP_RP_EXCESS      = 1.6   # phot_bp_rp_excess_factor; filters blends

# ── Membership box before sigma-clipping ────────────────────────────────────
# Wide enough to capture members with measurement errors but reject the field.
PM_BOX_PMRA_MIN   =  85.0   # mas/yr
PM_BOX_PMRA_MAX   = 125.0
PM_BOX_PMDEC_MIN  = -45.0
PM_BOX_PMDEC_MAX  = -10.0

# Parallax already pre-filtered to 12–35 mas by ADQL; tighten here.
PARALLAX_BOX_MIN  = 16.0    # mas
PARALLAX_BOX_MAX  = 27.0    # mas

# ── Sigma-clipping parameters ────────────────────────────────────────────────
SIGMA_CLIP_NSIGMA  = 3.0
SIGMA_CLIP_MAXITER = 10


def apply_quality_cuts(df: pd.DataFrame) -> pd.DataFrame:
    """Return rows that pass astrometric and photometric quality thresholds."""
    n0 = len(df)

    # Astrometric quality
    mask = (
        (df["parallax_over_error"] > MIN_PARALLAX_SNR)          # good parallax SNR
        & (df["astrometric_excess_noise"] < MAX_ASTR_EXCESS_NOISE)  # clean solution
        & (df["astrometric_excess_noise_sig"] < MAX_EXCESS_NOISE_SIG)
        & (df["visibility_periods_used"] >= MIN_VISIBILITY_PERIODS)  # enough epochs
    )

    # Photometric quality (flux SNR + BP/RP excess factor)
    mask &= (
        (df["phot_g_mean_flux_over_error"] > MIN_G_FLUX_SNR)
        & (df["phot_bp_mean_flux_over_error"] > MIN_BP_FLUX_SNR)
        & (df["phot_rp_mean_flux_over_error"] > MIN_RP_FLUX_SNR)
        & (df["phot_bp_rp_excess_factor"] < MAX_BP_RP_EXCESS)  # reject blends
        & df["phot_bp_mean_mag"].notna()
        & df["phot_rp_mean_mag"].notna()
    )

    result = df[mask].copy()
    print(f"[membership] Quality cuts: {n0:,} → {len(result):,} sources retained.")
    return result


def apply_pm_parallax_box(df: pd.DataFrame) -> pd.DataFrame:
    """Coarse box cut in proper-motion and parallax space."""
    n0 = len(df)
    mask = (
        (df["pmra"]     >= PM_BOX_PMRA_MIN)  & (df["pmra"]     <= PM_BOX_PMRA_MAX)
        & (df["pmdec"]  >= PM_BOX_PMDEC_MIN) & (df["pmdec"]    <= PM_BOX_PMDEC_MAX)
        & (df["parallax"] >= PARALLAX_BOX_MIN) & (df["parallax"] <= PARALLAX_BOX_MAX)
    )
    result = df[mask].copy()
    print(f"[membership] PM+parallax box:  {n0:,} → {len(result):,} candidates.")
    return result


def sigma_clip_members(df: pd.DataFrame) -> pd.DataFrame:
    """
    Iterative 3-sigma clipping in (pmra, pmdec, parallax).

    In each iteration the median and MAD-based sigma are recomputed from
    surviving members only, converging on the tight kinematic Hyades locus.
    """
    mask = np.ones(len(df), dtype=bool)
    cols = ["pmra", "pmdec", "parallax"]

    for iteration in range(SIGMA_CLIP_MAXITER):
        n_before = mask.sum()
        sub = df[mask]
        for col in cols:
            med = sub[col].median()
            # Robust sigma via median absolute deviation (MAD × 1.4826)
            mad = (sub[col] - med).abs().median()
            sigma = mad * 1.4826 if mad > 0 else sub[col].std()
            lo = med - SIGMA_CLIP_NSIGMA * sigma
            hi = med + SIGMA_CLIP_NSIGMA * sigma
            # Apply clip only to currently surviving rows
            idx = df.index[mask]
            drop = df.loc[idx, col].lt(lo) | df.loc[idx, col].gt(hi)
            mask[df.index.get_indexer(idx[drop])] = False

        n_after = mask.sum()
        print(f"[membership] Sigma-clip iter {iteration + 1}: {n_before} → {n_after}")
        if n_after == n_before:
            break  # converged

    result = df[mask].copy()
    print(f"[membership] Final member count: {len(result):,}")
    return result


def compute_absolute_magnitude(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add absolute G magnitude column using the measured parallax.

    M_G = G + 5 + 5 * log10(parallax_mas / 1000)
        = G + 5 * log10(parallax_mas) - 10
    """
    df = df.copy()
    df["abs_g_mag"] = (
        df["phot_g_mean_mag"]
        + 5.0 * np.log10(df["parallax"] / 1000.0)
        + 5.0
    )
    return df


def select_members(raw_df: pd.DataFrame) -> pd.DataFrame:
    """Full membership pipeline: quality cuts → box → sigma-clip → M_G."""
    df = apply_quality_cuts(raw_df)
    df = apply_pm_parallax_box(df)
    df = sigma_clip_members(df)
    df = compute_absolute_magnitude(df)
    return df
