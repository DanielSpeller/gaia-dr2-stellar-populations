# this file figures out which stars actually belong to the Hyades cluster
# and which ones just happen to be in the same patch of sky by coincidence
#
# the plan is:
# 1. throw away any stars with dodgy measurements (quality cuts)
# 2. keep only stars moving roughly the same direction as the Hyades (box cut)
# 3. tighten that up by repeatedly removing outliers until it stabilises (sigma clipping)

import numpy as np
import pandas as pd

# roughly where the Hyades cluster sits in parallax/proper-motion space
# these are just starting values - the sigma clipping will find the real centre
HYADES_PARALLAX   = 21.5   # mas - about 46 pc away
HYADES_PMRA       = 104.0  # mas/yr - moving fast to the right on the sky
HYADES_PMDEC      = -28.0  # mas/yr - and slightly downward

# how strict to be on measurement quality
# these thresholds come from the Gaia DR2 documentation recommendations
MIN_PARALLAX_SNR      = 5.0   # parallax needs to be 5x its own error - otherwise it's just noise
MAX_ASTR_EXCESS_NOISE = 1.0   # if the position fit left this much noise behind, something's off
MAX_EXCESS_NOISE_SIG  = 2.0   # how statistically significant that leftover noise is
MIN_VISIBILITY_PERIODS = 7    # needs at least 7 separate visits from Gaia

MIN_G_FLUX_SNR        = 50.0  # brightness needs to be well-measured in the main G band
MIN_BP_FLUX_SNR       = 10.0  # blue and red can be a bit noisier, especially for faint red stars
MIN_RP_FLUX_SNR       = 10.0
MAX_BP_RP_EXCESS      = 1.6   # if blue+red brightness adds up to way more than G, it's probably two stars on top of each other

# the rough box in proper motion / parallax space where we expect Hyades stars to live
# deliberately generous - the sigma clipping will clean it up afterwards
PM_BOX_PMRA_MIN   =  85.0   # mas/yr
PM_BOX_PMRA_MAX   = 125.0
PM_BOX_PMDEC_MIN  = -45.0
PM_BOX_PMDEC_MAX  = -10.0

# tighter parallax range now that we're in python (ADQL already filtered to 12-35 mas)
PARALLAX_BOX_MIN  = 16.0    # mas
PARALLAX_BOX_MAX  = 27.0    # mas

# sigma clipping settings - 3 sigma is standard, 10 iterations is way more than we'll need
SIGMA_CLIP_NSIGMA  = 3.0
SIGMA_CLIP_MAXITER = 10


def apply_quality_cuts(df: pd.DataFrame) -> pd.DataFrame:
    # throws away stars where the measurements aren't trustworthy
    # keeps anything that passes all the quality thresholds defined above
    n0 = len(df)

    # check the position/motion measurements are decent
    mask = (
        (df["parallax_over_error"] > MIN_PARALLAX_SNR)
        & (df["astrometric_excess_noise"] < MAX_ASTR_EXCESS_NOISE)
        & (df["astrometric_excess_noise_sig"] < MAX_EXCESS_NOISE_SIG)
        & (df["visibility_periods_used"] >= MIN_VISIBILITY_PERIODS)
    )

    # also check the brightness measurements are decent
    mask &= (
        (df["phot_g_mean_flux_over_error"] > MIN_G_FLUX_SNR)
        & (df["phot_bp_mean_flux_over_error"] > MIN_BP_FLUX_SNR)
        & (df["phot_rp_mean_flux_over_error"] > MIN_RP_FLUX_SNR)
        & (df["phot_bp_rp_excess_factor"] < MAX_BP_RP_EXCESS)
        & df["phot_bp_mean_mag"].notna()
        & df["phot_rp_mean_mag"].notna()
    )

    result = df[mask].copy()
    print(f"[membership] Quality cuts: {n0:,} -> {len(result):,} sources retained.")
    return result


def apply_pm_parallax_box(df: pd.DataFrame) -> pd.DataFrame:
    # keeps only stars that are moving in roughly the same direction as the Hyades
    # and are at roughly the right distance - cuts out most of the unrelated field stars
    n0 = len(df)
    mask = (
        (df["pmra"]     >= PM_BOX_PMRA_MIN)  & (df["pmra"]     <= PM_BOX_PMRA_MAX)
        & (df["pmdec"]  >= PM_BOX_PMDEC_MIN) & (df["pmdec"]    <= PM_BOX_PMDEC_MAX)
        & (df["parallax"] >= PARALLAX_BOX_MIN) & (df["parallax"] <= PARALLAX_BOX_MAX)
    )
    result = df[mask].copy()
    print(f"[membership] PM+parallax box:  {n0:,} -> {len(result):,} candidates.")
    return result


def sigma_clip_members(df: pd.DataFrame) -> pd.DataFrame:
    # iteratively removes stars that are too far from the cluster average
    # each round we recompute the centre and spread of the remaining stars,
    # kick out anything more than 3 sigma away, and repeat until nothing changes
    # uses median absolute deviation instead of std because it's less sensitive to outliers
    mask = np.ones(len(df), dtype=bool)
    cols = ["pmra", "pmdec", "parallax"]

    for iteration in range(SIGMA_CLIP_MAXITER):
        n_before = mask.sum()
        sub = df[mask]
        for col in cols:
            med = sub[col].median()
            # MAD * 1.4826 gives the same number as std would for a normal distribution
            # but it's not thrown off by the outliers we're trying to remove
            mad = (sub[col] - med).abs().median()
            sigma = mad * 1.4826 if mad > 0 else sub[col].std()
            lo = med - SIGMA_CLIP_NSIGMA * sigma
            hi = med + SIGMA_CLIP_NSIGMA * sigma
            # only look at stars still in the running
            idx = df.index[mask]
            drop = df.loc[idx, col].lt(lo) | df.loc[idx, col].gt(hi)
            mask[df.index.get_indexer(idx[drop])] = False

        n_after = mask.sum()
        print(f"[membership] Sigma-clip iter {iteration + 1}: {n_before} -> {n_after}")
        if n_after == n_before:
            break  # nothing got removed this round, we're done

    result = df[mask].copy()
    print(f"[membership] Final member count: {len(result):,}")
    return result


def compute_absolute_magnitude(df: pd.DataFrame) -> pd.DataFrame:
    # works out how bright each star actually is (not just how bright it looks from earth)
    # uses the parallax to figure out the distance, then corrects the apparent magnitude
    # the formula is: absolute mag = apparent mag + 5 + 5 * log10(parallax in arcsec)
    # since parallax here is in milliarcsec we divide by 1000 first
    df = df.copy()
    df["abs_g_mag"] = (
        df["phot_g_mean_mag"]
        + 5.0 * np.log10(df["parallax"] / 1000.0)
        + 5.0
    )
    return df


def select_members(raw_df: pd.DataFrame) -> pd.DataFrame:
    # runs the full pipeline: quality cuts, box filter, sigma clipping, then adds absolute magnitude
    df = apply_quality_cuts(raw_df)
    df = apply_pm_parallax_box(df)
    df = sigma_clip_members(df)
    df = compute_absolute_magnitude(df)
    return df
