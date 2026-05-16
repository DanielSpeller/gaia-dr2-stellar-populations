"""
query.py – Download Gaia DR2 data for the Hyades open cluster.

We use the Hyades (not the Pleiades) because its large parallax (~21.5 mas,
d ≈ 46.5 pc) and extreme proper motion (pmra ≈ +104, pmdec ≈ -28 mas/yr)
make membership selection unambiguous and the cone search comfortable at
10-degree radius without exceeding reasonable row counts.

Note on query strategy
----------------------
The ESA Gaia TAP supports both synchronous (fast, ≤ a few thousand rows) and
asynchronous jobs.  We push quality cuts into the ADQL WHERE clause so that
the synchronous endpoint can handle the result comfortably without timing out.
A retry with the async endpoint is attempted automatically if sync fails.
"""

from pathlib import Path
import time

import pandas as pd
from astroquery.gaia import Gaia

# ── Cluster parameters ────────────────────────────────────────────────────────
# Hyades centre from van Leeuwen (2009), propagated to J2015.5 (DR2 epoch)
HYADES_RA_DEG     = 66.75   # degrees
HYADES_DEC_DEG    = 16.87   # degrees
SEARCH_RADIUS_DEG = 10.0    # captures core + first-generation tidal tails

# Parallax window applied in ADQL (Hyades ≈ 21.5 mas; field stars are mostly
# much nearer or much further, so this still leaves interesting context stars)
PARALLAX_MIN_MAS = 12.0   # mas  (~83 pc upper distance limit)
PARALLAX_MAX_MAS = 35.0   # mas  (~29 pc lower distance limit)

RAW_OUTPUT = Path("data/raw/hyades_gaia_dr2_raw.csv")

# ── ADQL query ────────────────────────────────────────────────────────────────
# Quality cuts are pushed into WHERE so the result fits in a synchronous job.
# Inline comments explain every clause.
ADQL_QUERY = f"""
SELECT TOP 50000
    -- ── astrometry ───────────────────────────────────────────────────────
    source_id,
    ra, ra_error,
    dec, dec_error,
    parallax, parallax_error,
    parallax_over_error,               -- SNR on parallax; require > 5
    pmra, pmra_error,
    pmdec, pmdec_error,
    -- ── astrometric quality flags ─────────────────────────────────────
    astrometric_excess_noise,          -- residual noise after best-fit model
    astrometric_excess_noise_sig,      -- how significant the residual is
    visibility_periods_used,           -- independent observing windows
    astrometric_chi2_al,               -- chi-squared of along-scan residuals
    astrometric_n_good_obs_al,         -- number of good along-scan obs
    -- ── photometry ───────────────────────────────────────────────────────
    phot_g_mean_mag,
    phot_bp_mean_mag,
    phot_rp_mean_mag,
    bp_rp,                             -- pre-computed BP-RP colour index
    phot_g_mean_flux_over_error,       -- G-band photometric SNR
    phot_bp_mean_flux_over_error,      -- BP-band photometric SNR
    phot_rp_mean_flux_over_error,      -- RP-band photometric SNR
    phot_bp_rp_excess_factor,          -- (flux_BP+flux_RP)/flux_G; flags blends
    -- ── radial velocity (sparse in DR2) ──────────────────────────────
    radial_velocity,
    radial_velocity_error
FROM gaiadr2.gaia_source
WHERE
    -- Spatial cone centred on the Hyades (ADQL 2.0 geometry functions)
    CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {HYADES_RA_DEG}, {HYADES_DEC_DEG}, {SEARCH_RADIUS_DEG})
    ) = 1
    -- Require a 5-parameter astrometric solution (NULL ↔ no parallax fit)
    AND parallax IS NOT NULL
    AND pmra     IS NOT NULL
    AND pmdec    IS NOT NULL
    -- Coarse parallax window: keeps stars at Hyades distance + context field
    -- (~60 % row reduction vs no parallax cut, speeds up synchronous response)
    AND parallax BETWEEN {PARALLAX_MIN_MAS} AND {PARALLAX_MAX_MAS}
    -- Push basic quality cuts into ADQL to shrink the download further
    AND parallax_over_error > 5            -- relative parallax error < 20 %
    AND visibility_periods_used >= 7       -- at least 7 independent epochs
    AND astrometric_excess_noise < 2.0     -- tolerate slight excess at ADQL level
    AND phot_g_mean_flux_over_error > 20   -- require at least basic G-band SNR
"""


def _sync_query(verbose: bool = False) -> pd.DataFrame:
    """Run as a synchronous TAP job (fast, works for small result sets)."""
    Gaia.ROW_LIMIT = 50000
    job = Gaia.launch_job(ADQL_QUERY, verbose=verbose)
    return job.get_results().to_pandas()


def _async_query(verbose: bool = False) -> pd.DataFrame:
    """Run as an asynchronous TAP job with retries (robust for larger results)."""
    Gaia.ROW_LIMIT = -1
    max_attempts = 3
    for attempt in range(1, max_attempts + 1):
        try:
            job = Gaia.launch_job_async(ADQL_QUERY, dump_to_file=False, verbose=verbose)
            return job.get_results().to_pandas()
        except Exception as exc:
            print(f"[query] Async attempt {attempt}/{max_attempts} failed: {exc}")
            if attempt < max_attempts:
                wait = 15 * attempt
                print(f"[query] Waiting {wait}s before retry …")
                time.sleep(wait)
    raise RuntimeError("All async TAP attempts failed. Check ESA service status.")


def download_raw(output_path: Path = RAW_OUTPUT, overwrite: bool = False) -> pd.DataFrame:
    """
    Submit the ADQL query to the ESA Gaia TAP service and return a DataFrame.

    Tries the synchronous endpoint first (lower overhead); falls back to the
    asynchronous endpoint with retries if the server returns an error.

    The raw result is cached to *output_path* (gitignored data/raw/) so that
    subsequent runs skip the network entirely.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not overwrite:
        print(f"[query] Raw cache found at {output_path} - loading from disk.")
        return pd.read_csv(output_path)

    print("[query] Submitting ADQL query to ESA Gaia TAP …")
    print(f"[query] Cone: RA={HYADES_RA_DEG} deg, Dec={HYADES_DEC_DEG} deg, r={SEARCH_RADIUS_DEG} deg")

    df = None
    try:
        print("[query] Trying synchronous endpoint …")
        df = _sync_query()
        print(f"[query] Synchronous query succeeded: {len(df):,} rows.")
    except Exception as exc:
        print(f"[query] Synchronous query failed ({exc}); switching to async …")
        df = _async_query()
        print(f"[query] Async query succeeded: {len(df):,} rows.")

    df.to_csv(output_path, index=False)
    print(f"[query] Raw data saved to {output_path}")
    return df
