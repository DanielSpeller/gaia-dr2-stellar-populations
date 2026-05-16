"""
query.py – Download Gaia DR2 data for the Hyades open cluster.

We use the Hyades (not the Pleiades) because its large parallax (~21.5 mas,
d ≈ 46.5 pc) and extreme proper motion (pmra ≈ +104, pmdec ≈ -28 mas/yr)
make membership selection unambiguous and the cone search comfortable at
10-degree radius without exceeding reasonable row counts.
"""

from pathlib import Path

import pandas as pd
from astroquery.gaia import Gaia

# ── Cluster parameters ────────────────────────────────────────────────────────
# Hyades centre from van Leeuwen (2009), propagated to J2015.5 (DR2 epoch)
HYADES_RA_DEG   = 66.75   # degrees
HYADES_DEC_DEG  = 16.87   # degrees
SEARCH_RADIUS_DEG = 10.0  # captures core + first-generation tidal tails

# Wide parallax window used in the ADQL itself to reduce download volume;
# refined cuts are applied in Python after download.
PARALLAX_MIN_MAS = 12.0   # mas  (< 83 pc – well below Hyades distance)
PARALLAX_MAX_MAS = 35.0   # mas  (> 29 pc – comfortably above background)

RAW_OUTPUT = Path("data/raw/hyades_gaia_dr2_raw.csv")

# ── ADQL query ────────────────────────────────────────────────────────────────
# Columns are ordered for readability; comments explain each clause.
ADQL_QUERY = f"""
SELECT TOP 80000
    -- astrometry
    source_id,
    ra, ra_error,
    dec, dec_error,
    parallax, parallax_error,
    parallax_over_error,               -- SNR; we require > 5 after download
    pmra, pmra_error,
    pmdec, pmdec_error,
    -- astrometric quality
    astrometric_excess_noise,          -- residual noise in astrometric solution
    astrometric_excess_noise_sig,      -- significance of excess noise
    visibility_periods_used,           -- independent transits; require > 6
    astrometric_chi2_al,
    astrometric_n_good_obs_al,
    -- photometry
    phot_g_mean_mag,
    phot_bp_mean_mag,
    phot_rp_mean_mag,
    bp_rp,
    phot_g_mean_flux_over_error,       -- photometric SNR in G
    phot_bp_mean_flux_over_error,      -- photometric SNR in BP
    phot_rp_mean_flux_over_error,      -- photometric SNR in RP
    phot_bp_rp_excess_factor,          -- BP+RP flux vs G; flags blends/artefacts
    -- radial velocity (not always available in DR2)
    radial_velocity,
    radial_velocity_error
FROM gaiadr2.gaia_source
WHERE
    -- Spatial cone around the Hyades centre
    CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {HYADES_RA_DEG}, {HYADES_DEC_DEG}, {SEARCH_RADIUS_DEG})
    ) = 1
    -- Require valid astrometry (NULL parallax means no 5-param solution)
    AND parallax IS NOT NULL
    AND pmra    IS NOT NULL
    AND pmdec   IS NOT NULL
    -- Coarse parallax window to discard obvious foreground/background
    -- This cuts the download size by ~60 % before quality filtering
    AND parallax BETWEEN {PARALLAX_MIN_MAS} AND {PARALLAX_MAX_MAS}
"""


def download_raw(output_path: Path = RAW_OUTPUT, overwrite: bool = False) -> pd.DataFrame:
    """
    Submit the ADQL query to the Gaia TAP service and return a DataFrame.

    The raw result is saved to *output_path* (in data/raw/, which is gitignored)
    so subsequent runs can skip the network request.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not overwrite:
        print(f"[query] Raw cache found at {output_path} – loading from disk.")
        return pd.read_csv(output_path)

    print("[query] Submitting ADQL query to ESA Gaia TAP …")
    print(f"[query] Cone: RA={HYADES_RA_DEG}°, Dec={HYADES_DEC_DEG}°, r={SEARCH_RADIUS_DEG}°")

    # Use async job so the TAP server doesn't time-out on large results
    Gaia.ROW_LIMIT = -1  # lift the default 2000-row cap for async jobs
    job = Gaia.launch_job_async(ADQL_QUERY, dump_to_file=False, verbose=False)
    results = job.get_results()

    df = results.to_pandas()
    print(f"[query] Retrieved {len(df):,} rows from Gaia DR2.")

    df.to_csv(output_path, index=False)
    print(f"[query] Raw data saved to {output_path}")
    return df
