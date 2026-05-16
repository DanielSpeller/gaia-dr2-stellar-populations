# this file handles talking to the Gaia website and downloading the star data
# we're looking at the Hyades cluster because it's close enough that its parallax
# and how fast it's moving across the sky are really obvious - easy to spot in the data
#
# the Gaia website can be a bit flakey so we try the fast method first,
# and if that times out we fall back to the slower queuing method

from pathlib import Path
import time

import pandas as pd
from astroquery.gaia import Gaia

# where the Hyades cluster is in the sky, and how big a circle to search around it
# the coordinates come from van Leeuwen 2009, adjusted to match the DR2 time reference
HYADES_RA_DEG     = 66.75   # right ascension in degrees
HYADES_DEC_DEG    = 16.87   # declination in degrees
SEARCH_RADIUS_DEG = 10.0    # search radius - big enough to catch the outer edges of the cluster

# we only want stars that could plausibly be at the Hyades distance (~46 pc)
# anything outside this parallax range is clearly too close or too far away
PARALLAX_MIN_MAS = 12.0   # anything smaller than this is too far away
PARALLAX_MAX_MAS = 35.0   # anything bigger than this is too close

RAW_OUTPUT = Path("data/raw/hyades_gaia_dr2_raw.csv")

# this is the actual query we send to the Gaia database
# ADQL is basically SQL - SELECT picks the columns, WHERE filters the rows
ADQL_QUERY = f"""
SELECT TOP 50000
    -- the star's position and how accurate the position measurement is
    source_id,
    ra, ra_error,
    dec, dec_error,
    -- parallax is how much the star appears to shift as earth goes around the sun
    -- bigger parallax = closer star. the Hyades should be around 21.5 mas
    parallax, parallax_error,
    parallax_over_error,               -- parallax divided by its error, basically a signal-to-noise
    -- proper motion = how fast the star is drifting across the sky year on year
    -- the Hyades moves really fast (~104 mas/yr sideways) which makes it easy to find
    pmra, pmra_error,
    pmdec, pmdec_error,
    -- these tell us how good the position/motion measurements are
    astrometric_excess_noise,          -- leftover noise after fitting the star's path - should be small
    astrometric_excess_noise_sig,      -- how significant that leftover noise is
    visibility_periods_used,           -- how many separate times Gaia observed this star
    astrometric_chi2_al,               -- how well the measurements fit the model
    astrometric_n_good_obs_al,         -- number of usable observations
    -- brightness measurements in three colour bands: G (broad), BP (blue), RP (red)
    phot_g_mean_mag,
    phot_bp_mean_mag,
    phot_rp_mean_mag,
    bp_rp,                             -- blue minus red colour, gaia already worked this out for us
    phot_g_mean_flux_over_error,       -- how clearly we can measure the brightness in G
    phot_bp_mean_flux_over_error,      -- same but blue band
    phot_rp_mean_flux_over_error,      -- same but red band
    phot_bp_rp_excess_factor,          -- if this is high the star might be two stars blended together
    -- how fast the star is moving toward/away from us (only available for some stars in DR2)
    radial_velocity,
    radial_velocity_error
FROM gaiadr2.gaia_source
WHERE
    -- only grab stars inside a circle centred on the Hyades
    CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {HYADES_RA_DEG}, {HYADES_DEC_DEG}, {SEARCH_RADIUS_DEG})
    ) = 1
    -- skip stars where gaia couldn't measure a parallax or proper motion at all
    AND parallax IS NOT NULL
    AND pmra     IS NOT NULL
    AND pmdec    IS NOT NULL
    -- rough distance filter - cuts out most of the unrelated background stars before we download
    AND parallax BETWEEN {PARALLAX_MIN_MAS} AND {PARALLAX_MAX_MAS}
    -- basic quality filters in the query itself so we don't download junk we'd just throw away
    AND parallax_over_error > 5            -- parallax needs to be at least 5x its own error
    AND visibility_periods_used >= 7       -- star needs to have been seen at least 7 separate times
    AND astrometric_excess_noise < 2.0     -- don't want stars with really messy position fits
    AND phot_g_mean_flux_over_error > 20   -- need a half-decent brightness measurement
"""


def _sync_query(verbose: bool = False) -> pd.DataFrame:
    # fast method - sends the query and waits right there for the answer
    # works fine for smaller result sets but the server sometimes just kills it
    Gaia.ROW_LIMIT = 50000
    job = Gaia.launch_job(ADQL_QUERY, verbose=verbose)
    return job.get_results().to_pandas()


def _async_query(verbose: bool = False) -> pd.DataFrame:
    # slower method - submits the query as a job and polls until it's done
    # more reliable for bigger queries, and survives temporary server hiccups
    # tries up to 3 times with increasing waits between attempts
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
                print(f"[query] Waiting {wait}s before retry ...")
                time.sleep(wait)
    raise RuntimeError("All async TAP attempts failed. Check ESA service status.")


def download_raw(output_path: Path = RAW_OUTPUT, overwrite: bool = False) -> pd.DataFrame:
    # downloads the star data from Gaia and saves it to a csv file
    # if the csv already exists we just load that instead of hitting the server again
    # pass overwrite=True if you want fresh data
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and not overwrite:
        print(f"[query] Raw cache found at {output_path} - loading from disk.")
        return pd.read_csv(output_path)

    print("[query] Submitting ADQL query to ESA Gaia TAP ...")
    print(f"[query] Cone: RA={HYADES_RA_DEG} deg, Dec={HYADES_DEC_DEG} deg, r={SEARCH_RADIUS_DEG} deg")

    df = None
    try:
        print("[query] Trying synchronous endpoint ...")
        df = _sync_query()
        print(f"[query] Synchronous query succeeded: {len(df):,} rows.")
    except Exception as exc:
        # sync timed out, try the async queue instead
        print(f"[query] Synchronous query failed ({exc}); switching to async ...")
        df = _async_query()
        print(f"[query] Async query succeeded: {len(df):,} rows.")

    df.to_csv(output_path, index=False)
    print(f"[query] Raw data saved to {output_path}")
    return df
