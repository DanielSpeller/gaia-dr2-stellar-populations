# gaia-dr2-stellar-populations

Downloads real star data from the Gaia space telescope and figures out which stars belong to the Hyades open cluster.

We use the Hyades rather than the Pleiades because the Hyades is close enough (~47 pc away) that its parallax is really obvious, and it moves across the sky unusually fast compared to background stars - so it's pretty easy to pick out which stars are actually in the cluster vs just happening to be in the same direction.

---

## how to run it

```bash
# 1. create a virtual environment and activate it
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2. install the dependencies
pip install -r requirements.txt

# 3. run the analysis (takes 1-3 min the first time while it downloads data)
python run_analysis.py

# if you want to re-download fresh data instead of using the cached version
python run_analysis.py --fresh
```

## what it produces

| file | what it is |
|---|---|
| `data/raw/hyades_gaia_dr2_raw.csv` | the raw download from Gaia (not committed to git, it's big) |
| `data/processed/members.csv` | the final list of cluster members |
| `figures/fig1_cmd.png` | colour-magnitude diagram (like an HR diagram) |
| `figures/fig2_proper_motion.png` | shows the cluster moving as a clump against background stars |
| `figures/fig3_parallax_histogram.png` | shows the cluster as a spike in the distance distribution |

---

## the query we send to Gaia

This is ADQL, which is basically SQL for astronomy databases. We send this to the Gaia archive and it sends back a table of stars.

```sql
SELECT TOP 50000
    source_id,
    ra, ra_error, dec, dec_error,
    parallax, parallax_error, parallax_over_error,
    pmra, pmra_error, pmdec, pmdec_error,
    astrometric_excess_noise, astrometric_excess_noise_sig,
    visibility_periods_used,
    astrometric_chi2_al, astrometric_n_good_obs_al,
    phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag,
    bp_rp,
    phot_g_mean_flux_over_error,
    phot_bp_mean_flux_over_error,
    phot_rp_mean_flux_over_error,
    phot_bp_rp_excess_factor,
    radial_velocity, radial_velocity_error
FROM gaiadr2.gaia_source
WHERE
    -- grab everything inside a 10 degree circle centred on the Hyades
    CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', 66.75, 16.87, 10.0)
    ) = 1
    -- skip stars where Gaia couldn't measure a parallax or proper motion
    AND parallax IS NOT NULL
    AND pmra     IS NOT NULL
    AND pmdec    IS NOT NULL
    -- rough distance filter to cut out most unrelated stars before downloading
    AND parallax BETWEEN 12.0 AND 35.0
    -- basic quality filters so we don't download junk we'd throw away anyway
    AND parallax_over_error > 5
    AND visibility_periods_used >= 7
    AND astrometric_excess_noise < 2.0
    AND phot_g_mean_flux_over_error > 20
```

---

## how we decide which stars are actually in the cluster

We do it in three steps:

### step 1 - throw away bad measurements

Some stars have noisy or unreliable measurements. We cut those out first.

| what we check | threshold | why |
|---|---|---|
| `parallax_over_error` | > 5 | parallax needs to be at least 5x its own error, otherwise it's basically noise |
| `astrometric_excess_noise` | < 1.0 mas | if there's this much leftover noise in the position fit, something's off |
| `astrometric_excess_noise_sig` | < 2.0 | how statistically significant that noise is |
| `visibility_periods_used` | >= 7 | needs to have been observed at least 7 separate times |
| `phot_g_mean_flux_over_error` | > 50 | brightness in the main band needs to be well measured |
| `phot_bp_mean_flux_over_error` | > 10 | same for blue and red bands |
| `phot_rp_mean_flux_over_error` | > 10 | |
| `phot_bp_rp_excess_factor` | < 1.6 | if this is high the "star" is probably two stars on top of each other |

### step 2 - keep only stars moving like the Hyades

The Hyades moves really fast across the sky compared to most stars. We draw a box around the expected proper motion and parallax of the cluster and throw everything outside it away.

| measurement | range we keep |
|---|---|
| `pmra` (sideways motion) | 85 to 125 mas/yr |
| `pmdec` (up/down motion) | -45 to -10 mas/yr |
| `parallax` (distance) | 16 to 27 mas |

### step 3 - sigma clipping to tighten it up

We repeatedly compute the average and spread of the remaining stars, remove anything more than 3 standard deviations away, and repeat until nothing gets removed. Usually takes about 5 rounds. This is a standard way to iteratively remove outliers.

### working out actual brightness

To plot the colour-magnitude diagram we need to know how bright each star actually is, not just how bright it looks from Earth. We use the parallax (which tells us the distance) to correct for that:

```
absolute magnitude = apparent magnitude + 5 + 5 * log10(parallax / 1000)
```

---

## files in this repo

```
.
├── run_analysis.py       # run this to do the whole analysis
├── requirements.txt      # python packages needed
├── src/
│   ├── query.py          # builds the query and downloads data from Gaia
│   ├── membership.py     # quality cuts and membership selection
│   └── plots.py          # makes the three figures
├── data/
│   ├── raw/              # downloaded data goes here (gitignored)
│   └── processed/
│       └── members.csv   # the final member list
└── figures/
    ├── fig1_cmd.png
    ├── fig2_proper_motion.png
    └── fig3_parallax_histogram.png
```

---

## papers this is based on

- Gaia Collaboration (2016) - overview of the Gaia mission. *A&A* 595, A1
- Gaia Collaboration (2018a) - the DR2 data release. *A&A* 616, A1
- Gaia Collaboration (2018b) - HR diagrams from DR2. *A&A* 616, A10
- van Leeuwen (2009) - Hyades distance from Hipparcos. *A&A* 497, 209
- Lindegren et al. (2018) - how the DR2 astrometry was done. *A&A* 616, A2
- Arenou et al. (2018) - quality checks on DR2. *A&A* 616, A17
