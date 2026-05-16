# gaia-dr2-stellar-populations

Queries the ESA Gaia DR2 archive via ADQL and identifies members of the
**Hyades open cluster** through astrometric and photometric selection.

The Hyades is chosen over the Pleiades because its large parallax (~21.5 mas,
d ≈ 46.5 pc) and extreme proper motion (µα* ≈ +104, µδ ≈ −28 mas yr⁻¹) make
membership identification unambiguous, and its ~10° sky footprint returns a
comfortably-sized catalogue without row-limit issues.

---

## Quickstart

```bash
# 1. Create and activate a virtual environment (Python ≥ 3.11)
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run the full pipeline (~1–3 min depending on TAP latency)
python run_analysis.py

# Force a fresh download from ESA even if a local cache exists
python run_analysis.py --fresh
```

### Outputs

| Path | Description |
|---|---|
| `data/raw/hyades_gaia_dr2_raw.csv` | Raw TAP result (gitignored) |
| `data/processed/members.csv` | Cleaned member catalogue (committed) |
| `figures/fig1_cmd.png` | Colour-magnitude diagram |
| `figures/fig2_proper_motion.png` | Proper-motion vector diagram |
| `figures/fig3_parallax_histogram.png` | Parallax histogram |

---

## ADQL Query

```sql
SELECT TOP 80000
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
    -- Spatial cone around the Hyades centre (van Leeuwen 2009, J2015.5)
    CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', 66.75, 16.87, 10.0)
    ) = 1
    -- Require 5-parameter astrometric solution
    AND parallax IS NOT NULL
    AND pmra     IS NOT NULL
    AND pmdec    IS NOT NULL
    -- Pre-filter: keep only sources consistent with Hyades distance (12–35 mas)
    AND parallax BETWEEN 12.0 AND 35.0
```

The `CONTAINS/POINT/CIRCLE` geometry functions follow the ADQL 2.0 standard
(IVOA 2008).  The coarse parallax window reduces download volume by ~60 %
before Python-side cuts.

---

## Membership Selection Criteria

### 1. Quality cuts (applied in Python)

| Cut | Threshold | Rationale |
|---|---|---|
| `parallax_over_error` | > 5 | Parallax relative error < 20 % |
| `astrometric_excess_noise` | < 1.0 mas | Residual astrometric noise; large values indicate unresolved binaries or bad solutions |
| `astrometric_excess_noise_sig` | < 2.0 | Significance of the excess noise |
| `visibility_periods_used` | ≥ 7 | Minimum independent-epoch constraint for a reliable solution |
| `phot_g_mean_flux_over_error` | > 50 | Photometric SNR in G |
| `phot_bp_mean_flux_over_error` | > 10 | Photometric SNR in BP |
| `phot_rp_mean_flux_over_error` | > 10 | Photometric SNR in RP |
| `phot_bp_rp_excess_factor` | < 1.6 | Flags blended sources or photometric artefacts |

### 2. Proper-motion and parallax box

A generous bounding box isolates the Hyades locus before sigma-clipping:

| Parameter | Range |
|---|---|
| `pmra` | 85 – 125 mas yr⁻¹ |
| `pmdec` | −45 – −10 mas yr⁻¹ |
| `parallax` | 16 – 27 mas |

### 3. Iterative 3σ clipping

Starting from the box-selected candidates, the median and MAD-based σ of
(pmra, pmdec, parallax) are recomputed each iteration and sources outside
3σ are rejected.  The loop converges when no further sources are removed
(typically 3–5 iterations).

### Absolute magnitude

Computed from the measured parallax as:

```
M_G = G + 5 + 5 × log₁₀(π / 1000)
    = G + 5 × log₁₀(π) − 10
```

where π is the parallax in mas.

---

## Repository structure

```
.
├── run_analysis.py          # Entry-point – runs the full pipeline
├── requirements.txt
├── src/
│   ├── query.py             # ADQL query & TAP download
│   ├── membership.py        # Quality cuts & membership selection
│   └── plots.py             # Three diagnostic figures
├── data/
│   ├── raw/                 # gitignored – large downloads
│   └── processed/
│       └── members.csv      # Committed output
└── figures/
    ├── fig1_cmd.png
    ├── fig2_proper_motion.png
    └── fig3_parallax_histogram.png
```

---

## References

- Gaia Collaboration et al. (2016) – Gaia mission description. *A&A* 595, A1.
- Gaia Collaboration et al. (2018a) – DR2 summary. *A&A* 616, A1.
- Gaia Collaboration et al. (2018b) – Hertzsprung-Russell diagram. *A&A* 616, A10.
- van Leeuwen (2009) – Hipparcos Hyades distance. *A&A* 497, 209.
- Lindegren et al. (2018) – DR2 astrometry. *A&A* 616, A2.
- Arenou et al. (2018) – DR2 catalogue validation. *A&A* 616, A17.
