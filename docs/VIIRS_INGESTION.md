# VIIRS A2 ingestion for Dyer County

## Scope

- Product: Suomi NPP `VNP46A2`, Collection `002`, NASA archive set `5200`.
- Tile: `h09v05`, 15 arc-second native geographic grid (2400 × 2400 cells per 10° tile).
- Baseline: February 5–April 1, 2025 (56 product dates).
- Event/recovery observation window: April 2–17, 2025 (16 product dates).
- AOI: saved Census Dyer County boundary. This is the January 2026 boundary used in the preliminary TING coverage check.
- No TING labels, geometries, or other observations are used by this pipeline.

Date labels represent the product acquisition day. They are not exact overpass timestamps or local calendar days. Before matching TING, establish actual observation times from appropriate ancillary data. The baseline/event split should then be checked around the April 2 local-evening onset. The April 17 endpoint is an observation limit, not an assumed recovery date.

## Environment and access

Use the project environment:

```bash
.venv/bin/python -m pip install -r requirements-viirs.txt
```

Discovery uses public NASA listings. HDF5 downloads require an Earthdata or LAADS download token. Obtain one through [NASA LAADS](https://ladsweb.modaps.eosdis.nasa.gov/) or your [Earthdata account](https://urs.earthdata.nasa.gov/). Do not paste credentials in notebooks, chat, command arguments, or tracked files.

To save the token from an interactive Bash terminal on the server:

```bash
mkdir -p .secrets
chmod 700 .secrets
read -rs -p 'NASA download token: ' viirs_token
(umask 077; printf '%s' "$viirs_token" > .secrets/earthdata_token)
unset viirs_token
chmod 600 .secrets/earthdata_token
```

The token file must contain just the token, without quotes or a `Bearer` prefix.
The directory is ignored by Git. Alternatively, the downloader accepts an existing
`LAADS_TOKEN` or `EARTHDATA_TOKEN` environment variable (in that priority order).
File permissions are checked when using a token file. Tokens are only attached to
HTTPS requests to NASA's LAADS host; cross-host redirects do not retain the authorization header.

## Run

```bash
# Inventory the dates and sizes without credentials or data downloads.
.venv/bin/python scripts/ingest_viirs.py --discover-only

# Download, validate, and extract all 72 dates.
.venv/bin/python scripts/ingest_viirs.py
```

Use `--start YYYY-MM-DD --end YYYY-MM-DD` to change the window (both inclusive).
Use a separate `--output` directory for a separate experiment: manifests and coverage
reports describe the most recent run in that directory. Extraction files for older
dates are retained, so use `inventory.json` to determine the active run's dates.
Only the configured Dyer tile is currently supported; an AOI crossing its edge fails explicitly.

Discovery caches directory responses and selects the newest production timestamp
among matching files for each day. All candidates are recorded. `--refresh` checks
NASA again for reprocessed versions; without it, the cached listings freeze the
selection. Record which production versions were used in any subsequent analysis.

Verified raw downloads are reused on subsequent runs. Interrupted downloads are
retried from the beginning of that file; completed files are not redownloaded.
Files are checked against archive byte counts, local SHA-256 receipts, HDF5 format,
and required dataset shapes. Checksums detect local changes; they are not independent
NASA-published checksums. A corrupt cached file fails explicitly for investigation.

Extraction is rerun from cached sources, so changes in code or screening rules are
reflected without another download. Authentication errors and partial runs are
recorded in the manifest. Files unavailable in a valid listing remain missing dates,
not zero-radiance images. Network/listing errors fail instead of being labeled missing.

## Outputs

Under `data/viirs/VNP46A2.002/`:

| Path | Content |
| --- | --- |
| `listings/` | Public source listings with retrieval times |
| `inventory.json` | Selected files, alternate production versions, sizes, and missing dates |
| `raw/*.h5` | Full original NASA tiles |
| `raw/*.h5.json` | Source URL, download time, byte count, SHA-256 |
| `dyer/YYYY-MM-DD/*.tif` | Losslessly compressed, native-grid county bounding-box COGs |
| `dyer/YYYY-MM-DD/metadata.json` | Source attributes, crop offsets, transform, and output hashes |
| `coverage.json`, `coverage.csv` | Counts and percentages of cells passing descriptive quality checks |
| `manifest.json` | Configuration, boundary hash, source size, progress, and completion status |

Each date contains corrected `radiance`, `gap_filled_radiance`, `lunar_irradiance`,
`latest_high_quality_retrieval`, `mandatory_quality`, `cloud_mask`, and `snow_flag`
rasters. Source fill values and QA integer bit fields are retained. The reader expects
Collection 2 physical radiance values (scale 1, offset 0). Lunar irradiance retains
its source encoding with scale 0.1, recorded in the TIFF band scale; apply that scale
when reading raw array values for analysis. Unexpected scaling fails explicitly.
Pixel values are not resampled.

`county_mask.tif` identifies cells whose centers are inside the county;
`county_intersects_mask.tif` additionally includes boundary-touching cells. All crop
values remain available outside those masks. Coverage percentages use the center-based
mask, not fractional county area, population, device counts, or illuminated cells only.

`strict_usable.tif` is an initial descriptive screening rule, **not an outage detector**:

- Non-fill, finite corrected radiance (zero radiance remains valid).
- Mandatory quality flag = 0.
- Cloud detection = confidently clear; cloud-mask quality = high.
- Nighttime flag; no infrared cirrus flag; snow flag = 0.

All native fields remain available for alternative screening. Gap-filled radiance is
never substituted into this strict-usable mask. The exact rule may be refined after
reviewing observation quality, without using TING to tune the outage detector.

On the Mac, the mounted outputs appear under
`~/mnt/ntl/viirs/VNP46A2.002/dyer/`. Open `radiance.tif` with its corresponding mask
and QA layers in QGIS. Keep source and generated data out of Git.

## Validation and references

Initial live ingestion completed for all 72 dates (2,161,285,697 source bytes).
All raw-file receipts and output hashes were verified; every extracted scientific
layer was compared pixel-for-pixel with its source HDF5 crop. Source dates and
tile bounds were also checked. Each crop is 81 rows × 139 columns, with 7,868
cell centers inside the county.

Under the initial strict screening rule, the median county cell has 26 usable
baseline nights and 8 event-window nights. April 2–6 have no county cells passing
that rule. April 7 has 65.8% passing, April 8 has 92.9%. These are product-date
labels and screening results, not inferred outage states or attribution of every
rejection to cloud. The missing early-event observations are a major limitation
for detecting the short TING events on April 2–3.

The post-ingestion audit is saved in `validation_summary.json`; supplemental
`baseline_usable_nights.tif` and `event_usable_nights.tif` maps count passing
observations per cell. These audit outputs describe this completed run and should
be regenerated if the extraction or screening changes.

```bash
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python -m pip check
```

Tests cover date/grid selection, production-version selection, credential handling,
rejection of login HTML, and COG extraction with valid zeros, missing values, and
cloud contamination. A successful live run additionally verifies each source file.

- [NASA VNP46A2 Collection 2 description](https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/VNP46A2/)
- [NASA Collection 2 user guide](https://landweb.modaps.eosdis.nasa.gov/data/userguide/BlackMarbleUserGuide_Collection2.0_20241203.pdf)
- [LAADS API v2](https://ladsweb.modaps.eosdis.nasa.gov/tools-and-services/)
