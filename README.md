# ntl-outage-detection

Research on nighttime-light outage detection in Dyer County, Tennessee.

- [Draft analysis plan](ANALYSIS_PLAN_DRAFT.md)
- [Counterfactual detector: literature review and proposed methods](docs/COUNTERFACTUAL_DETECTOR_PLAN.md)
- [M0/M1 implementation, validation, and first results](docs/COUNTERFACTUAL_IMPLEMENTATION.md)
- [M0/M1 on the retained 10 m reallocation](docs/COUNTERFACTUAL_10M.md)
- [Original reallocation: native-grid sanity check](docs/REALLOCATION_NATIVE_CHECK.md)
- [TING ingestion and data caveats](docs/TING_INGESTION.md)
- [VIIRS A2 ingestion, authentication, and outputs](docs/VIIRS_INGESTION.md)
- [QGIS: recover a stalled Mac data mount](docs/QGIS_MOUNT_RECOVERY.md)
- [Colleague discussion checkpoint](notebooks/06_discussion_checkpoint.ipynb)

Use this project's `.venv` for dependencies and analysis. VIIRS ingestion:

```bash
.venv/bin/python -m pip install -r requirements-viirs.txt
.venv/bin/python scripts/ingest_viirs.py --discover-only
.venv/bin/python scripts/ingest_viirs.py
```

Downloads require a NASA token; see the VIIRS instructions above.

Inspect the cached data with [02_inspect_viirs.ipynb](notebooks/02_inspect_viirs.ipynb).
Select the project's `.venv` kernel in VS Code and run all cells. Notebook dependencies
are listed in `requirements-notebooks.txt`. The notebook exports baseline reference,
variability, observation-count, and selected daily comparison COGs to
`data/published/viirs_dyer/` (on the Mac: `~/mnt/ntl/published/viirs_dyer/`).
These are exploratory radiance layers, not outage classifications.

Fit and inspect the initial non-gap-filled A2 counterfactual models:

```bash
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python scripts/run_counterfactual.py
```

Open [03_inspect_counterfactual.ipynb](notebooks/03_inspect_counterfactual.ipynb)
with the project's `.venv` kernel. Model diagnostics and QGIS COGs are published to
`data/published/counterfactual_dyer/` (Mac: `~/mnt/ntl/published/counterfactual_dyer/`).
These use no TING data and do not yet classify outages; predictive intervals need
calibration review, particularly for bright cells.

The original structural reallocation has also been aggregated to native VIIRS
cells and run through matched M0/M1 comparisons. Open
[04_inspect_reallocation_check.ipynb](notebooks/04_inspect_reallocation_check.ipynb)
and see the [reproduction instructions and results](docs/REALLOCATION_NATIVE_CHECK.md).
QGIS outputs are in `data/published/reallocation_check/`. The original operator
mixes neighboring radiance and does not reproduce each native A2 cell exactly.

For open-ended exploration, [05_explore_nightly_radiance.ipynb](notebooks/05_explore_nightly_radiance.ipynb)
loads A2 and aggregated reallocation as time × height × width NumPy cubes and
nightly pandas dataframes with integer row/column indices. It also opens the
retained 10 m reallocation cube as a read-only memory map (~6.5 GB on disk).
Each input preserves its missing observations. It performs no modeling or
clustering. Use the project's `.venv` kernel with `requirements-notebooks.txt`
installed. To reconstruct the fine scenes from the cached original inputs, run
`.venv/bin/python scripts/export_reallocation_fine.py` with
`requirements-reallocation.txt` installed and the original reallocation checkout
available. Every nightly fine scene is verified against its published native
aggregate before the export is marked complete.

For a visual discussion of all six dataset/model combinations, open
[06_discussion_checkpoint.ipynb](notebooks/06_discussion_checkpoint.ipynb).
It includes common-area imagery, spatial-median diagnostics, missing-aware
residual clustering, and Dyersburg close-ups. Use `requirements-checkpoint.txt`
in the project `.venv`. Shareable figures, a PDF, cluster rasters, tables, and the
exported HTML are in `data/published/discussion_checkpoint/`; the notebook's first
section records the comparison masks, clustering choices, and limitations.
Radiance maps use shared log1p color normalization, and residual maps use a
zero-centered symmetric-log scale (linear within ±0.5 radiance units by default).
The display controls include automatic full-data limits and manual overrides;
the underlying values and statistical calculations stay in their original units.

Regenerate the executed discussion notebook and shareable outputs with:

```bash
.venv/bin/python -m pip install -r requirements-checkpoint.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/render_discussion_checkpoint.py
```

Download a new TING snapshot with Python 3.10+ (no third-party dependencies):

```bash
python3 scripts/ingest_ting.py
```

Run offline ingestion tests:

```bash
python3 -m unittest discover -s tests -v
```
