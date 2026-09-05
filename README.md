# ntl-outage-detection

Research on nighttime-light outage detection in Dyer County, Tennessee.

- [Draft analysis plan](ANALYSIS_PLAN_DRAFT.md)
- [TING ingestion and data caveats](docs/TING_INGESTION.md)
- [VIIRS A2 ingestion, authentication, and outputs](docs/VIIRS_INGESTION.md)

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

Download a new TING snapshot with Python 3.10+ (no third-party dependencies):

```bash
python3 scripts/ingest_ting.py
```

Run offline ingestion tests:

```bash
python3 -m unittest discover -s tests -v
```
