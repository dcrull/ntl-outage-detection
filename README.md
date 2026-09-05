# ntl-outage-detection

Research on nighttime-light outage detection in Dyer County, Tennessee.

- [Draft analysis plan](ANALYSIS_PLAN_DRAFT.md)
- [TING ingestion and data caveats](docs/TING_INGESTION.md)

Download a new TING snapshot with Python 3.10+ (no third-party dependencies):

```bash
python3 scripts/ingest_ting.py
```

Run offline ingestion tests:

```bash
python3 -m unittest discover -s tests -v
```
