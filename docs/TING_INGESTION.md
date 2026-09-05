# TING ingestion

## Run

From the repository root, using Python 3.10 or newer:

```bash
python3 scripts/ingest_ting.py
```

The default destination is a new timestamped directory under `data/ting/`. To name a snapshot:

```bash
python3 scripts/ingest_ting.py --output data/ting/initial
```

Network access to services8.arcgis.com is required. No API key or additional packages were needed for this public view during inspection. Existing output directories are never overwritten. If a run fails, retain its incomplete snapshot for diagnosis and rerun into a new directory. Downloads are excluded from Git.

## Source and extraction

[TING ArcGIS FeatureServer](https://services8.arcgis.com/S9R3NgKp66dTIzOU/arcgis/rest/services/Ting_Outages_TN_Floods_4_2025_AOI_View/FeatureServer), layer 0 (`WL Polygons Layer`).

The pipeline retrieves all records visible through this view with `where=1=1`. It saves service/layer metadata, obtains an object-ID inventory, checks the server count, and downloads bounded ID batches with all fields and polygon geometry in EPSG:4326. Transfer-limited batches are split; missing or duplicate IDs fail the run. The ID inventory is checked again at completion. This catches additions/deletions but cannot guarantee protection against in-place edits during a download; the service is not a transactionally isolated historical snapshot.

Transient network errors receive bounded retries. ArcGIS JSON errors fail explicitly. A manifest is marked complete only after extraction, validation, and output generation succeed.

## Outputs

| File | Contents |
| --- | --- |
| `raw/service.json`, `raw/layer.json` | Source metadata and schema |
| `raw/object_ids*.json`, `raw/count.json` | Completeness evidence |
| `raw/page_*.json` | Parsed source responses, including geometry; no record deduplication |
| `features.esri.json` | Combined Esri JSON polygons in EPSG:4326; not GeoJSON |
| `attributes.csv` | Original attributes plus ISO UTC columns for ArcGIS date fields; original epoch milliseconds retained |
| `profile.json` | Record/event counts, nulls, date ranges, states, county labels, basic checks and caveats |
| `manifest.json` | Source, extraction times, status, counts, and SHA-256 hashes |

## What these records support

The schema describes outage polygons with `eventId`, `start`, `modified`, `displaydate`, `powerState`, and `recoveredPercentage`. It also contains State/County/Zip labels and geometry area/length fields. It does **not** expose individual devices, eligible-device counts, or reporting-device counts.

- A polygon row is not necessarily a distinct outage. Repeated event IDs may reflect multiple spatial parts or updates. Preserve all rows until those semantics are established.
- `modified` is not yet established as restoration time. Do not subtract `start` from it and label the result outage duration.
- `displaydate` is the service's time field. Preserve it independently of event start and modification times.
- `recoveredPercentage` lacks a documented denominator here; do not equate it with the population or fraction of a 500 m cell with power.
- Source geometry bounds extend outside Tennessee. County labels may be null. The full extract is **not** a Dyer County clip, regardless of the service name.
- Geometry area fields are not assumed to be square meters. Spatial analysis needs appropriate projection and geometry validation.
- Missing outage polygons do not establish normal power or sensor coverage. Device-hour outage rates and device-density confidence cannot be derived from this schema alone.

## Next steps

Initial successful extract (`data/ting/initial`): 10,095 polygon records, 268 distinct event IDs, and 227 event IDs appearing in multiple records. All records have `powerState=Outage`; all State/County/Zip values are null. Actual `displaydate` values span March 28–April 17, 2025 (UTC), narrower than the advertised metadata time extent. No records have missing geometry or modification timestamps preceding start timestamps. These checks do not establish polygon validity or complete temporal coverage. No explicit recovered-state records appear in this extract, so restoration timing requires further investigation.

1. Review the profile and event histories to establish extraction coverage, duplication patterns, and timestamp semantics.
2. Obtain an authoritative Dyer County boundary, validate polygon geometry, and select/intersect polygons spatially while retaining original geometry and IDs. Do not rely solely on County labels.
3. Confirm polygon-generation rules and recovery semantics with the data provider.
4. Build daily and overpass-aligned summaries only for quantities supported by those semantics. Maintain unknown intervals and avoid treating absent records as zero outages.
5. Acquire coverage/denominator data if available before constructing outage fractions or calibrated confidence levels.
