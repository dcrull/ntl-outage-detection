# Dyer County TING footprint check

Calculated from `data/ting/initial/features.esri.json` using
`scripts/ting_dyer_coverage.py`. Requires project-local `shapely` and `pyproj`.

Boundary: [US Census TIGERweb Counties](https://tigerweb.geo.census.gov/arcgis/rest/services/TIGERweb/State_County/MapServer/1),
January 1, 2026 vintage, GEOID 47045, saved in
`data/boundaries/dyer_county.geojson` with query provenance alongside it.
This is a preliminary check using a current boundary, not an event-vintage boundary analysis.

All source polygons were valid single-ring geometries. Projected into EPSG:5070
(equal area), clipped to the county, and dissolved across records and dates to
avoid double counting overlaps. Denominator includes county land and water.

| Quantity | Result |
| --- | ---: |
| Full county area | 1,362.82 km² |
| Union of recorded TING polygons inside county | 201.72 km² |
| Share inside footprint | 14.80% |
| Share outside footprint | 85.20% |
| Intersecting source records | 359 |
| Distinct intersecting event IDs | 5 |

Daily union by `displaydate`, America/Chicago: April 2, 2025 = 14.78%;
April 3 = 4.02%. These are daily footprint unions, not daily mean outage area
or continuous outage duration. All five events start on April 2 local time;
the last display timestamp is April 3 at 17:30 CDT. No other dates in the
extract have positive-area intersections with Dyer County, including the
April 7–8 flood period. Absence is not evidence of restoration or no outages.

The 85.20% outside the footprint is not confirmed normal power and does not
measure TING sensor coverage. The geometric balance supports a Dyer-only pilot,
but adequacy depends on usable nighttime-light observations and how many
sufficiently illuminated cells lie inside/outside the polygons. This AOI
check must not be used to tune the detector with TING labels.

Saved outputs in `data/ting/dyer_coverage/`: `summary.json`,
`ting_dyer_clipped.geojson` (all intersecting records), and
`ting_dyer_union.geojson` (combined footprint). The GeoJSON files can be opened
in local QGIS through the mounted data directory.
