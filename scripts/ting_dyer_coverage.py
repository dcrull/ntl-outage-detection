"""Measure recorded TING polygon footprint, not sensor coverage or confirmed power status.

Run from repository root: .venv/bin/python scripts/ting_dyer_coverage.py
Requires shapely and pyproj in the project environment. Uses the saved Census boundary.
"""
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
from zoneinfo import ZoneInfo

from pyproj import Transformer
from shapely.geometry import Polygon, shape, mapping
from shapely.ops import transform, unary_union


def main():
    snapshot = Path('data/ting/initial')
    boundary_path = Path('data/boundaries/dyer_county.geojson')
    out = Path('data/ting/dyer_coverage')
    out.mkdir(parents=True, exist_ok=True)
    county_ll = shape(json.loads(boundary_path.read_text())['features'][0]['geometry'])
    assert county_ll.is_valid
    forward = Transformer.from_crs(4326, 5070, always_xy=True).transform
    reverse = Transformer.from_crs(5070, 4326, always_xy=True).transform
    county = transform(forward, county_ll)
    by_day = defaultdict(list)
    by_time = defaultdict(list)
    selected = []
    local = ZoneInfo('America/Chicago')
    features = json.loads((snapshot / 'features.esri.json').read_text())['features']
    for feature in features:
        rings = feature['geometry']['rings']
        # Verified for this extract; fail explicitly if future input needs hole/multipart handling.
        if len(rings) != 1:
            raise ValueError('This snapshot utility requires single-ring polygons')
        polygon = Polygon(rings[0])
        if not polygon.is_valid:
            raise ValueError('Invalid input polygon; investigate before calculating coverage')
        if not polygon.intersects(county_ll):
            continue
        clipped = transform(forward, polygon).intersection(county)
        if clipped.area == 0:
            continue
        attrs = feature['attributes']
        stamp = attrs['displaydate']
        day = datetime.fromtimestamp(stamp / 1000, timezone.utc).astimezone(local).date().isoformat()
        selected.append((attrs, clipped))
        by_day[day].append(clipped)
        by_time[stamp].append(clipped)
    footprint = unary_union([g for _, g in selected])
    def measure(geometries):
        area = unary_union(geometries).area
        return {'area_km2': area / 1e6, 'county_percent': 100 * area / county.area}
    daily = {day: measure(gs) for day, gs in sorted(by_day.items())}
    instant = {datetime.fromtimestamp(t / 1000, timezone.utc).isoformat(): measure(gs)
               for t, gs in sorted(by_time.items())}
    summary = {
        'metric': 'Union of recorded polygon footprints clipped to county; not confirmed outage prevalence',
        'source_snapshot': str(snapshot), 'boundary': str(boundary_path),
        'boundary_vintage': '2026-01-01', 'area_crs': 'EPSG:5070',
        'denominator': 'Full county polygon including water',
        'county_area_km2': county.area / 1e6,
        'intersecting_records': len(selected),
        'distinct_event_ids': len({a['eventId'] for a, _ in selected}),
        'all_extract_dates_union': measure([g for _, g in selected]),
        'outside_union_percent': 100 * (1 - footprint.area / county.area),
        'april_2_through_11_local_union': measure([g for day, gs in by_day.items()
                                                    if '2025-04-02' <= day <= '2025-04-11' for g in gs]),
        'daily_timezone': 'America/Chicago', 'daily_recorded_union': daily,
        'largest_same_displaydate_footprint': max(instant.items(), key=lambda item: item[1]['area_km2']) if instant else None,
        'caveats': [
            'All-time and daily unions do not imply simultaneous or continuous outage.',
            'Same-displaydate footprint also depends on unconfirmed source snapshot semantics.',
            'No record is not evidence of power availability or monitoring coverage.',
            'Land-area balance does not measure balance among sufficiently bright VIIRS cells.',
            'Days with no intersecting records are omitted, not labeled as zero outages.',
        ],
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (out / 'ting_dyer_clipped.geojson').write_text(json.dumps({'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': a, 'geometry': mapping(transform(reverse, g))}
        for a, g in selected]}))
    (out / 'ting_dyer_union.geojson').write_text(json.dumps({'type': 'FeatureCollection', 'features': [
        {'type': 'Feature', 'properties': summary['all_extract_dates_union'],
         'geometry': mapping(transform(reverse, footprint))}]}))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
