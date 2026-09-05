#!/usr/bin/env python3
"""Snapshot the public TING ArcGIS layer using only the Python standard library."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.parse
import urllib.request

SERVICE = 'https://services8.arcgis.com/S9R3NgKp66dTIzOU/arcgis/rest/services/Ting_Outages_TN_Floods_4_2025_AOI_View/FeatureServer'


def request(url, **params):
    query = urllib.parse.urlencode({'f': 'json', **params})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url + '?' + query, timeout=90) as response:
                result = json.load(response)
            if 'error' in result:
                raise RuntimeError(f'ArcGIS error: {result["error"]}')
            return result
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def iso(milliseconds):
    if milliseconds is None:
        return None
    return datetime.fromtimestamp(milliseconds / 1000, timezone.utc).isoformat()


def fetch_batch(query, ids, oid):
    result = request(query, objectIds=','.join(map(str, ids)), outFields='*',
                     returnGeometry='true', outSR=4326)
    features = result.get('features', [])
    if result.get('exceededTransferLimit'):
        if len(ids) == 1:
            raise RuntimeError('Transfer limit exceeded for a single feature')
        mid = len(ids) // 2
        return fetch_batch(query, ids[:mid], oid) + fetch_batch(query, ids[mid:], oid)
    received = [f['attributes'][oid] for f in features]
    if len(received) != len(set(received)) or set(received) != set(ids):
        raise RuntimeError('Requested/received object IDs differ; snapshot is incomplete')
    return [result]


def profile(features, fields):
    rows = [f['attributes'] for f in features]
    dates = [f['name'] for f in fields if f['type'] == 'esriFieldTypeDate']
    ranges = {}
    for name in dates:
        values = [r[name] for r in rows if r.get(name) is not None]
        ranges[name] = {'min_utc': iso(min(values)) if values else None,
                        'max_utc': iso(max(values)) if values else None}
    events = Counter(r.get('eventId') for r in rows)
    return {
        'records': len(rows),
        'distinct_nonnull_event_ids': len([e for e in events if e is not None]),
        'event_ids_with_multiple_records': sum(n > 1 for e, n in events.items() if e is not None),
        'null_counts': {f['name']: sum(r.get(f['name']) is None for r in rows) for f in fields},
        'date_ranges': ranges,
        'power_states': dict(Counter(str(r.get('powerState')) for r in rows)),
        'county_labels': dict(Counter(str(r.get('County')) for r in rows)),
        'missing_geometry': sum(not f.get('geometry') for f in features),
        'modified_before_start': sum(r.get('modified') is not None and r.get('start') is not None
                                     and r['modified'] < r['start'] for r in rows),
        'caveats': [
            'Polygon records are not devices or necessarily distinct outages.',
            'Repeated event IDs may represent spatial parts or temporal updates; retained without deduplication.',
            'No monitored-device denominator: outage percentages and coverage confidence cannot be calculated.',
            'modified is not assumed to be restoration time; displaydate is not assumed to be overpass time.',
            'Full service extract is not a Dyer County spatial clip.',
            'Source Shape__Area is not assumed to be square meters.',
        ],
    }


def ingest(output, service=SERVICE, layer=0, batch_size=250):
    if batch_size < 1:
        raise ValueError('batch_size must be positive')
    output.mkdir(parents=True, exist_ok=False)
    raw = output / 'raw'
    raw.mkdir()
    started = datetime.now(timezone.utc).isoformat()
    manifest = {'status': 'incomplete', 'started_utc': started, 'service_url': service,
                'layer': layer, 'where': '1=1', 'output_spatial_reference': 4326}
    save(output / 'manifest.json', manifest)
    service_meta = request(service)
    url = f'{service.rstrip("/")}/{layer}'
    metadata = request(url)
    save(raw / 'service.json', service_meta)
    save(raw / 'layer.json', metadata)
    oid = metadata['objectIdField']
    query = url + '/query'
    id_response = request(query, where='1=1', returnIdsOnly='true')
    save(raw / 'object_ids.json', id_response)
    if id_response.get('exceededTransferLimit') or 'objectIds' not in id_response:
        raise RuntimeError('Server did not return a complete object ID inventory')
    ids = sorted(id_response['objectIds'] or [])
    if len(ids) != len(set(ids)):
        raise RuntimeError('Duplicate object IDs in inventory')
    count = request(query, where='1=1', returnCountOnly='true')
    save(raw / 'count.json', count)
    if count['count'] != len(ids):
        raise RuntimeError('Count differs from object ID inventory')
    size = min(batch_size, metadata.get('maxRecordCount', batch_size))
    features = []
    page = 0
    for offset in range(0, len(ids), size):
        for result in fetch_batch(query, ids[offset:offset + size], oid):
            save(raw / f'page_{page:05d}.json', result)
            features.extend(result['features'])
            page += 1
        print(f'Downloaded {len(features)}/{len(ids)} records', flush=True)
    final_ids = request(query, where='1=1', returnIdsOnly='true')
    save(raw / 'object_ids_end.json', final_ids)
    if final_ids.get('exceededTransferLimit') or 'objectIds' not in final_ids or sorted(final_ids['objectIds'] or []) != ids:
        raise RuntimeError('Source object ID inventory changed during download')
    features.sort(key=lambda f: f['attributes'][oid])
    fields = metadata['fields']
    dates = [f['name'] for f in fields if f['type'] == 'esriFieldTypeDate']
    names = [f['name'] for f in fields]
    with (output / 'attributes.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=names + [d + '_utc' for d in dates])
        writer.writeheader()
        for feature in features:
            row = dict(feature['attributes'])
            row.update({d + '_utc': iso(row.get(d)) for d in dates})
            writer.writerow(row)
    # ArcGIS rings retained verbatim: this is Esri JSON, not GeoJSON.
    save(output / 'features.esri.json', {'geometryType': metadata['geometryType'],
         'spatialReference': {'wkid': 4326}, 'fields': fields, 'features': features})
    save(output / 'profile.json', profile(features, fields))
    manifest.update(status='complete', completed_utc=datetime.now(timezone.utc).isoformat(),
                    record_count=len(features), pages=page,
                    consistency_note='ID inventories checked before/after; in-place source edits are not transactionally isolated.',
                    sha256={str(p.relative_to(output)): hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'manifest.json'})
    save(output / 'manifest.json', manifest)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='New snapshot directory; existing directories are never overwritten')
    parser.add_argument('--service', default=SERVICE)
    parser.add_argument('--layer', type=int, default=0)
    parser.add_argument('--batch-size', type=int, default=250)
    args = parser.parse_args()
    output = args.output or Path('data/ting') / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    print(ingest(output, args.service, args.layer, args.batch_size))


if __name__ == '__main__':
    main()
