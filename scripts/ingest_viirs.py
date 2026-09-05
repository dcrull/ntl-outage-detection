#!/usr/bin/env python3
"""Discover, cache and subset VNP46A2 Collection 2 for Dyer County.

Run from the repo root. No TING data are read. Credentials are never written to manifests.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import date, datetime, timedelta, timezone
import hashlib
from html.parser import HTMLParser
import json
import math
import os
from pathlib import Path
import re
import time
from urllib.parse import urljoin, urlparse

import h5py
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from shapely.geometry import shape
import requests

HOST = 'ladsweb.modaps.eosdis.nasa.gov'
ARCHIVE = f'https://{HOST}/api/v2/content/archives/allData/5200/VNP46A2'
GROUP = 'HDFEOS/GRIDS/VIIRS_Grid_DNB_2d/Data Fields'
LAYERS = {
    'radiance': 'DNB_BRDF-Corrected_NTL',
    'gap_filled_radiance': 'Gap_Filled_DNB_BRDF-Corrected_NTL',
    'lunar_irradiance': 'DNB_Lunar_Irradiance',
    'latest_high_quality_retrieval': 'Latest_High_Quality_Retrieval',
    'mandatory_quality': 'Mandatory_Quality_Flag',
    'cloud_mask': 'QF_Cloud_Mask',
    'snow_flag': 'Snow_Flag',
}
PATTERN = re.compile(r'VNP46A2\.A(\d{7})\.(h\d{2}v\d{2})\.002\.(\d{13})\.h5$')


class AuthenticationRequired(RuntimeError):
    pass


def now():
    return datetime.now(timezone.utc).isoformat()


def save_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def dates_between(start, end):
    if end < start:
        raise ValueError('End must be on or after start')
    return [start + timedelta(days=i) for i in range((end - start).days + 1)]


class ListingParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows = []
        self.cells = []
        self.href = None
        self.in_cell = False

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.cells, self.href = [], None
        elif tag == 'td':
            self.in_cell = True
            self.cells.append('')
        elif tag == 'a':
            self.href = dict(attrs).get('href')

    def handle_data(self, data):
        if self.in_cell:
            self.cells[-1] += data

    def handle_endtag(self, tag):
        if tag == 'td':
            self.in_cell = False
        elif tag == 'tr' and self.href and self.cells:
            self.rows.append((self.href, self.cells))


def get(url, token=None, stream=False):
    # Requests strips Authorization on cross-host redirects. Never attach tokens to other hosts.
    if urlparse(url).scheme != 'https' or urlparse(url).hostname != HOST:
        raise ValueError('Unexpected archive host')
    headers = {'User-Agent': 'ntl-outage-detection/0.1'}
    if token:
        headers['Authorization'] = 'Bearer ' + token
    for attempt in range(4):
        try:
            response = requests.get(url, headers=headers, timeout=(20, 120), stream=stream)
            if response.status_code in (401, 403) or urlparse(response.url).hostname == 'urs.earthdata.nasa.gov':
                response.close()
                raise AuthenticationRequired('NASA download authentication required; configure a valid token.')
            if response.status_code == 429 or response.status_code >= 500:
                response.close()
                if attempt == 3:
                    raise RuntimeError('NASA service unavailable after retries')
                time.sleep(2 ** attempt)
                continue
            response.raise_for_status()
            return response
        except (requests.ConnectionError, requests.Timeout):
            if attempt == 3:
                raise RuntimeError('NASA connection failed after retries') from None
            time.sleep(2 ** attempt)


def discover_day(day, tile, root, refresh=False):
    stamp = day.strftime('%Y%j')
    url = f'{ARCHIVE}/{day.year}/{day.timetuple().tm_yday:03d}'
    cache = root / 'listings' / f'{stamp}.json'
    if cache.exists() and not refresh:
        listing = json.loads(cache.read_text())
    else:
        with get(url) as response:
            html = response.text
        if '<table' not in html.lower():
            raise RuntimeError(f'Unrecognized archive listing for {day}')
        listing = {'url': url, 'retrieved_utc': now(), 'html': html}
        save_json(cache, listing)
    parser = ListingParser()
    parser.feed(listing['html'])
    matches = []
    for href, cells in parser.rows:
        file_url = urljoin(url + '/', href)
        filename = Path(urlparse(file_url).path).name
        match = PATTERN.fullmatch(filename)
        if match and match[1] == stamp and match[2] == tile:
            size = int(cells[2].strip()) if len(cells) >= 3 else None
            matches.append({'filename': filename, 'url': file_url, 'bytes': size,
                            'production_timestamp': match[3]})
    matches.sort(key=lambda m: m['production_timestamp'])
    return {'date': day.isoformat(), 'status': 'available' if matches else 'missing',
            'selected': matches[-1] if matches else None, 'candidates': matches}


def load_token(path):
    token = os.environ.get('LAADS_TOKEN') or os.environ.get('EARTHDATA_TOKEN')
    if token:
        return token.strip()
    if path.exists():
        if path.stat().st_mode & 0o077:
            raise ValueError('Token file must be private: chmod 600 ' + str(path))
        return path.read_text().strip() or None
    return None


def validate_hdf(path):
    with h5py.File(path, 'r') as f:
        group = f[GROUP]
        for name in LAYERS.values():
            if group[name].shape != (2400, 2400):
                raise ValueError('Unexpected A2 grid shape')


def download(item, root, token):
    target = root / 'raw' / item['filename']
    target.parent.mkdir(parents=True, exist_ok=True)
    receipt = target.with_suffix('.h5.json')
    if target.exists() and receipt.exists():
        record = json.loads(receipt.read_text())
        if (target.stat().st_size == record['bytes']
                and (item['bytes'] is None or target.stat().st_size == item['bytes'])
                and digest(target) == record['sha256']):
            validate_hdf(target)
            return target, record
        raise RuntimeError(f'Cache integrity failure for {target}; inspect before replacing')
    if target.exists():
        raise RuntimeError(f'File lacks integrity receipt: {target}; inspect before replacing')
    part = target.with_suffix('.h5.part')
    for attempt in range(3):
        try:
            with get(item['url'], token=token, stream=True) as response:
                with part.open('wb') as f:
                    first = True
                    for chunk in response.iter_content(1024 * 1024):
                        if first:
                            first = False
                            if not chunk.startswith(b'\x89HDF\r\n\x1a\n'):
                                raise AuthenticationRequired('Expected HDF5 but received another format; check NASA access.')
                        f.write(chunk)
            if item['bytes'] is not None and part.stat().st_size != item['bytes']:
                raise RuntimeError('Downloaded size differs from archive listing')
            validate_hdf(part)
            record = {'source_url': item['url'], 'downloaded_utc': now(),
                      'bytes': part.stat().st_size, 'sha256': digest(part)}
            part.replace(target)
            save_json(receipt, record)
            return target, record
        except AuthenticationRequired:
            part.unlink(missing_ok=True)
            raise
        except (requests.RequestException, OSError, RuntimeError, ValueError):
            part.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def crop_grid(boundary, tile):
    h, v = int(tile[1:3]), int(tile[4:6])
    west, north, step = -180 + h * 10, 90 - v * 10, 1 / 240
    xmin, ymin, xmax, ymax = boundary.bounds
    if not (west <= xmin <= xmax <= west + 10 and north - 10 <= ymin <= ymax <= north):
        raise ValueError('AOI exceeds configured tile; multi-tile ingestion is not implemented')
    c0, c1 = math.floor((xmin - west) / step), math.ceil((xmax - west) / step)
    r0, r1 = math.floor((north - ymax) / step), math.ceil((north - ymin) / step)
    affine = from_origin(west + c0 * step, north - r0 * step, step, step)
    return (slice(r0, r1), slice(c0, c1)), affine


def json_value(value):
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, (tuple, list)):
        return [json_value(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def cog(path, values, affine, nodata, scale=1.0, offset=0.0, **tags):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.stem + '.tmp.tif')
    with rasterio.open(temp, 'w', driver='COG', width=values.shape[1], height=values.shape[0],
                       count=1, dtype=values.dtype, crs='EPSG:4326', transform=affine,
                       nodata=nodata, compress='DEFLATE', blocksize=256,
                       overview_resampling='NEAREST') as dst:
        dst.write(values, 1)
        dst.scales = (scale,)
        dst.offsets = (offset,)
        dst.update_tags(**tags)
    temp.replace(path)


def extract(path, day, out, boundary, tile):
    slices, affine = crop_grid(boundary, tile)
    arrays, valid, attrs = {}, {}, {}
    with h5py.File(path, 'r') as source:
        group = source[GROUP]
        metadata = {'source_filename': path.name, 'date': day,
                    'global_attributes': {k: json_value(v) for k, v in source.attrs.items()},
                    'layer_attributes': attrs, 'row_start': slices[0].start,
                    'column_start': slices[1].start, 'crs': 'EPSG:4326',
                    'transform': list(affine)[:6], 'tile': tile}
        for key, name in LAYERS.items():
            ds = group[name]
            values = ds[slices]
            attrs[key] = {k: json_value(v) for k, v in ds.attrs.items()}
            fill = np.asarray(ds.attrs.get('_FillValue')).ravel()
            if len(fill) != 1 or fill[0] is None:
                raise ValueError(f'Missing or nonscalar fill metadata for {name}')
            fill = fill[0]
            ok = np.isfinite(values) & (values != fill)
            # Collection 2 stores radiance in physical units. Fail rather than silently
            # guessing an alternate scale/offset convention on changed source schemas.
            scale = float(np.asarray(ds.attrs.get('scale_factor', 1)).ravel()[0])
            offset = float(np.asarray(ds.attrs.get('add_offset', 0)).ravel()[0])
            expected_scale = 0.1 if key == 'lunar_irradiance' else 1.0
            if not math.isclose(scale, expected_scale, rel_tol=1e-6) or offset != 0:
                raise ValueError(f'Unexpected Collection 2 scale/offset for {name}: {scale}, {offset}')
            valid[key], arrays[key] = ok, values
            cog(out / f'{key}.tif', values, affine, fill, scale=scale, offset=offset,
                source_filename=path.name, source_layer=name, date=day)
    dims = arrays['radiance'].shape
    mask = geometry_mask([boundary.__geo_interface__], dims, affine, invert=True, all_touched=False)
    touched = geometry_mask([boundary.__geo_interface__], dims, affine, invert=True, all_touched=True)
    cog(out / 'county_mask.tif', mask.astype('uint8'), affine, None,
        definition='1 if cell center is inside county, 0 otherwise; boundary values are retained in source rasters')
    cog(out / 'county_intersects_mask.tif', touched.astype('uint8'), affine, None,
        definition='1 for cells touched by county, 0 otherwise')
    cloud = arrays['cloud_mask']
    cloud_class, quality = (cloud >> 6) & 3, (cloud >> 4) & 3
    # Descriptive screening only. Native source data remain available for alternative rules.
    usable = (valid['radiance'] & valid['mandatory_quality'] & valid['cloud_mask']
              & valid['snow_flag'] & (arrays['mandatory_quality'] == 0)
              & (cloud_class == 0) & (quality == 3) & ((cloud & 1) == 0)
              & (((cloud >> 9) & 1) == 0) & (arrays['snow_flag'] == 0))
    cog(out / 'strict_usable.tif', usable.astype('uint8'), affine, None,
        definition='Valid non-gap-filled radiance; mandatory QA=0; high-quality confidently clear nighttime mask; no IR cirrus; snow flag=0')
    metadata['sha256'] = {p.name: digest(p) for p in sorted(out.glob('*.tif'))}
    save_json(out / 'metadata.json', metadata)
    n = int(mask.sum())
    return {'date': day, 'status': 'complete', 'county_cells': n,
            'radiance_present_cells': int((mask & valid['radiance']).sum()),
            'confident_clear_cells': int((mask & valid['cloud_mask'] & (cloud_class == 0)).sum()),
            'strict_usable_cells': int((mask & usable).sum()),
            'strict_usable_percent': float(100 * (mask & usable).sum() / n),
            'source_filename': path.name}


def run(args):
    root = args.output
    root.mkdir(parents=True, exist_ok=True)
    boundary = shape(json.loads(args.boundary.read_text())['features'][0]['geometry'])
    if not boundary.is_valid:
        raise ValueError('County boundary is invalid')
    crop_grid(boundary, args.tile)
    days = dates_between(args.start, args.end)
    manifest = {'product': 'VNP46A2', 'collection': '002', 'archive_set': '5200',
                'tile': args.tile, 'start': str(args.start), 'end': str(args.end),
                'event_start': '2025-04-02', 'boundary': str(args.boundary),
                'boundary_sha256': digest(args.boundary), 'status': 'discovering', 'started_utc': now(),
                'note': 'Product dates are acquisition-day labels, not exact local overpass timestamps.'}
    save_json(root / 'manifest.json', manifest)
    def discover(day):
        return discover_day(day, args.tile, root, args.refresh)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        inventory = []
        for item in executor.map(discover, days):
            inventory.append(item)
            print(f'Discovered {len(inventory)}/{len(days)}: {item["date"]} {item["status"]}', flush=True)
    save_json(root / 'inventory.json', inventory)
    manifest.update(status='discovered', available_days=sum(i['selected'] is not None for i in inventory),
                    expected_days=len(days), total_source_bytes=sum(i['selected']['bytes'] or 0 for i in inventory if i['selected']))
    save_json(root / 'manifest.json', manifest)
    if args.discover_only:
        return
    token = load_token(args.token_file)
    summaries = []
    try:
        for item in inventory:
            day = item['date']
            if item['selected'] is None:
                summaries.append({'date': day, 'status': 'missing_source'})
                continue
            path, receipt = download(item['selected'], root, token)
            out = root / 'dyer' / day
            summary = extract(path, day, out, boundary, args.tile)
            summaries.append(summary)
            manifest.update(status='running', completed_days=len([s for s in summaries if s['status'] == 'complete']))
            save_json(root / 'manifest.json', manifest)
            print(f'{day}: {summary["strict_usable_percent"]:.1f}% of county cells pass strict screening', flush=True)
    except AuthenticationRequired:
        manifest.update(status='authentication_required')
        save_json(root / 'manifest.json', manifest)
        raise
    except Exception:
        manifest.update(status='incomplete')
        save_json(root / 'manifest.json', manifest)
        raise
    finally:
        save_json(root / 'coverage.json', summaries)
    fields = ['date', 'status', 'county_cells', 'radiance_present_cells', 'confident_clear_cells',
              'strict_usable_cells', 'strict_usable_percent', 'source_filename']
    with (root / 'coverage.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summaries)
    manifest.update(status='complete' if all(s['status'] == 'complete' for s in summaries) else 'complete_with_missing_dates',
                    finished_utc=now(), inventory_sha256=digest(root / 'inventory.json'))
    save_json(root / 'manifest.json', manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=date.fromisoformat, default=date(2025, 2, 5))
    parser.add_argument('--end', type=date.fromisoformat, default=date(2025, 4, 17))
    parser.add_argument('--tile', default='h09v05', choices=['h09v05'])
    parser.add_argument('--boundary', type=Path, default=Path('data/boundaries/dyer_county.geojson'))
    parser.add_argument('--output', type=Path, default=Path('data/viirs/VNP46A2.002'))
    parser.add_argument('--token-file', type=Path, default=Path('.secrets/earthdata_token'))
    parser.add_argument('--workers', type=int, choices=range(1, 5), default=4)
    parser.add_argument('--refresh', action='store_true', help='Refresh cached NASA listings (may select new production versions)')
    parser.add_argument('--discover-only', action='store_true')
    args = parser.parse_args()
    try:
        run(args)
    except AuthenticationRequired as exc:
        parser.exit(2, str(exc) + '\nSave a token in .secrets/earthdata_token (mode 600), then rerun.\n')


if __name__ == '__main__':
    main()
