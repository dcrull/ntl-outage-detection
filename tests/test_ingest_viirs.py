from datetime import date
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import h5py
import numpy as np
import rasterio
from shapely.geometry import box

from scripts.ingest_viirs import (AuthenticationRequired, GROUP, LAYERS, crop_grid,
                                 dates_between, discover_day, download, extract, get, load_token)


class VIIRSTests(unittest.TestCase):
    def test_dates_and_native_grid(self):
        self.assertEqual(len(dates_between(date(2025, 2, 5), date(2025, 4, 17))), 72)
        slices, affine = crop_grid(box(-89.5, 36, -89.49, 36.01), 'h09v05')
        self.assertAlmostEqual(affine.a, 1 / 240)
        self.assertLess(affine.e, 0)
        self.assertLessEqual(affine.c, -89.5)
        self.assertGreaterEqual(affine.f, 36.01)
        with self.assertRaises(ValueError):
            crop_grid(box(-91, 35, -89, 36), 'h09v05')

    def test_discovery_selects_latest_production_and_filters_tile(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'listings').mkdir()
            names = ['VNP46A2.A2025092.h09v05.002.2025101024806.h5',
                     'VNP46A2.A2025092.h09v05.002.2026101024806.h5',
                     'VNP46A2.A2025092.h10v05.002.2026101024806.h5']
            html = '<table>' + ''.join(f'<tr><td><a href="{n}">{n}</a></td><td>File</td><td>123</td></tr>' for n in names) + '</table>'
            (root / 'listings/2025092.json').write_text(json.dumps({'html': html}))
            result = discover_day(date(2025, 4, 2), 'h09v05', root)
            self.assertEqual(result['selected']['filename'], names[1])
            self.assertEqual(len(result['candidates']), 2)

    def test_credentials_and_host(self):
        with tempfile.TemporaryDirectory() as tmp, patch.dict('os.environ', {}, clear=True):
            p = Path(tmp) / 'token'
            self.assertIsNone(load_token(p))
            p.write_text('test-token\n')
            p.chmod(0o644)
            with self.assertRaises(ValueError):
                load_token(p)
            p.chmod(0o600)
            self.assertEqual(load_token(p), 'test-token')
        with self.assertRaises(ValueError):
            get('https://example.com/file', token='test-token')

    def test_login_html_never_becomes_cached_hdf(self):
        with tempfile.TemporaryDirectory() as tmp, patch('scripts.ingest_viirs.get') as mocked:
            mocked.return_value.__enter__.return_value.iter_content.return_value = [b'<html>login</html>']
            with self.assertRaises(AuthenticationRequired):
                download({'filename': 'sample.h5', 'url': 'https://example.com', 'bytes': 18}, Path(tmp), None)
            self.assertFalse((Path(tmp) / 'raw/sample.h5').exists())
            self.assertFalse((Path(tmp) / 'raw/sample.h5.part').exists())

    def test_extract_preserves_dark_valid_values_and_quality(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / 'sample.h5'
            # Exactly two native cells within tile h09v05.
            boundary = box(-90 + 100.01 / 240, 40 - 101.99 / 240,
                           -90 + 101.99 / 240, 40 - 100.01 / 240)
            with h5py.File(src, 'w') as f:
                group = f.create_group(GROUP)
                for key, name in LAYERS.items():
                    dtype = 'float32' if 'radiance' in key else 'uint16' if key == 'cloud_mask' else 'uint8'
                    fill = 65535 if dtype != 'uint8' else 255
                    ds = group.create_dataset(name, (2400, 2400), dtype=dtype, chunks=(100, 100), fillvalue=fill)
                    ds.attrs['_FillValue'] = fill
                    ds.attrs['scale_factor'] = 0.1 if key == 'lunar_irradiance' else 1
                    ds.attrs['add_offset'] = 0
                    ds[100:102, 100:102] = 48 if key == 'cloud_mask' else 0
                group[LAYERS['radiance']][100, 101] = 65535
                group[LAYERS['cloud_mask']][101, 100] = 48 | (3 << 6)
            result = extract(src, '2025-04-02', root / 'out', boundary, 'h09v05')
            self.assertEqual(result['county_cells'], 4)
            self.assertEqual(result['strict_usable_cells'], 2)
            with rasterio.open(root / 'out/radiance.tif') as ds:
                self.assertEqual(ds.read(1)[0, 0], 0)
                self.assertEqual(ds.nodata, 65535)
                self.assertEqual(ds.crs.to_epsg(), 4326)
                self.assertEqual(ds.tags(ns='IMAGE_STRUCTURE')['LAYOUT'], 'COG')
            with rasterio.open(root / 'out/cloud_mask.tif') as ds:
                self.assertEqual(ds.read(1)[1, 0], 240)
            with rasterio.open(root / 'out/lunar_irradiance.tif') as ds:
                self.assertAlmostEqual(ds.scales[0], 0.1)


if __name__ == '__main__':
    unittest.main()
