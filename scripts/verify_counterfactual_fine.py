#!/usr/bin/env python3
"""Audit completed 10 m model outputs against inputs, model fits, and COG grids."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import rasterio
from rasterio.transform import Affine

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.counterfactual import Settings, fit_models, forecast, evidence
from scripts.ingest_viirs import digest, now, save_json
from scripts.run_counterfactual_fine import FIELDS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'data/published/counterfactual_reallocation_10m')
    args = parser.parse_args()
    output = args.output
    manifest = json.loads((output / 'manifest.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Run is not complete')
    print('Verifying output hashes...', flush=True)
    for name, sha in manifest['output_sha256'].items():
        if digest(output / name) != sha:
            raise ValueError(f'Output checksum mismatch: {name}')
    source = ROOT / 'data/published/reallocation_check/fine'
    if digest(source / 'manifest.json') != manifest['identity']['source_manifest_sha256']:
        raise ValueError('Input provenance changed')
    mask = np.load(source / 'native_aoi_mask.npy', mmap_mode='r')
    cube = np.load(source / 'reallocation_10m.npy', mmap_mode='r')
    cut = len(manifest['baseline_dates'])
    event_dates = manifest['event_dates']
    transform = Affine(*manifest['transform'])
    # Deterministic distributed cells, independent of the fitter's row chunks.
    aoi_indices = np.flatnonzero(mask)
    sample = aoi_indices[np.linspace(0, len(aoi_indices) - 1, 2048, dtype=int)]
    del aoi_indices
    selected = cube.reshape(len(cube), -1)[:, sample]
    pool = manifest['pools'][str(cut)]
    fit = fit_models(selected[:cut], Settings(**manifest['identity']['settings']),
                     scale_pool=(pool['edges'], pool['scales']))
    count_path = output / 'baseline_usable_nights.tif'
    with rasterio.open(count_path) as src:
        counts = src.read(1)
        if src.crs.to_string() != manifest['crs'] or src.transform != transform:
            raise ValueError('Baseline grid mismatch')
    supported = counts >= manifest['identity']['settings']['min_observations']
    assert int(supported.sum()) == manifest['supported_pixels']
    rasters_checked = 0
    for method in ('M0', 'M1'):
        prediction = forecast(fit, method, np.arange(1, len(event_dates) + 1))
        expected = {**prediction, **evidence(prediction, selected[cut:])}
        for field in FIELDS:
            print(f'Checking {method} {field}...', flush=True)
            array = np.load(output / method / f'{field}.npy', mmap_mode='r')
            assert array.shape == (len(event_dates), *mask.shape)
            np.testing.assert_equal(array.reshape(len(event_dates), -1)[:, sample], expected[field].astype('float32'))
            for start in range(0, mask.shape[0], 61):
                stop = min(start + 61, mask.shape[0])
                block = array[:, start:stop]
                expected_valid = np.broadcast_to(supported[start:stop], block.shape)
                if field in ('residual', 'standardized_residual', 'lower_tail'):
                    expected_valid = expected_valid & np.isfinite(cube[cut:, start:stop])
                np.testing.assert_array_equal(np.isfinite(block), expected_valid)
                assert not np.isinf(block).any()
                if field == 'lower_tail':
                    valid_values = block[expected_valid]
                    assert np.all((valid_values >= 0) & (valid_values <= 1))
            for j, day in enumerate(event_dates):
                path = output / method / 'daily' / day / f'{field}.tif'
                with rasterio.open(path) as src:
                    assert src.shape == mask.shape and src.transform == transform
                    assert src.crs.to_string() == manifest['crs']
                    assert src.tags(ns='IMAGE_STRUCTURE').get('LAYOUT') == 'COG'
                    assert np.isnan(src.nodata)
                    # Check full windows distributed across the image, including its edges.
                    for start in (0, mask.shape[0] // 2, mask.shape[0] - 3):
                        window = rasterio.windows.Window(0, start, mask.shape[1], 3)
                        np.testing.assert_equal(src.read(1, window=window), array[j, start:start + 3])
                rasters_checked += 1
    report = {'status': 'passed', 'verified_utc': now(),
              'manifest_sha256': digest(output / 'manifest.json'),
              'output_hashes_verified': len(manifest['output_sha256']),
              'cogs_checked': rasters_checked, 'sampled_model_pixels': len(sample),
              'full_cube_missingness_checked': True, 'supported_pixels': int(supported.sum())}
    save_json(ROOT / 'data/reallocation/counterfactual_10m_validation.json', report)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
