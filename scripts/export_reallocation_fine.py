#!/usr/bin/env python3
"""Retain nightly 10 m proxy scenes for exploration, verifying native aggregates."""
import argparse
from datetime import date
import json
import math
from pathlib import Path
import sys

import h5py
import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.ingest_viirs import GROUP, digest, now, save_json
from scripts.reallocation import upstream_modules, allocate_nightly, aggregate_complete
from scripts.run_counterfactual import load_cache


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream', type=Path, default=Path.home() / 'research/ntl-psf-disaggregation')
    parser.add_argument('--output', type=Path, default=ROOT / 'data/published/reallocation_check/fine')
    args = parser.parse_args()
    work = ROOT / 'data/reallocation'
    published = ROOT / 'data/published/reallocation_check'
    source_root = ROOT / 'data/viirs/VNP46A2.002'
    parent = json.loads((published / 'manifest.json').read_text())
    static = json.loads((work / 'static_manifest.json').read_text())
    contract = static['contract']
    if parent['status'] != 'complete' or parent['contract'] != contract:
        raise ValueError('Expected a complete matching native comparison')
    if digest(work / 'static_manifest.json') != parent['static_sha256']:
        raise ValueError('Static provenance changed')
    for name, sha in static['outputs'].items():
        if digest(work / name) != sha:
            raise ValueError(f'Static input changed: {name}')
    for name, sha in contract['upstream_sha256'].items():
        if digest(args.upstream / name) != sha:
            raise ValueError(f'Original implementation changed: {name}')
    for name in ('scripts/reallocation.py', 'scripts/run_counterfactual.py'):
        if digest(ROOT / name) != parent['code_sha256'][name]:
            raise ValueError(f'Comparison implementation changed: {name}')
    if digest(published / 'comparison_arrays.npz') != parent['output_sha256']['comparison_arrays.npz']:
        raise ValueError('Published arrays changed')

    args.output.mkdir(parents=True, exist_ok=True)
    receipt = args.output / 'manifest.json'
    identity = {'comparison_manifest_sha256': digest(published / 'manifest.json'),
                'export_script_sha256': digest(Path(__file__))}
    if receipt.exists():
        previous = json.loads(receipt.read_text())
        if (previous.get('status') == 'complete' and previous.get('identity') == identity
                and all((args.output / name).exists() and digest(args.output / name) == sha
                        for name, sha in previous['output_sha256'].items())):
            print('Verified existing fine-resolution export:', args.output)
            return
        raise ValueError('Existing export is incomplete or changed; use a separate --output directory')
    data = load_cache(source_root, date(2025, 4, 2))
    if data['hashes'] != parent['input_sha256']:
        raise ValueError('Nightly A2 inputs changed since the native comparison')
    with np.load(published / 'comparison_arrays.npz', allow_pickle=False) as saved:
        expected = saved['allocated'].copy()
        np.testing.assert_array_equal(saved['dates'], data['dates'])
        np.testing.assert_equal(saved['direct'], data['arrays']['observed'])

    modules = upstream_modules(args.upstream)
    kernel = modules['kernels'].circular_mean_kernel(radius_m=contract['kernel_radius_m'], resolution_m=10)
    corr = modules['operator']._correlate_constant_zero
    gain = np.load(work / 'gain.npy', mmap_mode='r')
    geom = np.load(work / 'geometric_support.npy', mmap_mode='r')
    mapping = np.load(work / 'fine_source_index.npy', mmap_mode='r')
    overlap = sparse.load_npz(work / 'native_overlap.npz')
    height, width = contract['shape']
    bbox = contract['bbox']
    row0, col0 = math.floor((40 - bbox[3]) * 240) - 1, math.floor((bbox[0] + 90) * 240) - 1
    row1, col1 = math.ceil((40 - bbox[1]) * 240) + 1, math.ceil((bbox[2] + 90) * 240) + 1
    sl = (slice(row0, row1), slice(col0, col1))
    cells = data['cells']
    stage = {'status': 'running', 'started_utc': now(), 'identity': identity,
             'dates': data['dates'], 'shape': [len(data['dates']), height, width],
             'axes': ['time', 'row', 'column'], 'dtype': 'float32',
             'units': 'nW cm-2 sr-1', 'resolution_m': 10,
             'crs': contract['crs'], 'bounds': contract['bounds'],
             'extent': 'Full original projected processing workspace, including halo outside Dyer.',
             'mask': 'native_aoi_mask.npy identifies fine pixels with positive area overlap with the selected native county cells; it is not applied to the cube.',
             'caveat': 'Structural reallocation proxy, not independent 10 m observations or validated fine-scale outages.',
             'raw_sha256': {}}
    save_json(receipt, stage)
    cube_path = args.output / 'reallocation_10m.npy'
    if cube_path.exists():
        raise FileExistsError(f'Refusing to overwrite an untracked cube: {cube_path}')
    cube = np.lib.format.open_memmap(cube_path, mode='w+', dtype='float32', shape=tuple(stage['shape']))
    for i, (day, src) in enumerate(zip(data['dates'], data['sources'])):
        if src['source_available']:
            path = source_root / 'raw' / src['source_filename']
            raw_sha = digest(path)
            cached = json.loads((work / 'nightly' / f'{day}.json').read_text())
            if raw_sha != cached['key']['source_sha256']:
                raise ValueError(f'Raw source differs from the native comparison: {day}')
            stage['raw_sha256'][src['source_filename']] = raw_sha
            with h5py.File(path, 'r') as source:
                g = source[GROUP]
                rad = g['DNB_BRDF-Corrected_NTL'][sl].astype('float32')
                fill = float(np.asarray(g['DNB_BRDF-Corrected_NTL'].attrs['_FillValue']).ravel()[0])
                cm, mq, snow = g['QF_Cloud_Mask'][sl], g['Mandatory_Quality_Flag'][sl], g['Snow_Flag'][sl]
                good = (np.isfinite(rad) & (rad != fill) & (mq == 0) & (cm != 65535)
                        & (((cm >> 6) & 3) == 0) & (((cm >> 4) & 3) == 3)
                        & ((cm & 1) == 0) & (((cm >> 9) & 1) == 0) & (snow == 0))
                rad = np.where(good, rad, np.nan)
            np.testing.assert_equal(rad[cells.source_row.to_numpy() - row0, cells.source_col.to_numpy() - col0],
                                    data['arrays']['observed'][i].astype('float32'))
            fine = rad.ravel()[mapping]
            allocated, uniform = allocate_nightly(fine, gain, kernel, corr, geom)
            del fine, uniform
        else:
            allocated = np.full((height, width), np.nan, dtype='float32')
        coarse, _ = aggregate_complete(overlap, allocated)
        np.testing.assert_equal(coarse, expected[i], err_msg=f'Native reaggregation mismatch on {day}')
        cube[i] = allocated
        del allocated
        print(f'{day}: saved 10 m scene; native aggregate matches exactly ({i + 1}/{len(data["dates"])})', flush=True)
    cube.flush()
    del cube
    # Positive polygon overlap preserves boundary pixels required for aggregation.
    mask = np.zeros(height * width, dtype=bool)
    mask[overlap.indices] = True
    np.save(args.output / 'native_aoi_mask.npy', mask.reshape(height, width))
    stage.update(status='complete', finished_utc=now(),
                 validation='All nightly native aggregates exactly match the published allocated array, including NaNs.',
                 output_sha256={name: digest(args.output / name)
                                for name in ('reallocation_10m.npy', 'native_aoi_mask.npy')})
    save_json(receipt, stage)
    print('Complete:', args.output, flush=True)


if __name__ == '__main__':
    main()
