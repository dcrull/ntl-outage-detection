#!/usr/bin/env python3
"""Export daily allocated-minus-direct A2 radiance on common native support."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from rasterio.transform import Affine

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.ingest_viirs import cog, digest, now, save_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'data/published/reallocation_check')
    args = parser.parse_args()
    source = args.input / 'comparison_arrays.npz'
    source_manifest = args.input / 'manifest.json'
    manifest = json.loads(source_manifest.read_text())
    if manifest['status'] != 'complete' or digest(source) != manifest['output_sha256'][source.name]:
        raise ValueError('A complete, verified reallocation comparison is required')
    output = args.input / 'differences'
    output.mkdir(parents=True, exist_ok=True)
    provenance = {
        'status': 'running', 'created_utc': now(),
        'source_manifest_sha256': digest(source_manifest),
        'source_arrays_sha256': digest(source), 'script_sha256': digest(Path(__file__)),
        'formula': 'allocated minus direct A2',
        'units': 'nW cm-2 sr-1',
        'support': 'Shared valid cell-night mask from the reallocation comparison',
    }
    save_json(output / 'manifest.json', provenance)
    hashes = {}
    with np.load(source, allow_pickle=False) as arrays:
        county = arrays['county']
        transform = Affine(*arrays['transform'])
        differences = np.where(arrays['common'], arrays['allocated'] - arrays['direct'], np.nan)
        for day, values in zip(arrays['dates'], differences):
            grid = np.full(county.shape, np.nan, dtype='float32')
            grid[county] = values
            path = output / f'{day}.tif'
            cog(path, grid, transform, np.nan, units='nW cm-2 sr-1',
                date=str(day), formula='allocated minus direct A2',
                description='Native-grid input difference; not an outage classification')
            hashes[path.name] = digest(path)
    (output / 'README.txt').write_text(
        'Daily reallocated-minus-direct A2 radiance, on the original VIIRS grid.\n'
        'Positive: allocation is brighter. Negative: allocation is dimmer.\n'
        'Zero: equal radiance. NoData: outside county or missing common support.\n'
        'Units: nW cm^-2 sr^-1. Acquisition-day labels, not exact local overpass times.\n'
        'QGIS: Singleband pseudocolor, diverging blue-white-red, zero at white; use symmetric limits.\n'
        'These compare radiance inputs, not M0/M1 residuals or outage states.\n')
    hashes['README.txt'] = digest(output / 'README.txt')
    provenance.update(status='complete', output_sha256=hashes)
    save_json(output / 'manifest.json', provenance)
    print(f'Exported {len(differences)} daily difference COGs to {output}')


if __name__ == '__main__':
    main()
