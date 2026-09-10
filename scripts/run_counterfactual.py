#!/usr/bin/env python3
"""Fit and evaluate M0/M1 from the local non-gap-filled VNP46A2 cache.

Run: .venv/bin/python scripts/run_counterfactual.py
No network requests, TING inputs, gap-filled radiance, or outage classifications.
"""
import argparse
from dataclasses import asdict
from datetime import date
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin, xy
import scipy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.counterfactual import Settings, evidence, fit_models, forecast, metrics, validation_folds
from scripts.ingest_viirs import cog, dates_between, digest, now, save_json

INPUT_LAYERS = ('radiance', 'strict_usable', 'cloud_mask', 'mandatory_quality', 'snow_flag', 'county_mask')


def load_cache(root, event_start):
    manifest = json.loads((root / 'manifest.json').read_text())
    if manifest.get('status') not in ('complete', 'complete_with_missing_dates'):
        raise ValueError('Ingestion must be complete before analysis')
    if (manifest['product'], manifest['collection'], manifest['tile']) != ('VNP46A2', '002', 'h09v05'):
        raise ValueError('Expected the Dyer VNP46A2.002 h09v05 cache')
    days = dates_between(date.fromisoformat(manifest['start']), date.fromisoformat(manifest['end']))
    labels = [d.isoformat() for d in days]
    cutoff = sum(d < event_start for d in days)
    if cutoff < 42 or cutoff == len(days):
        raise ValueError('Need at least 42 pre-event calendar days and an event period')
    coverage = pd.read_csv(root / 'coverage.csv').set_index('date')
    if not coverage.index.is_unique:
        raise ValueError('Duplicate coverage dates')
    if set(coverage.index) != set(labels):
        raise ValueError('Coverage dates must match the active ingestion manifest')
    first_day = next(d for d in labels if coverage.loc[d, 'status'] == 'complete')
    first = root / 'dyer' / first_day
    meta = json.loads((first / 'metadata.json').read_text())
    with rasterio.open(first / 'county_mask.tif') as src:
        county, transform, crs = src.read(1) == 1, src.transform, src.crs
    expected = from_origin(-90 + meta['column_start'] / 240, 40 - meta['row_start'] / 240, 1 / 240, 1 / 240)
    if crs.to_epsg() != 4326 or not transform.almost_equals(expected):
        raise ValueError('Crop is not aligned to the native h09v05 grid')
    if meta['tile'] != 'h09v05' or not county.any():
        raise ValueError('Invalid tile metadata or county mask')
    rows, cols = np.where(county)
    shape = (len(days), len(rows))
    arrays = {name: np.full(shape, np.nan if name == 'radiance' else 65535,
                            dtype=float if name == 'radiance' else 'uint16')
              for name in ('radiance', 'cloud_mask', 'mandatory_quality', 'snow_flag')}
    arrays['strict_usable'] = np.zeros(shape, dtype=bool)
    records, hashes = [], {}
    for name in ('manifest.json', 'inventory.json', 'coverage.csv'):
        hashes[name] = digest(root / name)
    for i, day in enumerate(labels):
        status = coverage.loc[day, 'status']
        if status != 'complete':
            if status != 'missing_source':
                raise ValueError(f'Unexpected input status {status} on {day}')
            records.append({'date': day, 'source_available': False, 'source_filename': ''})
            continue
        directory = root / 'dyer' / day
        metadata = json.loads((directory / 'metadata.json').read_text())
        if any(metadata[k] != meta[k] for k in ('row_start', 'column_start', 'tile', 'crs', 'transform')):
            raise ValueError(f'Metadata grid changed on {day}')
        if metadata['date'] != day or metadata['source_filename'] != coverage.loc[day, 'source_filename']:
            raise ValueError(f'Source date/provenance mismatch on {day}')
        hashes[str((directory / 'metadata.json').relative_to(root))] = digest(directory / 'metadata.json')
        for name in INPUT_LAYERS:
            path = directory / f'{name}.tif'
            sha = digest(path)
            if sha != metadata['sha256'][path.name]:
                raise ValueError(f'Input hash mismatch: {path}')
            hashes[str(path.relative_to(root))] = sha
            with rasterio.open(path) as src:
                if src.shape != county.shape or src.transform != transform or src.crs != crs:
                    raise ValueError(f'Raster grid changed: {path}')
                if name == 'radiance':
                    if src.tags().get('source_layer') != 'DNB_BRDF-Corrected_NTL':
                        raise ValueError('Expected non-gap-filled corrected radiance source tag')
                    values = src.read(1, masked=True).astype(float).filled(np.nan)
                    values = values * src.scales[0] + src.offsets[0]
                else:
                    values = src.read(1)
                    if name in ('county_mask', 'strict_usable') and not np.isin(values, [0, 1]).all():
                        raise ValueError(f'Expected binary mask: {path}')
                if name == 'county_mask':
                    if not np.array_equal(values == 1, county):
                        raise ValueError('County mask changed between dates')
                else:
                    arrays[name][i] = values[county]
        cm = arrays['cloud_mask'][i]
        strict = (np.isfinite(arrays['radiance'][i]) & (arrays['mandatory_quality'][i] == 0)
                  & (cm != 65535) & (((cm >> 6) & 3) == 0) & (((cm >> 4) & 3) == 3)
                  & ((cm & 1) == 0) & (((cm >> 9) & 1) == 0) & (arrays['snow_flag'][i] == 0))
        if not np.array_equal(strict, arrays['strict_usable'][i]):
            raise ValueError(f'Cached screening disagrees with independent QA check: {day}')
        records.append({'date': day, 'source_available': True, 'source_filename': metadata['source_filename']})
    arrays['observed'] = np.where(arrays['strict_usable'], arrays['radiance'], np.nan)
    source_rows, source_cols = rows + meta['row_start'], cols + meta['column_start']
    lon, lat = xy(transform, rows, cols)
    cells = pd.DataFrame({'cell_id': [f'h09v05_r{r:04d}_c{c:04d}' for r, c in zip(source_rows, source_cols)],
                          'source_row': source_rows, 'source_col': source_cols,
                          'crop_row': rows, 'crop_col': cols, 'longitude': lon, 'latitude': lat})
    return {'arrays': arrays, 'dates': labels, 'cutoff': cutoff, 'county': county, 'transform': transform,
            'cells': cells, 'sources': records, 'hashes': hashes, 'manifest': manifest}


def evaluate(observed, cutoff, settings, labels):
    summary, daily = [], []
    archive = {}
    for fold in validation_folds(cutoff):
        end, width = fold['train_days'], fold['test_days']
        fold_id = f"{fold['design']}_origin_{labels[end-1]}"
        fit = fit_models(observed[:end], settings)
        actual = observed[end:end + width]
        strata = {'all': np.ones(observed.shape[1], bool),
                  'dim_lt1': fit['baseline_median'] < 1,
                  'moderate_1to5': (fit['baseline_median'] >= 1) & (fit['baseline_median'] < 5),
                  'bright_ge5': fit['baseline_median'] >= 5}
        for method in ('M0', 'M1'):
            pred = forecast(fit, method, np.arange(1, width + 1))
            prefix = f'{fold_id}_{method}'
            for name, array in pred.items():
                archive[f'{prefix}_{name}'] = array.astype('float32')
            for stratum, mask in strata.items():
                common = {**fold, 'fold': fold_id, 'method': method, 'stratum': stratum,
                          'train_end': labels[end-1], 'test_start': labels[end], 'test_end': labels[end+width-1],
                          'eligible_cells': int((fit['supported'] & mask).sum())}
                summary.append({**common, **metrics(pred, actual, mask)})
                for j in range(width):
                    daily.append({**common, 'date': labels[end+j], 'horizon_days': j+1,
                                  **metrics({k: v[j] for k, v in pred.items()}, actual[j], mask)})
    return pd.DataFrame(summary), pd.DataFrame(daily), archive


def run(input_root, output, event_start, settings):
    print('Reading and verifying non-gap-filled A2 inputs...', flush=True)
    data = load_cache(input_root, event_start)
    arrays, labels, cutoff = data['arrays'], data['dates'], data['cutoff']
    config = {'input_root': str(input_root.resolve()), 'event_start': event_start.isoformat(),
              'settings': asdict(settings), 'source_manifest_sha256': data['hashes']['manifest.json']}
    manifest_path = output / 'manifest.json'
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text())
        if old.get('config') != config:
            raise ValueError('Existing output uses a different configuration. Choose a separate --output directory.')
    elif output.exists() and any(output.iterdir()):
        raise ValueError('Output directory is nonempty and has no analysis manifest')
    output.mkdir(parents=True, exist_ok=True)
    manifest = {'status': 'running', 'started_utc': now(), 'config': config,
                'input_sha256': data['hashes'], 'code_sha256': {
                    str(p.relative_to(ROOT)): digest(p) for p in
                    [Path(__file__), ROOT / 'scripts/counterfactual.py', ROOT / 'scripts/ingest_viirs.py']},
                'versions': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__,
                             'scipy': scipy.__version__, 'rasterio': rasterio.__version__},
                'radiance_layer': 'DNB_BRDF-Corrected_NTL', 'gap_filled_used': False, 'ting_used': False,
                'timing': 'Product acquisition dates; exact local overpass times are not yet resolved.'}
    save_json(manifest_path, manifest)
    print('Evaluating frozen 7- and 14-day pre-event forecasts...', flush=True)
    summary, daily_validation, validation_arrays = evaluate(arrays['observed'], cutoff, settings, labels)
    outputs = []

    def csv(name, frame):
        path = output / name
        frame.to_csv(path, index=False, float_format='%.8g')
        outputs.append(path)

    csv('validation_summary.csv', summary)
    csv('validation_daily.csv', daily_validation)
    np.savez_compressed(output / 'validation_predictions.npz', **validation_arrays)
    outputs.append(output / 'validation_predictions.npz')
    fit = fit_models(arrays['observed'][:cutoff], settings)
    cells = data['cells'].copy()
    for name in ('count', 'supported', 'baseline_median', 'baseline_mad_sigma', 'pool_group', 'pooled_sigma'):
        cells[name] = fit[name]
    cells['M0_reference'] = np.where(fit['supported'], fit['M0']['center'], np.nan)
    cells['M1_reference'] = np.where(fit['supported'], fit['M1']['center'], np.nan)
    csv('cells.csv', cells)
    csv('source_dates.csv', pd.DataFrame(data['sources']))
    saved = {**arrays, 'dates': np.array(labels), 'county': data['county'],
             'crop_rows': cells['crop_row'].to_numpy(), 'crop_cols': cells['crop_col'].to_numpy(),
             'transform': np.array(list(data['transform'])[:6]), 'cutoff': np.array(cutoff)}

    def raster(name, values, units='nW cm-2 sr-1'):
        grid = np.full(data['county'].shape, np.nan, dtype='float32')
        grid[data['county']] = values
        path = output / name
        cog(path, grid, data['transform'], np.nan, units=units,
            description='Counterfactual radiance diagnostic; not an outage classification',
            source_layer='DNB_BRDF-Corrected_NTL', event_start=event_start.isoformat())
        outputs.append(path)

    raster('baseline_usable_nights.tif', fit['count'], 'nights')
    raster('baseline_supported.tif', fit['supported'].astype(float), '0/1')
    raster('baseline_median.tif', fit['baseline_median'])
    event_records, daily_records = [], []
    actual = arrays['observed'][cutoff:]
    horizons = np.arange(1, len(actual) + 1)
    print(f"Fitting {int(fit['supported'].sum())}/{len(cells)} cells and exporting {len(actual)} nights...", flush=True)
    for method in ('M0', 'M1'):
        prediction = forecast(fit, method, horizons)
        scores = evidence(prediction, actual)
        for name, array in {**prediction, **scores}.items():
            saved[f'{method}_{name}'] = array
        for j, day in enumerate(labels[cutoff:]):
            source_idx = cutoff + j
            supported = fit['supported']
            observed = np.isfinite(actual[j])
            evaluable = supported & observed
            status = np.where(~supported, 'insufficient_baseline',
                              np.where(observed, 'observed', 'no_usable_observation'))
            record = pd.DataFrame({'cell_id': cells['cell_id'], 'method': method, 'date': day,
                                   'forecast_horizon_days': j+1, 'training_nights': fit['count'],
                                   'supported': supported, 'usable': observed, 'status': status,
                                   'source_available': data['sources'][source_idx]['source_available'],
                                   'source_filename': data['sources'][source_idx]['source_filename'],
                                   'raw_corrected_radiance': arrays['radiance'][source_idx],
                                   'observed_radiance': actual[j],
                                   'cloud_mask': arrays['cloud_mask'][source_idx],
                                   'mandatory_quality': arrays['mandatory_quality'][source_idx],
                                   'snow_flag': arrays['snow_flag'][source_idx],
                                   **{name: v[j] for name, v in {**prediction, **scores}.items()}})
            event_records.append(record)
            daily_records.append({'date': day, 'method': method, 'county_cells': len(cells),
                                  'supported_cells': int(supported.sum()), 'usable_cells': int(observed.sum()),
                                  'evaluable_cells': int(evaluable.sum()),
                                  'mean_residual': float(np.mean(scores['residual'][j, evaluable])) if evaluable.any() else None,
                                  'fraction_below_q05': float(np.mean(scores['lower_tail'][j, evaluable] < .05)) if evaluable.any() else None})
            base = f'{method}/daily/{day}'
            for key in ('center', 'q05', 'q95'):
                raster(f'{base}/{key}.tif', prediction[key][j])
            for key in ('residual', 'standardized_residual', 'lower_tail'):
                raster(f'{base}/{key}.tif', scores[key][j],
                       'nW cm-2 sr-1' if key == 'residual' else 'dimensionless')
            if method == 'M0':
                raster(f'observed/{day}.tif', actual[j])
    csv('event_cell_nights.csv.gz', pd.concat(event_records, ignore_index=True))
    csv('event_daily_summary.csv', pd.DataFrame(daily_records))
    np.savez_compressed(output / 'analysis_arrays.npz', **saved)
    outputs.append(output / 'analysis_arrays.npz')
    model_info = {'settings': asdict(settings), 'baseline_start': labels[0], 'baseline_end': labels[cutoff-1],
                  'forecast_start': labels[cutoff], 'forecast_end': labels[-1],
                  'supported_cells': int(fit['supported'].sum()), 'county_cells': len(cells),
                  'brightness_pool_edges': fit['pool_edges'].tolist(), 'pooled_sigma': fit['pool_scales'].tolist(),
                  'validation_folds': validation_folds(cutoff),
                  'intervals': 'Central 90% uncalibrated Student-t predictive references; may extend below zero.',
                  'M0': 'Median + MAD variance shrunk toward training brightness-group scale; approximate median uncertainty.',
                  'M1': 'Fixed q/theta=.001 by default; no trend; inverse-gamma scale prior; diffuse initial location. Empirical-Bayes hyperparameter uncertainty omitted.',
                  'negative_controls': 'None; no causal attribution or outage state inferred.'}
    save_json(output / 'model_details.json', model_info)
    outputs.append(output / 'model_details.json')
    (output / 'README.txt').write_text(
        'M0/M1 non-gap-filled A2 counterfactual radiance diagnostics; NOT outage maps.\n'
        'Native h09v05 grid; county cell centers. Product dates are not local overpass timestamps.\n'
        'center/q05/q95: model reference and nominal central 90% predictive interval, not calibrated outage confidence.\n'
        'residual: observed minus predicted radiance; negative = dimmer. Blank = insufficient evidence.\n'
        'lower_tail: uncalibrated reference-model CDF at observed radiance; not outage probability.\n'
        'Cloudy nights may have a forecast, but no residual/tail evidence. No gap-filled data or TING used.\n'
        'See docs/COUNTERFACTUAL_IMPLEMENTATION.md and notebooks/03_inspect_counterfactual.ipynb in the repo.\n'
        'The 7d and 14d validation designs reuse dates across designs; analyze them separately.\n')
    outputs.append(output / 'README.txt')
    manifest.update(status='complete', finished_utc=now(),
                    output_sha256={str(p.relative_to(output)): digest(p) for p in outputs})
    save_json(manifest_path, manifest)
    print(f'Completed. Outputs: {output}', flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'data/viirs/VNP46A2.002')
    parser.add_argument('--output', type=Path, default=ROOT / 'data/published/counterfactual_dyer')
    parser.add_argument('--event-start', type=date.fromisoformat, default=date(2025, 4, 2))
    parser.add_argument('--min-observations', type=int, default=10)
    parser.add_argument('--scale-floor', type=float, default=.1)
    parser.add_argument('--scale-prior-weight', type=float, default=20)
    parser.add_argument('--process-ratio', type=float, default=.001)
    args = parser.parse_args()
    settings = Settings(args.min_observations, args.scale_floor, args.scale_prior_weight, args.process_ratio)
    run(args.input, args.output, args.event_start, settings)


if __name__ == '__main__':
    main()
