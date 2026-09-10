#!/usr/bin/env python3
"""Chunked M0/M1 on retained 10 m allocation; no TING or outage classification."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.counterfactual import Settings, fit_models, forecast, evidence, validation_folds
from scripts.ingest_viirs import digest, now, save_json

FIELDS = ('center', 'q05', 'q95', 'residual', 'standardized_residual', 'lower_tail')
METRICS = ('mae', 'mse', 'bias', 'coverage90', 'below_q05', 'above_q95',
           'mean_interval_width', 'interval_score90', 'mean_negative_log_density')


def blocks(cube, mask, rows):
    for start in range(0, mask.shape[0], rows):
        stop = min(start + rows, mask.shape[0])
        selected = mask[start:stop].ravel()
        yield start, stop, selected, np.asarray(cube[:, start:stop]).reshape(len(cube), -1)[:, selected]


def training_pool(cube, mask, end, settings, rows):
    """Exact AOI-wide training quartiles/MAD pools, independent of chunk layout."""
    medians, deviations = [], []
    for _, _, _, values in blocks(cube, mask, rows):
        y = values[:end].astype(float)
        if np.isinf(y).any():
            raise ValueError('Infinite training input')
        supported = np.isfinite(y).sum(axis=0) >= settings.min_observations
        if not supported.any():
            continue
        y = y[:, supported]
        center = np.nanmedian(y, axis=0)
        mad = 1.4826 * np.nanmedian(np.abs(y - center), axis=0)
        medians.append(center)
        deviations.append(mad)
    if not medians:
        return np.zeros(3), np.full(4, settings.scale_floor)
    center, mad = np.concatenate(medians), np.concatenate(deviations)
    edges = np.quantile(center, [.25, .5, .75])
    groups = np.searchsorted(edges, center, side='right')
    global_scale = max(float(np.median(mad)), settings.scale_floor)
    scales = [max(float(np.median(mad[groups == group])), settings.scale_floor)
              if np.any(groups == group) else global_scale for group in range(4)]
    return edges, np.asarray(scales)


def score_sums(pred, actual, brightness):
    """Additive held-out scores; counts are descriptive, not independent trials."""
    use = np.isfinite(actual) & np.isfinite(pred['center'])
    row, col = np.where(use)
    actual = actual[use]
    center, scale, df, lo, hi = (pred[k][use] for k in ('center', 'scale', 'df', 'q05', 'q95'))
    error, width = actual - center, hi - lo
    components = (np.abs(error), error**2, error, (actual >= lo) & (actual <= hi),
                  actual < lo, actual > hi, width,
                  width + 20 * np.maximum(lo - actual, 0) + 20 * np.maximum(actual - hi, 0),
                  -student_t.logpdf(error / scale, df) + np.log(scale))
    masks = {'all': np.ones(len(col), bool), 'dim_lt1': brightness[col] < 1,
             'moderate_1to5': (brightness[col] >= 1) & (brightness[col] < 5),
             'bright_ge5': brightness[col] >= 5}
    return {name: {'n': int(m.sum()), **{key: float(v[m].sum()) for key, v in zip(METRICS, components)}}
            for name, m in masks.items()}


def write_cog(path, values, transform, crs, units, **tags):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.tmp.tif')
    with rasterio.open(temp, 'w', driver='COG', width=values.shape[1], height=values.shape[0],
                       count=1, dtype='float32', crs=crs, transform=transform, nodata=np.nan,
                       compress='DEFLATE', blocksize=256, overview_resampling='NEAREST',
                       NUM_THREADS='2') as dst:
        dst.write(values.astype('float32'), 1)
        dst.update_tags(units=units, description='10 m allocation counterfactual diagnostic; not an outage classification', **tags)
    temp.replace(path)


def run(source, output, rows):
    manifest_path = source / 'manifest.json'
    source_manifest = json.loads(manifest_path.read_text())
    if source_manifest['status'] != 'complete' or source_manifest['resolution_m'] != 10:
        raise ValueError('Expected a complete 10 m export')
    print('Verifying the fine-resolution cube and AOI mask...', flush=True)
    for name, sha in source_manifest['output_sha256'].items():
        if digest(source / name) != sha:
            raise ValueError(f'Fine input checksum mismatch: {name}')
    settings = Settings()
    identity = {'source_manifest_sha256': digest(manifest_path), 'settings': asdict(settings),
                'event_start': '2025-04-02', 'code_sha256': {name: digest(ROOT / name) for name in
                    ('scripts/run_counterfactual_fine.py', 'scripts/counterfactual.py')}}
    receipt = output / 'manifest.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        if (old.get('status') == 'complete' and old.get('identity') == identity
                and all((output / name).exists() and digest(output / name) == sha
                        for name, sha in old['output_sha256'].items())):
            print('Existing run verified:', output)
            return
        raise ValueError('Existing run is incomplete or changed; choose a separate --output directory')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output directory must be empty')
    output.mkdir(parents=True, exist_ok=True)
    cube = np.load(source / 'reallocation_10m.npy', mmap_mode='r')
    mask = np.load(source / 'native_aoi_mask.npy', mmap_mode='r')
    labels = source_manifest['dates']
    dates = pd.DatetimeIndex(labels)
    if not dates.equals(pd.date_range('2025-02-05', '2025-04-17')):
        raise ValueError('Expected the original complete calendar sequence')
    if tuple(source_manifest['shape']) != cube.shape or mask.shape != cube.shape[1:]:
        raise ValueError('Fine input grid mismatch')
    cut = labels.index(identity['event_start'])
    event_dates = labels[cut:]
    h, w = mask.shape
    bounds = source_manifest['bounds']
    transform = from_origin(bounds[0], bounds[3], 10, 10)
    stage = {'status': 'running', 'started_utc': now(), 'identity': identity,
             'baseline_dates': labels[:cut], 'event_dates': event_dates,
             'shape': list(cube.shape), 'axes': ['time', 'row', 'column'],
             'crs': source_manifest['crs'], 'transform': list(transform)[:6],
             'aoi_pixels': int(mask.sum()), 'rows_per_chunk': rows, 'ting_used': False,
             'support': 'Fine pixels overlapping selected native county cells; own nightly fine support, not the native complete-footprint mask.',
             'caveat': 'Correlated structural proxy pixels, not independent observations; uncalibrated intervals/tails are not outage probabilities.',
             'pools': {}}
    save_json(receipt, stage)
    baseline_grids = {key: np.full((h, w), np.nan, dtype='float32')
                      for key in ('baseline_usable_nights', 'baseline_supported', 'baseline_median')}
    model_arrays = {}
    for method in ('M0', 'M1'):
        folder = output / method
        folder.mkdir()
        model_arrays[method] = {key: np.lib.format.open_memmap(folder / f'{key}.npy', mode='w+',
                                dtype='float32', shape=(len(event_dates), h, w)) for key in FIELDS}
        # Explicit NaN outside the analysis footprint, including blocks with no AOI cells.
        for value in model_arrays[method].values():
            value[:] = np.nan
            value.flush()
    folds = validation_folds(cut)
    summaries, daily = [], {}
    for end in sorted({f['train_days'] for f in folds} | {cut}):
        print(f'Training through {labels[end - 1]}: computing full-AOI brightness pools...', flush=True)
        pool = training_pool(cube, mask, end, settings, rows)
        stage['pools'][str(end)] = {'edges': pool[0].tolist(), 'scales': pool[1].tolist()}
        save_json(receipt, stage)
        relevant = [f for f in folds if f['train_days'] == end]
        width = len(event_dates) if end == cut else max(f['test_days'] for f in relevant)
        accum, eligible = {}, {}
        supported_total = 0
        for start, stop, selected, values in blocks(cube, mask, rows):
            if not selected.any():
                continue
            fitted = fit_models(values[:end], settings, scale_pool=pool)
            supported_total += int(fitted['supported'].sum())
            strata = {'all': np.ones(len(fitted['count']), bool),
                      'dim_lt1': fitted['baseline_median'] < 1,
                      'moderate_1to5': (fitted['baseline_median'] >= 1) & (fitted['baseline_median'] < 5),
                      'bright_ge5': fitted['baseline_median'] >= 5}
            for name, m in strata.items():
                eligible[name] = eligible.get(name, 0) + int((m & fitted['supported']).sum())
            if end == cut:
                for key, values_1d in (('baseline_usable_nights', fitted['count']),
                                       ('baseline_supported', fitted['supported']),
                                       ('baseline_median', fitted['baseline_median'])):
                    baseline_grids[key][start:stop].reshape(-1)[selected] = values_1d
            for method in ('M0', 'M1'):
                pred = forecast(fitted, method, np.arange(1, width + 1))
                actual = values[end:end + width]
                if end == cut:
                    scores = evidence(pred, actual)
                    combined = {**pred, **scores}
                    for field in FIELDS:
                        block = model_arrays[method][field][:, start:stop].reshape(width, -1)
                        block[:, selected] = combined[field].astype('float32')
                    for j, day in enumerate(event_dates):
                        key = (method, day)
                        record = daily.setdefault(key, {'method': method, 'date': day,
                            'usable_pixels': 0, 'evaluable_pixels': 0, 'sum_residual': 0., 'below_q05_pixels': 0})
                        use = np.isfinite(scores['residual'][j])
                        record['usable_pixels'] += int(np.isfinite(actual[j]).sum())
                        record['evaluable_pixels'] += int(use.sum())
                        record['sum_residual'] += float(scores['residual'][j, use].sum())
                        record['below_q05_pixels'] += int((scores['lower_tail'][j, use] < .05).sum())
                else:
                    for fold in relevant:
                        n = fold['test_days']
                        scored = score_sums({k: v[:n] for k, v in pred.items()}, actual[:n], fitted['baseline_median'])
                        for stratum, sums in scored.items():
                            key = (method, fold['design'], stratum)
                            total = accum.setdefault(key, {k: 0 for k in sums})
                            for name, value in sums.items():
                                total[name] += value
                del pred
            if start % (rows * 8) == 0:
                print(f'  train={end} rows {start}:{stop}/{h}', flush=True)
        if end != cut:
            for (method, design, stratum), total in accum.items():
                n = total.pop('n')
                means = {k: value / n if n else None for k, value in total.items()}
                mse = means.pop('mse')
                means['rmse'] = float(np.sqrt(mse)) if mse is not None else None
                days = int(design[:-1])
                summaries.append({'method': method, 'design': design, 'stratum': stratum,
                    'train_days': end, 'test_days': days, 'train_end': labels[end - 1],
                    'test_start': labels[end], 'test_end': labels[end + days - 1],
                    'eligible_cells': eligible[stratum], 'n': n, **means})
            pd.DataFrame(summaries).to_csv(output / 'validation_summary.csv', index=False)
        else:
            stage['supported_pixels'] = supported_total
        print(f'Completed training origin {end}: {supported_total:,} supported fine pixels.', flush=True)
    for method in model_arrays.values():
        for value in method.values():
            value.flush()
    for key, values in baseline_grids.items():
        write_cog(output / f'{key}.tif', values, transform, source_manifest['crs'],
                  'nights' if key.endswith('nights') else ('0/1' if key.endswith('supported') else 'nW cm-2 sr-1'))
    for method in ('M0', 'M1'):
        for j, day in enumerate(event_dates):
            for field in FIELDS:
                write_cog(output / method / 'daily' / day / f'{field}.tif', model_arrays[method][field][j],
                          transform, source_manifest['crs'],
                          'dimensionless' if field in ('lower_tail', 'standardized_residual') else 'nW cm-2 sr-1',
                          method=method, date=day, event_start=event_dates[0])
            print(f'Exported {method} {day} COGs.', flush=True)
    for record in daily.values():
        n = record['evaluable_pixels']
        record['supported_pixels'] = stage['supported_pixels']
        record['mean_residual'] = record.pop('sum_residual') / n if n else None
        record['fraction_below_q05'] = record['below_q05_pixels'] / n if n else None
    pd.DataFrame(daily.values()).to_csv(output / 'event_daily_summary.csv', index=False)
    (output / 'README.txt').write_text(
        'M0/M1 on nightly 10 m structural allocation; NOT validated 10 m outage maps.\n'
        'Same baseline, event period, minimum observations, scale floor, prior weight, and M1 process ratio as native runs.\n'
        'AOI: fine pixels overlapping selected native Dyer cells; no processing halo in fits or pooling.\n'
        'Own nightly fine support, not the more restrictive native aggregation support.\n'
        'M0 and M1 share observation masks; each validation origin has training-only full-AOI brightness pools.\n'
        'M0/M1/*.npy: memory-mapped float32 (16 event dates, row, column) fields, matching full fine workspace.\n'
        'M0/M1/daily/YYYY-MM-DD/*.tif: same fields as EPSG:32616 10 m COGs for QGIS.\n'
        'center/q05/q95: forecast and nominal central 90% reference interval.\n'
        'residual: observed minus predicted; standardized_residual: residual / predictive SD.\n'
        'lower_tail: uncalibrated predictive CDF, not outage probability. Missing observations yield no evidence.\n'
        'Fine pixels share native radiance and structural assumptions, so validation counts are not independent samples.\n'
        'Contemporary 2026 structural inputs applied to 2025 radiance; no TING, controls, or causal attribution.\n'
        'Analyze 7d and 14d validation designs separately; their target dates overlap across designs.\n')
    print('Hashing completed outputs...', flush=True)
    stage.update(status='complete', finished_utc=now(), output_sha256={str(p.relative_to(output)): digest(p)
                 for p in sorted(output.rglob('*')) if p.is_file() and p != receipt})
    save_json(receipt, stage)
    print('Complete:', output, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'data/published/reallocation_check/fine')
    parser.add_argument('--output', type=Path, default=ROOT / 'data/published/counterfactual_reallocation_10m')
    parser.add_argument('--chunk-rows', type=int, default=32)
    args = parser.parse_args()
    if args.chunk_rows < 1:
        parser.error('--chunk-rows must be positive')
    run(args.input, args.output, args.chunk_rows)


if __name__ == '__main__':
    main()
