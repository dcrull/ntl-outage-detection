import unittest
from datetime import date, timedelta
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from rasterio.transform import from_origin

from scripts.counterfactual import Settings, evidence, fit_models, forecast, metrics, validation_folds
from scripts.run_counterfactual import evaluate, load_cache
from scripts.ingest_viirs import cog, digest, save_json


class CounterfactualTests(unittest.TestCase):
    def test_m0_median_resists_single_large_training_value(self):
        values = np.array([1, 2, 2, 3, 4, 100.])[:, None]
        fitted = fit_models(values, Settings(min_observations=3))
        predicted = forecast(fitted, 'M0', [1, 7])
        np.testing.assert_array_equal(predicted['center'], [[2.5], [2.5]])
        self.assertTrue(np.all(predicted['q05'] < predicted['center']))
        self.assertTrue(np.all(predicted['q95'] > predicted['center']))

    def test_zero_is_observed_but_sparse_and_missing_cells_are_unsupported(self):
        y = np.column_stack([np.zeros(12), [1] + [np.nan] * 11, np.full(12, np.nan)])
        fitted = fit_models(y)
        np.testing.assert_array_equal(fitted['supported'], [True, False, False])
        for model in ('M0', 'M1'):
            p = forecast(fitted, model, [1])
            self.assertEqual(p['center'][0, 0], 0)
            self.assertGreater(p['scale'][0, 0], 0)
            self.assertTrue(np.isnan(p['center'][0, 1:]).all())

    def test_static_m1_matches_analytic_unknown_variance_posterior(self):
        y = np.array([1., 3., 4., 2., 5.])[:, None]
        settings = Settings(min_observations=3, process_ratio=0)
        fitted = fit_models(y, settings)
        p = fitted['M1']
        self.assertAlmostEqual(p['center'][0], np.mean(y))
        self.assertAlmostEqual(p['state_factor'][0], 1 / len(y))
        self.assertAlmostEqual(p['df'][0], settings.scale_prior_weight + len(y) - 1)
        expected_sum = settings.scale_prior_weight * fitted['pooled_sigma'][0]**2 + np.sum((y - y.mean())**2)
        self.assertAlmostEqual(p['variance'][0], expected_sum / p['df'][0])

    def test_missing_calendar_days_propagate_state_uncertainty(self):
        y = np.arange(12.)[:, None]
        settings = Settings(process_ratio=.001)
        a = fit_models(y, settings)
        b = fit_models(np.concatenate([y, np.full((5, 1), np.nan)]), settings)
        np.testing.assert_allclose(a['M1']['center'], b['M1']['center'])
        np.testing.assert_allclose(b['M1']['state_factor'] - a['M1']['state_factor'], .005)
        # Five unobserved training dates followed by day 1 equals a six-day forecast.
        pa, pb = forecast(a, 'M1', [6]), forecast(b, 'M1', [1])
        np.testing.assert_allclose(pa['scale'], pb['scale'])
        p = forecast(a, 'M1', [1, 16])
        self.assertGreater(p['scale'][1, 0], p['scale'][0, 0])
        self.assertEqual(p['center'][1, 0], p['center'][0, 0])

    def test_missing_observations_never_produce_dimming_evidence(self):
        fitted = fit_models(np.full((12, 2), 5.))
        for method in ('M0', 'M1'):
            prediction = forecast(fitted, method, [1, 2])
            result = evidence(prediction, np.array([[np.nan, np.nan], [0., 5.]]))
            for values in result.values():
                self.assertTrue(np.isnan(values[0]).all())
            self.assertLess(result['lower_tail'][1, 0], .05)
            self.assertAlmostEqual(result['lower_tail'][1, 1], .5)
            self.assertEqual(metrics(prediction, np.full((2, 2), np.nan))['n'], 0)

    def test_validation_forecasts_do_not_use_their_holdouts_or_event(self):
        rng = np.random.default_rng(7)
        y = rng.normal(5, 1, (72, 4))
        labels = [str(d) for d in range(72)]
        _, _, original = evaluate(y, 56, Settings(), labels)
        changed = y.copy()
        changed[28:] = 10000
        _, _, alternative = evaluate(changed, 56, Settings(), labels)
        for key in original:
            if 'origin_27_' in key:
                np.testing.assert_equal(original[key], alternative[key])
        changed = y.copy()
        changed[56:] = 0
        _, _, event_changed = evaluate(changed, 56, Settings(), labels)
        for key in original:
            np.testing.assert_equal(original[key], event_changed[key])

    def test_fold_targets_disjoint_within_design_and_pre_event(self):
        folds = validation_folds(56)
        self.assertEqual(len(folds), 6)
        for design in ('7d', '14d'):
            targets = []
            for f in folds:
                if f['design'] == design:
                    targets.extend(range(f['train_days'], f['train_days'] + f['test_days']))
            self.assertEqual(len(targets), len(set(targets)))
            self.assertLess(max(targets), 56)

    def test_nominal_m1_intervals_on_stationary_simulated_data(self):
        rng = np.random.default_rng(21)
        y = rng.normal(5, 1, (76, 500))
        fit = fit_models(y[:60], Settings(process_ratio=0))
        score = metrics(forecast(fit, 'M1', np.arange(1, 17)), y[60:])
        self.assertGreater(score['coverage90'], .86)
        self.assertLess(score['coverage90'], .94)
        self.assertLess(abs(score['bias']), .1)

    def test_cache_reader_keeps_missing_dates_checks_qa_and_rejects_gap_filled_source(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            start = date(2025, 2, 5)
            dates = [(start + timedelta(days=d)).isoformat() for d in range(43)]
            save_json(root / 'manifest.json', {'status': 'complete_with_missing_dates',
                      'product': 'VNP46A2', 'collection': '002', 'tile': 'h09v05',
                      'start': dates[0], 'end': dates[-1]})
            save_json(root / 'inventory.json', {})
            pd.DataFrame([{'date': day, 'status': 'complete' if i == 0 else 'missing_source',
                           'source_filename': 'test.h5' if i == 0 else ''}
                          for i, day in enumerate(dates)]).to_csv(root / 'coverage.csv', index=False)
            directory = root / 'dyer' / dates[0]
            transform = from_origin(-90 + 64 / 240, 40 - 908 / 240, 1/240, 1/240)
            layers = {'radiance': np.array([[0., 2.], [3., np.nan]], dtype='float32'),
                      'strict_usable': np.array([[1, 0], [1, 0]], dtype='uint8'),
                      'cloud_mask': np.array([[48, 240], [48, 48]], dtype='uint16'),
                      'mandatory_quality': np.zeros((2, 2), dtype='uint8'),
                      'snow_flag': np.zeros((2, 2), dtype='uint8'),
                      'county_mask': np.ones((2, 2), dtype='uint8')}
            hashes = {}
            for name, values in layers.items():
                path = directory / f'{name}.tif'
                cog(path, values, transform, np.nan if name == 'radiance' else None,
                    source_layer='DNB_BRDF-Corrected_NTL' if name == 'radiance' else name)
                hashes[path.name] = digest(path)
            metadata = {'row_start': 908, 'column_start': 64, 'tile': 'h09v05', 'crs': 'EPSG:4326',
                        'transform': list(transform)[:6], 'date': dates[0], 'source_filename': 'test.h5',
                        'sha256': hashes}
            save_json(directory / 'metadata.json', metadata)
            # Deliberately unreadable alternate source: reader must not touch it.
            (directory / 'gap_filled_radiance.tif').write_bytes(b'not a raster')
            loaded = load_cache(root, start + timedelta(days=42))
            observed = loaded['arrays']['observed']
            self.assertEqual(observed.shape, (43, 4))
            np.testing.assert_equal(observed[0], [0, np.nan, 3, np.nan])
            self.assertTrue(np.isnan(observed[1:]).all())
            self.assertEqual(loaded['cells'].iloc[0]['cell_id'], 'h09v05_r0908_c0064')
            # Even with a valid updated hash, incorrect QA is rejected.
            mask_path = directory / 'strict_usable.tif'
            cog(mask_path, np.ones((2, 2), dtype='uint8'), transform, None)
            metadata['sha256'][mask_path.name] = digest(mask_path)
            save_json(directory / 'metadata.json', metadata)
            with self.assertRaisesRegex(ValueError, 'screening disagrees'):
                load_cache(root, start + timedelta(days=42))
            # A radiance filename alone cannot disguise the gap-filled source.
            path = directory / 'radiance.tif'
            cog(path, layers['radiance'], transform, np.nan, source_layer='Gap_Filled_DNB_BRDF-Corrected_NTL')
            metadata['sha256'][path.name] = digest(path)
            save_json(directory / 'metadata.json', metadata)
            with self.assertRaisesRegex(ValueError, 'non-gap-filled'):
                load_cache(root, start + timedelta(days=42))


if __name__ == '__main__':
    unittest.main()
