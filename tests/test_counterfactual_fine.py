import unittest

import numpy as np
from scipy.stats import t as student_t

from scripts.counterfactual import Settings, fit_models, forecast, metrics
from scripts.run_counterfactual_fine import training_pool, score_sums


class FineCounterfactualTests(unittest.TestCase):
    def test_global_pool_and_chunked_predictions_match_whole_aoi(self):
        rng = np.random.default_rng(45)
        cube = rng.lognormal(1, .6, (56, 7, 9))
        cube[rng.random(cube.shape) < .3] = np.nan
        cube[:, 0, 0] = 0
        cube[:, 1, 1] = np.nan
        cube[:8, 1, 1] = 2
        mask = np.ones((7, 9), bool)
        mask[:, -1] = False
        cube[:, :, -1] = 1e9  # Halo must not enter pooling.
        settings = Settings()
        full = fit_models(cube[:, mask], settings)
        for rows in (1, 3, 7):
            pool = training_pool(cube, mask, 56, settings, rows)
            np.testing.assert_equal(pool[0], full['pool_edges'])
            np.testing.assert_equal(pool[1], full['pool_scales'])
            chunks = [fit_models(cube[:, mask][:, start:start + 11], settings, scale_pool=pool)
                      for start in range(0, int(mask.sum()), 11)]
            for method in ('M0', 'M1'):
                expected = forecast(full, method, [1, 7, 16])
                for key, value in expected.items():
                    actual = np.concatenate([forecast(f, method, [1, 7, 16])[key] for f in chunks], axis=1)
                    np.testing.assert_equal(actual, value)
                np.testing.assert_equal(expected['q05'], expected['center'] + expected['scale'] * student_t.ppf(.05, expected['df']))
                np.testing.assert_equal(expected['q95'], expected['center'] + expected['scale'] * student_t.ppf(.95, expected['df']))

    def test_pool_uses_only_training_and_handles_no_supported_pixels(self):
        cube = np.ones((56, 2, 3))
        mask = np.ones((2, 3), bool)
        before = training_pool(cube, mask, 28, Settings(), 1)
        cube[28:] = 1000
        after = training_pool(cube, mask, 28, Settings(), 2)
        np.testing.assert_equal(before, after)
        cube[:] = np.nan
        edges, scales = training_pool(cube, mask, 28, Settings(), 1)
        self.assertTrue(np.isfinite(edges).all())
        np.testing.assert_equal(scales, np.full(4, .1))
        fit = fit_models(cube.reshape(56, -1), scale_pool=(edges, scales))
        self.assertFalse(fit['supported'].any())

    def test_additive_validation_scores_match_existing_metrics(self):
        rng = np.random.default_rng(4)
        y = rng.normal(3, 1, (42, 20))
        y[34:, :3] = np.nan
        fit = fit_models(y[:28])
        pred = forecast(fit, 'M1', np.arange(1, 15))
        reference = metrics(pred, y[28:])
        totals = score_sums(pred, y[28:], fit['baseline_median'])['all']
        self.assertEqual(totals['n'], reference['n'])
        for name, value in totals.items():
            if name == 'n':
                continue
            if name == 'mse':
                self.assertAlmostEqual(np.sqrt(value / totals['n']), reference['rmse'])
            else:
                self.assertAlmostEqual(value / totals['n'], reference[name])


if __name__ == '__main__':
    unittest.main()
