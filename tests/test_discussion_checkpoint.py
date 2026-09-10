import unittest

import numpy as np

from scripts.discussion_checkpoint import missing_kmeans


class DiscussionCheckpointTests(unittest.TestCase):
    def test_missingness_does_not_become_a_zero_signal(self):
        y = np.array([[-10, np.nan, -10, np.nan], [np.nan, -10, np.nan, -10],
                      [10, np.nan, 10, np.nan], [np.nan, 10, np.nan, 10],
                      [np.nan, np.nan, np.nan, np.nan], [3, np.nan, np.nan, np.nan]])
        labels, info = missing_kmeans(y, k=2, min_dates=2)
        np.testing.assert_equal(labels, [0, 0, 1, 1, -1, -1])
        self.assertTrue(info['converged'])
        self.assertEqual(info['eligible_pixels'], 4)
        # Adding an entirely missing date must have no effect.
        expanded, _ = missing_kmeans(np.column_stack([y, np.full(len(y), np.nan)]), k=2, min_dates=2)
        np.testing.assert_equal(labels, expanded)

    def test_chunk_size_does_not_change_separated_clusters(self):
        rng = np.random.default_rng(2)
        y = np.concatenate([rng.normal(level, .05, (30, 6)) for level in (-4, 0, 4)])
        y[rng.random(y.shape) < .15] = np.nan
        a, _ = missing_kmeans(y, k=3, chunk=7)
        b, _ = missing_kmeans(y, k=3, chunk=100)
        np.testing.assert_equal(a, b)

    def test_no_eligible_observations_remain_unclassified(self):
        labels, info = missing_kmeans(np.full((10, 9), np.nan))
        np.testing.assert_equal(labels, np.full(10, -1))
        self.assertEqual(info['clusters'], 0)


if __name__ == '__main__':
    unittest.main()
