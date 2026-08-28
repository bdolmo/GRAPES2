import unittest

import numpy as np

from modules.plot import _genomewide_scatter_style, _genomewide_y_limits


class GenomewidePlotTestCase(unittest.TestCase):
    def test_exome_density_uses_small_transparent_markers(self):
        marker_size, alpha = _genomewide_scatter_style(200000)

        self.assertLessEqual(marker_size, 3.0)
        self.assertLessEqual(alpha, 0.2)

    def test_rare_extreme_outlier_does_not_flatten_display(self):
        ratios = np.concatenate([np.linspace(-0.8, 0.8, 10000), [25.0]])

        lower_limit, upper_limit = _genomewide_y_limits(ratios, -0.6, 0.4)

        self.assertLessEqual(lower_limit, -1.5)
        self.assertGreaterEqual(upper_limit, 1.5)
        self.assertLessEqual(upper_limit, 3.5)


if __name__ == "__main__":
    unittest.main()
