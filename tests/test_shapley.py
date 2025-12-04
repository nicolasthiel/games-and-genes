import unittest
import numpy as np

from games_and_genes.shapley import *


class test_shapley_values(unittest.TestCase):
    """Test suite for the shapley_values function."""

    def test_example1(self):
        """Example 1 in thesis paper."""

        coalitions = [
            {2},
            {0, 2, 3},
            {1, 2, 3}
        ]
        unanimity_coefficients = np.array([1/3, 1/3, 1/3])

        expected_shapley_values = np.array([1/9, 1/9, 5/9, 2/9])
        shapley_values = shapley_value(4, coalitions, unanimity_coefficients)

        np.testing.assert_array_equal(shapley_values, expected_shapley_values, "Example 1 shapley values are incorrect.")
        self.assertEqual(np.sum(shapley_values), 1.0, "Shapley values do not sum to v(N).")