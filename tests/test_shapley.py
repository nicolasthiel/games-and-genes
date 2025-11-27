import unittest
import numpy as np

from python.shapley import boolean_diff_expressed_matrix
from python.shapley import support_of_binary_matrix

class test_boolean_diff_expressed_matrix(unittest.TestCase):
    """
    Test suite for the boolean_diff_expressed_matrix function.
    Uses numpy.testing.assert_array_equal for robust array comparisons.
    """

    def setUp(self):
        """Set up standard input arrays for testing."""
        # A_SR (Reference): Used to define the percentile bounds
        self.A_SR = np.array([
            [10, 20, 30, 40],
            [1, 2, 3, 100],
            [5, 5, 5, 5]
        ])
        
        # A_SD (Diseased): Used to test against the bounds
        self.A_SD = np.array([
            [5, 20, 45, 100],  
            [0, 50, 100, 1],   
            [4, 5, 6, 7]       
        ])

    def test_default_100_percentile_range(self):
        """
        Test the default 0th to 100th percentile range.
        The bounds are the MIN and MAX of A_SR for each row.
        B is True if A_SD value is strictly outside the MIN/MAX range of A_SR.
        """
        # Default: lower_prctl=0, upper_prctl=100
        B = boolean_diff_expressed_matrix(self.A_SD, self.A_SR)

        # Expected Bounds based on A_SR:
        # R0 Bounds: [10, 40]
        # R1 Bounds: [1, 100]
        # R2 Bounds: [5, 5] (Since 5 <= 5 and 5 >= 5 are both True, 5 is 'outside' in this strict case)

        expected_B = np.array([
            # A_SD [5, 20, 45, 100] vs [10, 40]
            [True, False, True, True],
            # A_SD [0, 50, 100, 1] vs [1, 100]
            [True, False, True, True],
            # A_SD [4, 5, 6, 7] vs [5, 5]
            [True, True, True, True]
        ])

        np.testing.assert_array_equal(B, expected_B, "Test failed for 0-100 percentile bounds.")

    def test_interquartile_range_check(self):
        """
        Test a middle-range percentile (25th to 75th).
        B should be True for values outside the IQR.
        """
        lower_prctl = 25
        upper_prctl = 75
        B = boolean_diff_expressed_matrix(self.A_SD, self.A_SR, lower_prctl, upper_prctl)

        # Expected Bounds (calculated on A_SR):
        # R0 Bounds: [17.5, 32.5]
        # R1 Bounds: [1.75, 27.25]
        # R2 Bounds: [5.0, 5.0]

        expected_B = np.array([
            # R0 [5, 20, 45, 100] vs [17.5, 32.5]
            [True, False, True, True],
            # R1 [0, 50, 100, 1] vs [1.75, 27.25]
            [True, True, True, True],
            # R2 [4, 5, 6, 7] vs [5.0, 5.0] (As 5 <= 5 and 5 >= 5, all are True)
            [True, True, True, True]
        ])

        np.testing.assert_array_equal(B, expected_B, "Test failed for 25-75 percentile bounds.")

    def test_single_median_bound_edge_case(self):
        """
        Test the extreme edge case where lower_prctl == upper_prctl (e.g., 50).
        Since the function uses <= and >=, any value in A_SD will satisfy 
        (A_SD <= median) OR (A_SD >= median), resulting in an all-True matrix.
        """
        lower_prctl = 50
        upper_prctl = 50
        B = boolean_diff_expressed_matrix(self.A_SD, self.A_SR, lower_prctl, upper_prctl)
        
        # The result must be all True due to the inclusive inequality operators (<= and >=).
        expected_B = np.full(self.A_SD.shape, True)
        
        np.testing.assert_array_equal(B, expected_B, "Test failed for 50-50 percentile bounds.")

    def test_output_properties(self):
        """Ensure the output array has the correct shape and is a boolean type."""
        B = boolean_diff_expressed_matrix(self.A_SD, self.A_SR)
        self.assertEqual(B.shape, self.A_SD.shape, "Output shape is incorrect.")
        self.assertEqual(B.dtype, np.dtype(bool), "Output type is not boolean.")

    def test_example1(self):
        """Example 1 in thesis paper."""
        A_SR = np.array([
            [3, 1, 0.3],
            [2, 12, 4],
            [6, 1.3, 3],
            [0.8, 2.5, 2]
        ])
        A_SD = np.array([
            [1, 3, 0.7],
            [4.1, 8.5, 1.0],
            [6, 25, 6],
            [1, 0.6, 0.8]
        ])
        B = boolean_diff_expressed_matrix(A_SD, A_SR)
        expected_B = np.array([
            [False, True, False],
            [False, False, True],
            [True, True, True],
            [False, True, True]
        ])
        np.testing.assert_array_equal(B, expected_B, "Example 1 test case failed.")


class test_support_of_binary_matrix(unittest.TestCase):
    """
    Test suite for the support_of_binary_matrix function.
    Uses numpy.testing.assert_array_equal for robust array comparisons.
    """

    def test_support_basic(self):
        """Test basic functionality of support_of_binary_matrix."""
        from python.shapley import support_of_binary_matrix

        B = np.array([
            [True, False, True],
            [False, True, False],
            [True, True, False]
        ])

        expected_support = [
            {0, 2},  # Column 0 has True at rows 0 and 2
            {1, 2},  # Column 1 has True at rows 1 and 2
            {0}      # Column 2 has True at row 0
        ]

        support = support_of_binary_matrix(B)
        self.assertEqual(support, expected_support, "Support calculation is incorrect.")

    def test_support_all_false_matrix(self):
        """Test support_of_binary_matrix with an all-false matrix."""
        from python.shapley import support_of_binary_matrix

        B = np.array([
            [False, False, False],
            [False, False, False],
            [False, False, False]
        ])
        expected_support = [set(), set(), set()] # list of empty sets

        support = support_of_binary_matrix(B)
        self.assertEqual(support, expected_support, "Support for all-false matrix should be a list of empty sets.")

    def test_example1(self):
        """Example 1 in thesis paper."""
        B = np.array([
            [False, True, False],
            [False, False, True],
            [True, True, True],
            [False, True, True]
        ])
        expected_support = [
            {2},      # Column 0 has True at row 2
            {0, 2, 3},# Column 1 has True at rows 0, 2, and 3
            {1, 2, 3} # Column 2 has True at rows 1, 2, and 3
        ]

        support = support_of_binary_matrix(B)
        self.assertEqual(support, expected_support, "Example 1 support calculation is incorrect.")