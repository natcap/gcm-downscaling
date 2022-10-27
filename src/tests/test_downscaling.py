import unittest

import numpy
import pandas
import xarray


class TestKNN(unittest.TestCase):
    """Tests knn.py."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_compute_delta_jp_matrices_prediction_dates(self):
        """Test compute_delta_jp_matrices raises ValueError."""
        from .. import knn

        lons = [0.0]
        lats = [20.0]
        gcm_start_date = '1980-01-01'
        gcm_end_date = '1990-12-31'
        times = pandas.date_range(gcm_start_date, gcm_end_date)

        lower_bound = 3
        upper_bound = 7
        data = numpy.reshape(
            numpy.random.uniform(
                low=0, high=10, size=len(times)*len(lons)*len(lats)),
            (len(lons), len(lats), len(times)))

        ds = xarray.Dataset(
            {
                'pr': (['lon', 'lat', 'time'], data)
            },
            coords={
                'lon': lons,
                'lat': lats,
                'time': times
            }
        )

        reference_start_date = gcm_start_date
        reference_end_date = '1985-01-01'

        prediction_start_date = '1986-01-01'
        prediction_end_date = '1995-01-01'  # after gcm_end_date
        with self.assertRaises(ValueError) as cm:
            knn.compute_delta_jp_matrices(
                ds, reference_start_date, reference_end_date,
                prediction_start_date, prediction_end_date,
                lower_bound, upper_bound)
        self.assertTrue(
            'is outside the time-range of the gcm' in str(cm.exception))

        prediction_start_date = '1970-01-01'  # before gcm_start_date
        prediction_end_date = '1990-01-01'
        with self.assertRaises(ValueError) as cm:
            knn.compute_delta_jp_matrices(
                ds, reference_start_date, reference_end_date,
                prediction_start_date, prediction_end_date,
                lower_bound, upper_bound)
        self.assertTrue(
            'is outside the time-range of the gcm' in str(cm.exception))


class TestJointProbability(unittest.TestCase):
    """Tests joint_probability.py."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_tri_state_joint_probability(self):
        """Test tri_state_joint_probability.py."""
        from .. import joint_probability

        timeseries = numpy.array([0, 0, 0, 5, 2, 1, 0, 0, 9, 8, 0])
        lower_bound = 1.5
        upper_bound = 7.5

        # This is the result from the R version of this function
        expected_matrix = numpy.array([
            [0.4, 0.1, 0.1],
            [0.1, 0.1, 0.0],
            [0.1, 0.0, 0.1]
        ])

        actual_matrix = joint_probability.tri_state_joint_probability(
            timeseries, lower_bound, upper_bound)

        numpy.testing.assert_array_equal(actual_matrix, expected_matrix)
