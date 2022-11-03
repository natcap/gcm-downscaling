import unittest

import numpy
import xarray


class TestKNN(unittest.TestCase):
    """Tests knn.py."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_state_transition_sequence(self):
        """Test create sequence of state transitions."""
        from .. import knn

        array = numpy.array([0, 0, 0, 5, 2, 1, 0, 0, 9, 8, 0])
        expected_array = [
            'AA', 'AA', 'AB', 'BA', 'AA',
            'AA', 'AA', 'AC', 'CC', 'CA']
        lower_bound = 3
        upper_bound = 7
        transitions = knn.state_transition_table(
            array, lower_bound, upper_bound)
        numpy.testing.assert_array_equal(transitions, expected_array)

    def test_jp_matrix_from_transitions_sequence(self):
        """Test joint probablity matrix from sequence of transitions."""
        from .. import knn

        n = 900
        reference_start_date = '1980-01-01'
        reference_end_date = '1989-12-31'
        ref_dates = xarray.date_range(
            reference_start_date, reference_end_date, calendar='noleap')
        ref_dates = ref_dates[:n]
        month = 6
        day = 20
        near_window = 4

        sample_transitions = [
            'AA', 'AB', 'AC', 'BA', 'BB', 'BC', 'CA', 'CB', 'AA']
        transitions = numpy.empty(ref_dates.shape, dtype='U2')
        transitions[:] = 100 * sample_transitions

        # Our sample_transitions pattern repeats every 9 days, and a full
        # window sequence is 9 days, so the expected frequencies of
        # transitions in the computed matrix match the frequencies of
        # occurrence in sample_transitions:
        expected_matrix = numpy.array([
            [0.222222, 0.111111, 0.111111],
            [0.111111, 0.111111, 0.111111],
            [0.111111, 0.111111, 0.0]
        ])
        jp_matrix = knn.jp_matrix_from_transitions_sequence(
            ref_dates, transitions, month, day, near_window)
        numpy.testing.assert_array_almost_equal(jp_matrix, expected_matrix)

    def test_compute_delta_jp_matrices_prediction_dates(self):
        """Test compute_delta_jp_matrices raises ValueError."""
        from .. import knn

        lons = [0.0]
        lats = [20.0]
        gcm_start_date = '1980-01-01'
        gcm_end_date = '1990-12-31'
        times = xarray.date_range(
            gcm_start_date, gcm_end_date, use_cftime=True)

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
                ds.isel(lon=0, lat=0), reference_start_date, reference_end_date,
                prediction_start_date, prediction_end_date,
                lower_bound, upper_bound)
        self.assertTrue(
            'is outside the time-range of the gcm' in str(cm.exception))

        prediction_start_date = '1970-01-01'  # before gcm_start_date
        prediction_end_date = '1990-01-01'
        with self.assertRaises(ValueError) as cm:
            knn.compute_delta_jp_matrices(
                ds.isel(lon=0, lat=0), reference_start_date, reference_end_date,
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
