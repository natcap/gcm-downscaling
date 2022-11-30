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

    def test_shift_longitude_error_on_missing_dimension(self):
        """Test shift longitude from 0-360."""
        from .. import knn

        precip = 10 * numpy.random.rand(2, 2)
        lon = [[-99.83, -99.32], [-99.79, -99.23]]
        lat = [[42.25, 42.21], [42.63, 42.59]]
        dataset = xarray.Dataset(
            {
                'pr': (['x', 'y'], precip),
            },
            coords={
                'long': (['x', 'y'], lon),  # deliberately mislabeled coord.
                'lat': (['x', 'y'], lat),
            }
        )

        with self.assertRaises(ValueError) as cm:
            _ = knn.shift_longitude_from_360(dataset)
        self.assertTrue(
            'Expected dimension "lon" but found coordinates' in str(cm.exception))

    def test_shift_longitude_correct_given_180(self):
        """Test shift longitude is correct given coords from -180-180."""
        from .. import knn

        precip = 10 * numpy.random.rand(2, 2)
        lon = [[-180.0, -5.0], [0.0, 179.9]]
        lat = [[42.25, 42.21], [42.63, 42.59]]
        dataset = xarray.Dataset(
            {
                'pr': (['x', 'y'], precip),
            },
            coords={
                'lon': (['x', 'y'], lon),
                'lat': (['x', 'y'], lat),
            }
        )

        shifted_dataset = knn.shift_longitude_from_360(dataset)
        numpy.testing.assert_array_almost_equal(
            shifted_dataset.lon.to_numpy(), lon)

    def test_tri_state_joint_probability(self):
        """Test tri_state_joint_probability."""
        from .. import knn

        timeseries = numpy.array(
            [0, 0, 2, 5, 2, 1, 0, 0, 9, 8, 0])
        lower_bound = 1.5
        upper_bound = 7.5

        expected_matrix = numpy.array([
            [0.3, 0.1, 0.1],
            [0.1, 0.2, 0.0],
            [0.1, 0.0, 0.1]
        ])

        actual_matrix = knn.tri_state_joint_probability(
            timeseries, lower_bound, upper_bound)

        numpy.testing.assert_array_equal(actual_matrix, expected_matrix)

    def test_slice_dates_around_dayofyear(self):
        """Test slicing dates from an xarray.Dataset time coordinate"""
        from .. import knn

        dates_index = knn.date_range_no_leap('1979-01-01', '1981-12-31')
        dataset = xarray.Dataset({'time': dates_index})
        month = 10
        day = 16
        near_window = 4
        expected_dates = [
            '1979-10-12', '1979-10-13', '1979-10-14', '1979-10-15',
            '1979-10-16', '1979-10-17', '1979-10-18', '1979-10-19', '1979-10-20',
            '1980-10-12', '1980-10-13', '1980-10-14', '1980-10-15',
            '1980-10-16', '1980-10-17', '1980-10-18', '1980-10-19', '1980-10-20',
            '1981-10-12', '1981-10-13', '1981-10-14', '1981-10-15',
            '1981-10-16', '1981-10-17', '1981-10-18', '1981-10-19', '1981-10-20',
        ]
        idx = knn.slice_dates_around_dayofyear(
            dataset.time, month, day, near_window)
        self.assertListEqual(
            [d.strftime(format='%Y-%m-%d') for d in dates_index[idx]],
            expected_dates)

    def test_state_transition_series(self):
        """Test create series of state transitions."""
        from .. import knn

        array = numpy.array([0, 0, 0, 5, 2, 1, 0, 0, 9, 8, 0])
        expected_array = [
            'AA', 'AA', 'AB', 'BA', 'AA',
            'AA', 'AA', 'AC', 'CC', 'CA']
        lower_bound = 3
        upper_bound = 7
        transitions = knn.state_transition_series(
            array, lower_bound, upper_bound)
        numpy.testing.assert_array_equal(transitions, expected_array)

    def test_jp_matrix_from_transitions_sequence(self):
        """Test joint probablity matrix from sequence of transitions."""
        from .. import knn

        n = 900
        reference_start_date = '1980-01-01'
        reference_end_date = '1989-12-31'
        ref_dates = knn.date_range_no_leap(
            reference_start_date, reference_end_date)
        ref_dates = ref_dates[:n]
        dataset = xarray.Dataset({'time': ref_dates})
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
            dataset.time, transitions, month, day, near_window)
        numpy.testing.assert_array_almost_equal(jp_matrix, expected_matrix)

    def test_marginal_probabilities(self):
        """Test marginal probability calc when JP shifts to <= 0."""
        from .. import knn

        # These matrices occurred when running the program.
        observed_jp_matrix = numpy.array([
            [0.00222469, 0.01001112, 0.        ],
            [0.00889878, 0.3003337 , 0.15906563],
            [0.00111235, 0.16129032, 0.3570634 ]])

        # These changes in probabilities will shift the entire first
        # row of observed joint-probabilities to <= 0.
        delta_jp_matrix = numpy.array([
            [-0.05576567, -0.02788283,  0.        ],
            [ 0.16544308, -0.45057471,  0.09885057],
            [ 0.        ,  0.02662217,  0.24330738]])

        result = knn.marginal_probability_of_transitions(
            observed_jp_matrix + delta_jp_matrix)

        sums = numpy.sum(result, axis=1)
        self.assertTrue(numpy.all(sums > 0))
        self.assertTrue(numpy.all(numpy.isclose(sums, sums[0])))
