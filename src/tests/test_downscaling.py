import unittest

import numpy
import xarray


class TestKNN(unittest.TestCase):
    """Tests knn.py."""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_tri_state_joint_probability(self):
        """Test tri_state_joint_probability.py."""
        from .. import knn

        timeseries = numpy.array([0, 0, 0, 5, 2, 1, 0, 0, 9, 8, 0])
        lower_bound = 1.5
        upper_bound = 7.5

        # This is the result from the R version of this function
        expected_matrix = numpy.array([
            [0.4, 0.1, 0.1],
            [0.1, 0.1, 0.0],
            [0.1, 0.0, 0.1]
        ])

        actual_matrix = knn.tri_state_joint_probability(
            timeseries, lower_bound, upper_bound)

        numpy.testing.assert_array_equal(actual_matrix, expected_matrix)

    def test_slice_dates_around_dayofyear(self):
        """"""
        from .. import knn

        dates_index = knn.date_range_no_leap('1979-01-01', '1981-12-31')
        month = 10
        day = 16
        near_window = 10
        idx = knn.slice_dates_around_dayofyear(
            dates_index, month, day, near_window)

        # Given the size of the window, and the day, all dates
        # in the index should be in the same month:
        self.assertEqual(
            numpy.count_nonzero(dates_index[idx].month == month),
            len(idx))

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

