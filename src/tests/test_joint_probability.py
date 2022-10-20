import unittest
import numpy


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
