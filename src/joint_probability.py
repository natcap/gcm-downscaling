import numpy


def tri_state_joint_probability(timeseries, lower_bound, upper_bound):
    """Calculate probabilities of state transitions for consecutive
       values in timeseries.

    Args:
        timeseries (numpy.array): a 1-dimensional array.
        lower_bound (float): the lower boundary of the middle bin.
        upper_bound (float): the upper boundary (inclusive) of the middle bin.

    Returns:
        (numpy.array): a 2-dimensional array of joint probabilities.
    """
    jp_matrix = numpy.zeros([3, 3])

    a = timeseries[:-1]
    b = timeseries[1:]

    a_low_mask = a <= lower_bound
    jp_matrix[0, 0] = numpy.count_nonzero(b[a_low_mask] <= lower_bound)
    jp_matrix[0, 2] = numpy.count_nonzero(b[a_low_mask] > upper_bound)
    jp_matrix[0, 1] = numpy.count_nonzero(
        numpy.count_nonzero(a_low_mask) - jp_matrix[0, 0] - jp_matrix[0, 2])

    a_high_mask = a > upper_bound
    jp_matrix[2, 0] = numpy.count_nonzero(b[a_high_mask] <= lower_bound)
    jp_matrix[2, 2] = numpy.count_nonzero(b[a_high_mask] > upper_bound)
    jp_matrix[2, 1] = numpy.count_nonzero(
        numpy.count_nonzero(a_high_mask) - jp_matrix[2, 0] - jp_matrix[2, 2])

    a_med_mask = (a > lower_bound) & (a <= upper_bound)
    jp_matrix[1, 0] = numpy.count_nonzero(b[a_med_mask] <= lower_bound)
    jp_matrix[1, 2] = numpy.count_nonzero(b[a_med_mask] > upper_bound)
    jp_matrix[1, 1] = numpy.count_nonzero(
        numpy.count_nonzero(a_med_mask) - jp_matrix[1, 0] - jp_matrix[1, 2])

    return (jp_matrix / numpy.sum(jp_matrix))
