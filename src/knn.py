from collections import defaultdict
import datetime
import logging
import os
import sys

import numpy
from numpy.lib.stride_tricks import sliding_window_view
import pandas
import xarray

LOGGER = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
logging.basicConfig(level=logging.INFO, handlers=[handler])

PERCENTILE_THRESHOLDS = (25, 75)
DRY = 0  # joint-probability matrices are (3, 3) and index at 0.
WET = 1
VERY_WET = 2
TRANSITION_TYPES = numpy.array([
    ['AA', 'AB', 'AC'],
    ['BA', 'BB', 'BC'],
    ['CA', 'CB', 'CC']
], dtype="U2")
DAYS_IN_YEAR = 365
NEAR_WINDOW = 15  # days
REF_PERIOD_START_DATE = '1985-01-01'
REF_PERIOD_END_DATE = '2014-12-31'
PREDICTION_PERIOD_START_DATE = '2015-01-01'
PREDICTION_PERIOD_END_DATE = '2016-01-31'
DATA_STORE_PATH = 'H://Shared drives/GCM_Climate_Tool/required_files'


def shift_longitude_from_360(dataset):
    """Reassign longitude coordinates from 0:360 scale to -180:180.

    Args:
        dataset (xarray.Dataset): an in-memory representation of a netCDF
            with a `coord` named `lon`.

    Raises:
        ValueError if the dataset does not include a `lon` coordinate.

    Returns:
        xarray.Dataset
    """
    if 'lon' in dataset.dims:
        return dataset.assign_coords(lon=(((dataset.lon + 180) % 360) - 180))
    raise ValueError('Could not reassign longitude coordinates,'
                     'No dimension "lon" in ')


def date_range_no_leap(start, stop):
    """Create a series of dates from start to stop, skipping leap-days.

    calendar='noleap' always creates 365-day years. It also forces
    cftime-format dates, which is not so convenient. isoformat is more
    interoperable, and pandas DatetimeIndex is most convenient for indexing.

    If start & stop include a leap year, xarray will still skip
    the leap-day, but pandas DatetimeIndex will number the date.dayofyear
    from 0 - 366, for the leap year, skipping over day 60.

    Args:
        start (string): of format 'YYYY-MM-DD'
        stop (string): of format 'YYYY-MM-DD'

    Returns:
        (pandas.DatetimeIndex) series of dates from start to stop.
    """
    return pandas.DatetimeIndex(
        [d.isoformat() for d in
            xarray.date_range(start, stop, calendar='noleap')])


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


def state_transition_series(array, lower_bound, upper_bound):
    def func(a):
        if a[0] <= lower_bound:
            if a[1] <= lower_bound:
                return 'AA'
            elif a[1] > upper_bound:
                return 'AC'
            return 'AB'
        elif a[0] > upper_bound:
            if a[1] <= lower_bound:
                return 'CA'
            elif a[1] > upper_bound:
                return 'CC'
            return 'CB'
        if a[1] <= lower_bound:
            return 'BA'
        elif a[1] > upper_bound:
            return 'BC'
        return 'BB'

    pairs = sliding_window_view(array, window_shape=2)
    transitions = numpy.apply_along_axis(func, 1, pairs)
    return transitions


def slice_dates_around_dayofyear(dates_index, month, day, near_window):
    """Get index of dates from within an n-day window of a given day-of-year.

    For example, get the index of all dates from all years included in
    `dates_index`, within `near_window` days of July 1st.

    Sometimes the `dates_index` comes from a historical record, which
    could have a 366-day year. So rather than getting a window
    around a day-of-year (e.g. day 65), we always get a window
    around a month-day date (e.g. March 3).

    Args:
        dates_index (pandas.DatetimeIndex): daily frequency
        month (integer): from 1 - 12
        day (integer): day of the month
        near_window (integer): number of days before & after month-day

    Returns:
        (numpy.array): a 1-d array of indexes of size `dates_index`
    """
    idx = numpy.nonzero(
        (dates_index.month == month)
        & (dates_index.day == day))[0]
    # mask to exclude dates with window extending beyond reference period
    mask = (idx - near_window >= 0) & (idx + near_window < dates_index.size)
    ranges = [range(x - near_window, x + near_window + 1) for x in idx[mask]]
    return numpy.array([i for r in ranges for i in r])


def marginal_probability_of_transitions(observed_matrix, delta_matrix):
    """Add joint-probability matrices; calculate marginal probability matrix.

    Args:
        observed_matrix (numpy.array): 2-d array of shape (3,3)
        delta_matrix (numpy.array): 2-d array of shape (3,3)

    Returns:
        (numpy.array): 2-d array of shape (3,3)
    """
    projected_jp_matrix = observed_matrix + delta_matrix
    projected_jp_matrix[projected_jp_matrix < 0] = 0

    def func(x):
        if x.sum() > 0:
            return x / x.sum()
        return x

    return numpy.apply_along_axis(func, 1, projected_jp_matrix)


def jp_matrix_from_transitions_sequence(
        dates_index, transitions_array, month, day, near_window):
    flat_idx = slice_dates_around_dayofyear(
        dates_index, month, day, near_window)
    frequencies = [numpy.count_nonzero(transitions_array[flat_idx] == x)
                   for x in TRANSITION_TYPES.flatten()]
    proportions = numpy.array(frequencies) / len(flat_idx)
    return proportions.reshape(3, 3)


def compute_historical_jp_matrices(
        transitions_array, reference_period_date_index):
    LOGGER.info(
        f'computing JP matrices for historical reference period '
        f'{reference_period_date_index.min()} : {reference_period_date_index.max()}')
    # Need a calendar for an arbitrary year
    # to collect all Jan 1sts, Jan 2nds, etc in the reference period
    caldates = date_range_no_leap('1980-01-01', '1980-12-31')
    # We don't know the transition for the final day in the reference period
    reference_period_date_index = reference_period_date_index[:-1]
    assert(len(reference_period_date_index) == len(transitions_array))

    jp_matrix_dict = defaultdict(dict)
    for date in caldates:
        jp_matrix_dict[date.month][date.day] = jp_matrix_from_transitions_sequence(
            reference_period_date_index, transitions_array,
            date.month, date.day, NEAR_WINDOW)

    return jp_matrix_dict


def downscale_precipitation(
        observed_data_path, var, simulation_dates_index, reference_start_date,
        reference_end_date, gcm_dataset):
    dates_lookup = {}

    obs_df = pandas.read_csv(
        observed_data_path,
        usecols=['Year', 'Month', 'Day', var],
        parse_dates={'date': ['Year', 'Month', 'Day']})
    obs_df.set_index('date', inplace=True)
    obs_df = obs_df[reference_start_date: reference_end_date]

    lower_bound_obs, upper_bound_obs = numpy.percentile(
        obs_df[var].values, q=PERCENTILE_THRESHOLDS)
    transitions_array_obs = state_transition_series(
        obs_df[var].values, lower_bound_obs, upper_bound_obs)
    historic_obs_jp_matrix_lookup = compute_historical_jp_matrices(
        transitions_array_obs, obs_df.index)

    gcm_reference_period_ds = gcm_dataset.sel(
        time=slice(reference_start_date, reference_end_date))
    lower_bound_gcm, upper_bound_gcm = numpy.percentile(
        gcm_reference_period_ds.pr.values, q=PERCENTILE_THRESHOLDS)
    transitions_array_gcm = state_transition_series(
        gcm_reference_period_ds.pr.values, lower_bound_gcm, upper_bound_gcm)
    reference_period_date_index = pandas.DatetimeIndex(
        gcm_reference_period_ds.time.values)
    historic_gcm_jp_matrix_lookup = compute_historical_jp_matrices(
        transitions_array_gcm, reference_period_date_index)

    day_one = simulation_dates_index[0]
    window_idx = slice_dates_around_dayofyear(
            obs_df.index, day_one.month, day_one.day, NEAR_WINDOW)

    a_date = pandas.Timestamp(obs_df.index[numpy.random.choice(window_idx)])
    for sim_date in simulation_dates_index:
        # TODO: don't need to re-determine this here, we know the wetstate
        # because we chose it deliberately. Though having this revealed
        # a bug in our indexing....
        precip = obs_df[obs_df.index == a_date][var].values[0]
        current_wet_state = WET
        if precip <= lower_bound_obs:
            current_wet_state = DRY
        if precip > upper_bound_obs:
            current_wet_state = VERY_WET

        # TODO: not sure we need to track all this, but for diagnostic purposes
        sim_dict = {
            'historic_date': a_date,
            'historic_precip': precip,
            'wet_state': current_wet_state,
        }

        date_offset = pandas.DateOffset(days=NEAR_WINDOW)
        window_start = sim_date - date_offset
        window_end = sim_date + date_offset
        array = gcm_dataset.sel(time=slice(window_start, window_end)).pr.values
        # TODO: Is it okay for the window to overlap the historical period?
        # This occurs in the first NEAR_WINDOW days of the prediction
        # period.

        jp_matrix = tri_state_joint_probability(
            array, lower_bound_gcm, upper_bound_gcm)
        gcm_ref_jp_matrix = historic_gcm_jp_matrix_lookup[sim_date.month][sim_date.day]
        delta_jp_doy_yr = jp_matrix - gcm_ref_jp_matrix

        observed_jp_doy = historic_obs_jp_matrix_lookup[sim_date.month][sim_date.day]
        margins_matrix = marginal_probability_of_transitions(
            observed_jp_doy, delta_jp_doy_yr)
        next_wet_state = numpy.random.choice(
            [DRY, WET, VERY_WET], p=margins_matrix[current_wet_state])
        sim_dict['next_wet_state'] = next_wet_state

        valid_mask = numpy.array([False])
        search_window = NEAR_WINDOW
        while not valid_mask.any():
            # transitions array is one day shorter than historical record
            # so trim the last day off the historical record before indexing.
            window_idx = slice_dates_around_dayofyear(
                obs_df.index[:-1], sim_date.month, sim_date.day, search_window)
            valid_mask = transitions_array_obs[window_idx] == \
                TRANSITION_TYPES[current_wet_state][next_wet_state]
            search_window += 1
        if search_window > NEAR_WINDOW + 1:
            LOGGER.info(
                f'Needed a window of {search_window - 1} days to find a '
                f'transition matching {current_wet_state}->{next_wet_state} '
                f'for simulation date {sim_date}')

        # Which leading-day from the matching transitions is most similar to
        # the current day's precip?
        neighbors = obs_df.iloc[window_idx[valid_mask], ]
        distances = numpy.abs(neighbors[var].values - precip)
        nearest = distances.argmin()  # TODO: what if there's a tie?
        # we want the day after the nearest-neighbor date, but it's not safe
        # to just add 1 day, we might hit Feb 29th, which might not exist in
        # observed record. So find the index of the NN date from the observed
        # record, and get the subsequent record's date.
        chosen_idx = obs_df.index.get_loc(neighbors.index[nearest]) + 1
        a_date = obs_df.index[chosen_idx]
        sim_dict['next_historic_date'] = a_date

        dates_lookup[sim_date] = sim_dict

    dataframe = pandas.DataFrame.from_dict(dates_lookup, orient='index')
    dataframe.to_csv('bootstrapped_dates.csv')


if __name__ == "__main__":
    precip_var = 'pr_Beni01'
    site_name = 'Beni01'
    observed_precip_path = os.path.join(
        DATA_STORE_PATH, 'OBSERVATIONS/LLdM_AOI2/series_pr_diario.csv')
    observed_location = os.path.join(
        DATA_STORE_PATH, 'OBSERVATIONS/LLdM_AOI2/catalogo.csv')
    historical_gcm_path = os.path.join(
        DATA_STORE_PATH, 'GCMs',
        'Amazon__pr_day_CanESM5_historical_r1i1p1f1_gn_18500101-20141231.nc')
    future_gcm_path = os.path.join(
        DATA_STORE_PATH, 'GCMs',
        'Amazon__pr_day_CanESM5_ssp126_r1i1p1f1_gn_20150101-21001231.nc')

    locations = pandas.read_csv(observed_location)
    lon, lat = locations[locations['CODIGO_CAT'] == site_name][['longitud', 'latitud']].values[0]

    simulation_dates_index = date_range_no_leap(
        PREDICTION_PERIOD_START_DATE,
        PREDICTION_PERIOD_END_DATE)
    LOGGER.info(
        f'Simulating for period {PREDICTION_PERIOD_START_DATE} : '
        f'{PREDICTION_PERIOD_END_DATE}')

    with xarray.open_mfdataset(
            [historical_gcm_path, future_gcm_path]) as mfds:

        gcm_start_date = mfds.time.values.min()
        gcm_end_date = mfds.time.values.max()
        date_offset = datetime.timedelta(days=NEAR_WINDOW)
        if ((pandas.Timestamp(PREDICTION_PERIOD_END_DATE) + date_offset) > gcm_end_date or
                (pandas.Timestamp(PREDICTION_PERIOD_START_DATE) - date_offset) < gcm_start_date):
            raise ValueError(
                f'the requested prediction period {PREDICTION_PERIOD_START_DATE} : '
                f'{PREDICTION_PERIOD_END_DATE} is outside the time-range of the gcm'
                f'({gcm_start_date} : {gcm_end_date})')

        mfds = shift_longitude_from_360(mfds)

        downscale_precipitation(
            observed_precip_path,
            precip_var,
            simulation_dates_index,
            REF_PERIOD_START_DATE,
            REF_PERIOD_END_DATE,
            mfds.sel(lon=lon, lat=lat, method='nearest'))
