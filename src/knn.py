from collections import defaultdict
import datetime
import logging
import os

import numpy
from numpy.lib.stride_tricks import sliding_window_view
import pandas
import pygeoprocessing
import xarray

from .joint_probability import tri_state_joint_probability

LOGGER = logging.getLogger(__name__)

LOWER_BOUND = 2.0e-06  # units of the GCM pr variable
UPPER_BOUND = 4.0e-06  # upper is included in middle bin, lower is not
DRY = 1
WET = 2
VERY_WET = 3
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


def state_transition_table(array, lower_bound, upper_bound):
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


def jp_matrix_from_transitions_sequence(
        dates_index, transitions_array, month, day, near_window):
    idx = numpy.nonzero(
        (dates_index.month == month)
        & (dates_index.day == day))[0]
    # mask to exclude dates with window extending beyond reference period
    mask = (idx - near_window >= 0) & (idx + near_window < dates_index.size)
    ranges = [range(x - near_window, x + near_window + 1) for x in idx[mask]]
    flat_idx = [i for r in ranges for i in r]
    frequencies = [numpy.count_nonzero(transitions_array[flat_idx] == x)
                   for x in TRANSITION_TYPES.flatten()]
    proportions = numpy.array(frequencies) / len(flat_idx)
    return proportions.reshape(3, 3)


def compute_historical_jp_matrices(
        transitions_array, reference_period_time_array):
    LOGGER.info(
        f'computing JP matrices for historical reference period '
        f'{reference_period_time_array.min()} : {reference_period_time_array.max()}')
    # Need a calendar for an arbitrary year
    # to collect all Jan 1sts, Jan 2nds, etc in the reference period
    caldates = date_range_no_leap('1980-01-01', '1980-12-31')
    # We don't know the transition for the final day in the reference period
    reference_period_time_array = reference_period_time_array[:-1]
    assert(len(reference_period_time_array) == len(transitions_array))

    jp_matrix_dict = defaultdict(dict)
    dates_index = pandas.DatetimeIndex(reference_period_time_array)
    for date in caldates:
        jp_matrix_dict[date.month][date.day] = jp_matrix_from_transitions_sequence(
            dates_index, transitions_array,
            date.month, date.day, NEAR_WINDOW)

    return jp_matrix_dict


def compute_delta_jp_matrices(
        gcm_dataset, reference_start_date, reference_end_date,
        simulation_dates_index, lower_bound, upper_bound):
    jp_delta_matrix_lookup = {}

    gcm_reference_period_ds = gcm_dataset.sel(
        time=slice(reference_start_date, reference_end_date))
    transitions_array = state_transition_table(
        gcm_reference_period_ds.pr.values, lower_bound, upper_bound)
    historic_gcm_jp_matrix_lookup = compute_historical_jp_matrices(
        transitions_array, gcm_reference_period_ds.time.values)

    LOGGER.info(
        f'computing GCM delta matrices for period '
        f'{simulation_dates_index.min()} : {simulation_dates_index.max()}')
    date_offset = pandas.DateOffset(days=NEAR_WINDOW)
    for base_date in simulation_dates_index:
        window_start = base_date - date_offset
        window_end = base_date + date_offset
        array = gcm_dataset.sel(time=slice(window_start, window_end)).pr.values
        # TODO: Is it okay for the window to overlap the historical period?
        # This occurs in the first NEAR_WINDOW/2 days of the prediction
        # period.

        jp_matrix = tri_state_joint_probability(
            array, lower_bound, upper_bound)
        gcm_ref_jp_matrix = historic_gcm_jp_matrix_lookup[base_date.month][base_date.day]
        jp_delta_matrix_lookup[base_date] = jp_matrix - gcm_ref_jp_matrix

    return jp_delta_matrix_lookup


def bootstrap_pairs_of_dates(
        observed_data_path, historical_start_date, historical_end_date,
        simulation_dates_index, gcm_jp_delta_matrix_dict):
    dates_lookup = {}
    var = 'pr_Beni01'  # sample data has many precip vars, pick one for now
    historical_dates = date_range_no_leap(
        historical_start_date, historical_end_date)

    obs_df = pandas.read_csv(
        observed_data_path,
        usecols=['Year', 'Month', 'Day', var],
        parse_dates={'date': ['Year', 'Month', 'Day']})

    transitions_array = state_transition_table(
        obs_df[var].values, LOWER_BOUND, UPPER_BOUND)
    historic_obs_jp_matrix_lookup = compute_historical_jp_matrices(
        transitions_array, obs_df.date.values)

    a_date = pandas.Timestamp(numpy.random.choice(historical_dates))

    if obs_df[obs_df['date'] == a_date][var].values[0] <= LOWER_BOUND:
        init_wet_state = DRY
    if obs_df[obs_df['date'] == a_date][var].values[0] > UPPER_BOUND:
        init_wet_state = VERY_WET
    else:
        init_wet_state = WET

    init_sim_dict = {
        'historical_date': a_date,
        'wet_state': init_wet_state,
        'day_of_yr': a_date.dayofyear
    }
    dates_lookup[simulation_dates_index[0]] = init_sim_dict
    print(dates_lookup)


def shift_longitude_from_360(dataset):
    if 'lon' in dataset.dims:
        return dataset.assign_coords(lon=(((dataset.lon + 180) % 360) - 180))
    raise ValueError('Could not reassign longitude coordinates,'
                     'No dimension "lon" in ')


def slice_from_bbox(dataset, minx, miny, maxx, maxy):
    return dataset.isel(
        lon=(dataset.lon > minx) & (dataset.lon < maxx),
        lat=(dataset.lat > miny) & (dataset.lat < maxy))


def date_range_no_leap(start, stop):
    # "noleap" always creates 365-day years. It also forces cftime-format
    # dates, which is not so convenient. isoformat is more interoperable,
    # and pandas DatetimeIndex is most convenient for indexing.
    return pandas.DatetimeIndex(
        [d.isoformat() for d in
            xarray.date_range(start, stop, calendar='noleap')])


if __name__ == "__main__":
    observed_precip_path = os.path.join(
        DATA_STORE_PATH,
        'OBSERVATIONS/LLdM_AOI2/series_pr_diario.csv')
    historical_gcm_path = os.path.join(
        DATA_STORE_PATH, 'GCMs',
        'Amazon__pr_day_CanESM5_historical_r1i1p1f1_gn_18500101-20141231.nc')
    future_gcm_path = os.path.join(
        DATA_STORE_PATH, 'GCMs',
        'Amazon__pr_day_CanESM5_ssp126_r1i1p1f1_gn_20150101-21001231.nc')
    aoi_path = os.path.join(
        DATA_STORE_PATH, 'OBSERVATIONS/LLdM_AOI2/SHP/Basin_LldM.shp')

    bbox_xyxy = pygeoprocessing.get_vector_info(aoi_path)['bounding_box']

    simulation_dates_index = date_range_no_leap(
        PREDICTION_PERIOD_START_DATE,
        PREDICTION_PERIOD_END_DATE)

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
        gcm_jp_delta_matrix_dict = compute_delta_jp_matrices(
            mfds.isel(lon=0, lat=0),  # for now, take first location
            REF_PERIOD_START_DATE,
            REF_PERIOD_END_DATE,
            simulation_dates_index,
            LOWER_BOUND,
            UPPER_BOUND)

    bootstrap_pairs_of_dates(
        observed_precip_path, REF_PERIOD_START_DATE, REF_PERIOD_END_DATE,
        simulation_dates_index, gcm_jp_delta_matrix_dict)
