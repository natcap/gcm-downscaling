"""
This work is an adaptation of the GCMClimTool Library
by Angarita H., Yates D., Depsky N. 2014-2021
"""

from collections import defaultdict
from datetime import datetime, timedelta
import logging
import os
from pprint import pformat

from affine import Affine
import gcsfs
import geopandas
import google.auth
import numpy
from numpy.lib.stride_tricks import sliding_window_view
import pandas
import pyextremes
import pygeoprocessing
from rasterio.features import rasterize
import scipy
import taskgraph
import xarray

from . import plot

LOGGER = logging.getLogger(__name__)
LOG_FMT = (
    "%(asctime)s "
    "(%(name)s) "
    "%(module)s.%(funcName)s(%(lineno)d) "
    "%(levelname)s %(message)s")

NL = '\n'  # for use in f-strings where \ is not permitted
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
EPSILON = 0.1
KG_M2S_TO_MM = (60 * 60 * 24)  # Kg / (m^2 * s) to mm
EXTREME_PRECIP_QUANTILE = 0.98
MODEL_LIST = [
    'CanESM5',
    'CESM2',
    'CESM2-WACCM',
    'CMCC-CM2-HR4',
    'CMCC-CM2-SR5',
    'CMCC-ESM2',
    'FGOALS-g3',
    'GFDL-EMS4',
    'MIROC6',
    'MPI-ESM1-2-LR',
    # 'CESM-FV2',  # never created zarrs; https://github.com/natcap/gcm-downscaling/issues/19
    # 'IPSL-CM6A-LR', # unreadable by xarray; https://github.com/h5netcdf/h5netcdf/issues/94
    # 'MPI-ESM1-2-HR'  # has a bad file with 0 bytes
]
GCM_EXPERIMENT_LIST = [
    'ssp119',
    'ssp126',
    'ssp245',
    'ssp370',
    'ssp460',
    'ssp585'
]
GCM_PRECIP_VAR = 'pr'
GCM_TEMPERATURE_VAR = 'tas'
GCM_VAR_LIST = [GCM_PRECIP_VAR, GCM_TEMPERATURE_VAR]

GCS_PROTOCOL = 'gcs://'
BUCKET = 'natcap-climate-data'
GCM_PREFIX = 'cmip6'
MSWEP_STORE_PATH = f'{GCS_PROTOCOL}{BUCKET}/mswep_1980_2020.zarr'
MSWEP_DATE_RANGE = ('1980-01-01', '2020-12-31')
MSWEP_VAR = 'precipitation'

# Chunk sizes used to create the zarr stores
# See scripts/preprocessing/*rechunk_to_zarr.py
MSWEP_ZARR_CHUNKS = {
    'lon': 90,
    'lat': 90,
    'time': -1
}
CMIP_ZARR_CHUNKS = {
    'lon': 10,
    'lat': 10,
    'time': -1
}


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
    if 'lon' in dataset.coords:
        return dataset.assign_coords(lon=(((dataset.lon + 180) % 360) - 180))
    raise ValueError(f'Could not reassign longitude coordinates,'
                     f'Expected dimension "lon" but found coordinates {dataset.coords}')


def date_range_no_leap(start, stop):
    """Create a series of dates from start to stop, skipping leap-days.

    calendar='noleap' always creates 365-day years. It also forces
    cftime-format dates.

    Args:
        start (string): of format 'YYYY-MM-DD'
        stop (string): of format 'YYYY-MM-DD'

    Returns:
        (CFTimeIndex) series of cftime objects.
    """
    return xarray.date_range(start, stop, calendar='noleap')


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
    jp_matrix[0, 1] = numpy.count_nonzero(a_low_mask) - jp_matrix[0, 0] - jp_matrix[0, 2]

    a_high_mask = a > upper_bound
    jp_matrix[2, 0] = numpy.count_nonzero(b[a_high_mask] <= lower_bound)
    jp_matrix[2, 2] = numpy.count_nonzero(b[a_high_mask] > upper_bound)
    jp_matrix[2, 1] = numpy.count_nonzero(a_high_mask) - jp_matrix[2, 0] - jp_matrix[2, 2]

    a_med_mask = (a > lower_bound) & (a <= upper_bound)
    jp_matrix[1, 0] = numpy.count_nonzero(b[a_med_mask] <= lower_bound)
    jp_matrix[1, 2] = numpy.count_nonzero(b[a_med_mask] > upper_bound)
    jp_matrix[1, 1] = numpy.count_nonzero(a_med_mask) - jp_matrix[1, 0] - jp_matrix[1, 2]

    return (jp_matrix / numpy.sum(jp_matrix))


def state_transition_series(array, lower_bound, upper_bound):
    def func(a):
        if a[0] <= lower_bound:
            if a[1] <= lower_bound:
                return TRANSITION_TYPES[0, 0]
            elif a[1] > upper_bound:
                return TRANSITION_TYPES[0, 2]
            return TRANSITION_TYPES[0, 1]
        elif a[0] > upper_bound:
            if a[1] <= lower_bound:
                return TRANSITION_TYPES[2, 0]
            elif a[1] > upper_bound:
                return TRANSITION_TYPES[2, 2]
            return TRANSITION_TYPES[2, 1]
        if a[1] <= lower_bound:
            return TRANSITION_TYPES[1, 0]
        elif a[1] > upper_bound:
            return TRANSITION_TYPES[1, 2]
        return TRANSITION_TYPES[1, 1]

    pairs = sliding_window_view(array, window_shape=2)
    transitions = numpy.apply_along_axis(func, 1, pairs)
    return transitions


def slice_dates_around_dayofyear(dates_array, month, day, near_window):
    """Get index of dates from within an n-day window of a given day-of-year.

    For example, get the index of all dates from all years included in
    `dates_array`, within `near_window` days of July 1st.

    Sometimes the `dates_array` comes from a historical record, which
    could have a 366-day year. So rather than getting a window
    around a day-of-year (e.g. day 65), we always get a window
    around a month-day date (e.g. March 3).

    Args:
        dates_array (xarray.DataArray): has a `.dt` accessor; daily frequency
        month (integer): from 1 - 12
        day (integer): day of the month
        near_window (integer): number of days before & after month-day

    Returns:
        (numpy.array): a 1-d array of indexes of size `dates_array`
    """
    idx = numpy.nonzero(
        (dates_array.dt.month == month).data
        & (dates_array.dt.day == day).data)[0]
    # mask to exclude dates with window extending beyond reference period
    mask = (idx - near_window >= 0) & (idx + near_window < dates_array.size)
    ranges = [range(x - near_window, x + near_window + 1) for x in idx[mask]]
    return numpy.array([i for r in ranges for i in r])


def marginal_probability_of_transitions(array):
    """Calculate marginal probability matrix.

    Args:
        array (numpy.array): shape (1, 3)

    Returns:
        (numpy.array): shape (1, 3)
    """
    array[array < 0] = 0
    if array.sum() <= 0:
        # force marginal probabilities to be uniform
        LOGGER.info('all joint probabilities were <= 0')
        array = numpy.ones_like(array)
    return array / array.sum()


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


def bootstrap_dates_precip(
        observed_data_path, prediction_dates, reference_period_dates,
        lower_precip_threshold, upper_precip_percentile, target_csv_path,
        gcm_netcdf_path=None, hindcast=False):
    """Create a sequence of paired dates, matching each date in the
       simulation period with a date from the observed, historical record.

    Args:
        observed_data_path (str): Path to netCDF with dims (`time`) and
            variable ``MSWEP_VAR``, which represents the mean value in
            the area of interest.
        prediction_dates (sequence): ('YYYY-MM-DD', 'YYYY-MM-DD')
            first and last day in the simulation period.
        reference_period_dates (sequence): ('YYYY-MM-DD', 'YYYY-MM-DD')
            first and last day in the reference period, which is used to
            calculate climate "normals".
        lower_precip_threshold (float): the lower boundary of the middle bin
            of precipitation states.
        upper_precip_percentile (float): a percentile (from 0:100) with which
            to extract the absolute precipitation value that will be the upper
            boundary (inclusive) of the middle bin of precipitation states.
        target_csv_path (str): Path to a CSV file created by this function. It
            will contain a sequence of paired dates, matching each date in the
            simulation period with a date from the observed record.
        gcm_netcdf_path (str, optional): Path to a GCM netCDF with dims (`time`)
            and variable ``GCM_PRECIP_VAR``. Required if `hindcast=False`.
        hindcast (bool, optional): If True, prediction dates must be within
            the range of the observed data (i.e. ``MSWEP_DATE_RANGE``) and
            observed data is substituted for GCM data to create
            joint-probability matrices.
    """
    dates_lookup = {}
    gcm_var = GCM_PRECIP_VAR

    with xarray.open_dataset(observed_data_path) as obs_dataset:
        LOGGER.info(
            f'computing observed JP matrices for reference period '
            f'{reference_period_dates[0]} : {reference_period_dates[1]}')
        obs_dataset = obs_dataset.sortby('time')
        obs_reference_period_ds = obs_dataset.sel(
                time=slice(*reference_period_dates))
        lower_bound_obs = lower_precip_threshold
        upper_bound_obs = numpy.percentile(
            obs_reference_period_ds[MSWEP_VAR].values,
            q=[upper_precip_percentile])
        LOGGER.info(
            f'Threshold precip values used for observational record: '
            f'Lower bound: {lower_bound_obs}, '
            f'Upper Bound: {upper_bound_obs}')
        transitions_array_obs = state_transition_series(
            obs_reference_period_ds[MSWEP_VAR].values,
            lower_bound_obs, upper_bound_obs)
        historic_obs_jp_matrix_lookup = compute_historical_jp_matrices(
            transitions_array_obs, obs_reference_period_ds.time)

    # If hindcast, mock the GCM with the observational record
    if hindcast:
        gcm_dataset = obs_dataset
        lower_bound_gcm = lower_bound_obs
        upper_bound_gcm = upper_bound_obs
        historic_gcm_jp_matrix_lookup = historic_obs_jp_matrix_lookup
        gcm_var = MSWEP_VAR
    else:
        with xarray.open_dataset(
            gcm_netcdf_path, decode_times=xarray.coders.CFDatetimeCoder(
                use_cftime=True)) as gcm_dataset:
            gcm_dataset = gcm_dataset.sortby('time')
            LOGGER.info(
                f'computing GCM JP matrices for reference period '
                f'{reference_period_dates[0]} : {reference_period_dates[1]}')
            gcm_reference_period_ds = gcm_dataset.sel(
                time=slice(*reference_period_dates))
            lower_bound_gcm = lower_precip_threshold / KG_M2S_TO_MM
            upper_bound_gcm = numpy.percentile(
                gcm_reference_period_ds[gcm_var].values, q=[upper_precip_percentile])
            LOGGER.info(
                f'Threshold precip values used for GCM record: '
                f'Lower bound: {lower_bound_gcm}, '
                f'Upper Bound: {upper_bound_gcm}')
            transitions_array_gcm = state_transition_series(
                gcm_reference_period_ds[gcm_var].values,
                lower_bound_gcm, upper_bound_gcm)
            historic_gcm_jp_matrix_lookup = compute_historical_jp_matrices(
                transitions_array_gcm, gcm_reference_period_ds.time)

    simulation_dates_index = date_range_no_leap(*prediction_dates)
    LOGGER.info(
        f'Bootstrapping dates for period {simulation_dates_index.min()} : '
        f'{simulation_dates_index.max()}')
    day_one = simulation_dates_index[0]
    window_idx = slice_dates_around_dayofyear(
            obs_reference_period_ds.time, day_one.month, day_one.day, NEAR_WINDOW)
    a_date = obs_reference_period_ds.time[numpy.random.choice(window_idx)].values
    for sim_date in simulation_dates_index:
        # TODO: don't need to re-determine this here, we know the wetstate
        # because we chose it deliberately. Though having this revealed
        # a bug in our indexing....
        precip = obs_reference_period_ds.sel(time=a_date)[MSWEP_VAR].values
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

        date_offset = timedelta(days=NEAR_WINDOW)
        window_start = sim_date - date_offset
        window_end = sim_date + date_offset
        # While most GCM calendars are 365-day "noleap", some are not.
        # isoformat dates are a universal way to slice into these calendars.
        # And it's okay if the result includes a leap-day for some models
        # but not others. https://github.com/natcap/gcm-downscaling/issues/3
        array = gcm_dataset.sel(
            time=slice(window_start.isoformat(), window_end.isoformat())
            )[gcm_var].values
        # TODO: Is it okay for the window to overlap the historical period?
        # This occurs in the first NEAR_WINDOW days of the prediction
        # period.

        jp_matrix = tri_state_joint_probability(
            array, lower_bound_gcm, upper_bound_gcm)
        gcm_ref_jp_matrix = historic_gcm_jp_matrix_lookup[sim_date.month][sim_date.day]
        delta_jp_doy_yr = jp_matrix - gcm_ref_jp_matrix
        observed_jp_doy = historic_obs_jp_matrix_lookup[sim_date.month][sim_date.day]
        projected_jp_matrix = observed_jp_doy + delta_jp_doy_yr
        jp_array = projected_jp_matrix[current_wet_state]
        # some projected probabilities could be negative after adding
        # the delta matrix; clamp at zero
        jp_array[jp_array < 0] = 0
        if jp_array.sum() <= 0:
            # force marginal probabilities to be uniform
            LOGGER.info(
                f'all joint probabilities of transitions from precip state'
                f'[{current_wet_state}] are <= 0 for sim_date {sim_date}')
            jp_array = numpy.ones_like(jp_array)
        margins_array = jp_array / jp_array.sum()
        next_wet_state = numpy.random.choice(
            [DRY, WET, VERY_WET], p=margins_array)
        sim_dict['next_wet_state'] = next_wet_state

        search_window = NEAR_WINDOW
        valid_mask = numpy.array([False])
        while not valid_mask.any():
            # transitions array is one day shorter than historical record
            # so trim the last day off the historical record before indexing.
            window_idx = slice_dates_around_dayofyear(
                obs_reference_period_ds.time[:-1],
                sim_date.month, sim_date.day, search_window)
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
        neighbors = obs_reference_period_ds.isel(time=window_idx[valid_mask])
        inverse_distances = 1 / (
            numpy.abs(neighbors[MSWEP_VAR].values - precip) + EPSILON)
        weights = inverse_distances / inverse_distances.sum()
        nearest = numpy.random.choice(window_idx[valid_mask], p=weights)
        # we want the day after the nearest-neighbor date. But it's not safe
        # to just add 1 day, we might hit Feb 29th, which might not exist in
        # observed record. So find the index of the `nearest` date from the
        # observed record and get the subsequent record's date. (There will
        # always be a 'next day'; see above where window_idx is defined.)
        a_date = obs_reference_period_ds.time[nearest + 1].values
        sim_dict['next_historic_date'] = a_date
        dates_lookup[sim_date] = sim_dict

    dataframe = pandas.DataFrame.from_dict(dates_lookup, orient='index')
    dataframe.to_csv(target_csv_path)
    LOGGER.info(f'Bootstrapped dates complete. Created file: {target_csv_path}')


def sample_from_gpd(series, threshold, n_samples):
    model = pyextremes.EVA(data=series)
    model.get_extremes(
        method='POT',
        extremes_type='high',
        threshold=threshold
    )
    # model.plot_extremes()
    model.fit_model(distribution='genpareto')
    LOGGER.info(model)
    c = model.distribution.mle_parameters['c']
    scale = model.distribution.mle_parameters['scale']
    loc = model.distribution.fixed_parameters['floc']
    sample = scipy.stats.genpareto.rvs(
        c, loc, scale, size=n_samples)
    return sample


def synthesize_extreme_values(
        historical_gcm_path, reference_period_dates,
        forecast_gcm_path, prediction_period_dates, target_csv_path):
    LOGGER.info(f'extreme values analysis for {forecast_gcm_path}')

    def extract_extreme_series(dataset_path, period_dates):
        with xarray.open_dataset(dataset_path) as dataset:
            # pyextremes expects to operate on a pandas series with a
            # DatetimeIndex some GCM have this type already while others
            # have CFTimeIndex
            if not isinstance(dataset.indexes['time'], pandas.DatetimeIndex):
                dataset['time'] = dataset.indexes['time'].to_datetimeindex()
            series = (dataset.sortby('time')
                             .sel(time=slice(*period_dates))
                             .max(['lon', 'lat'])[GCM_PRECIP_VAR]
                             .to_series() * KG_M2S_TO_MM)
        return series

    historical_series = extract_extreme_series(
        historical_gcm_path, reference_period_dates)
    forecast_series = extract_extreme_series(
        forecast_gcm_path, prediction_period_dates)

    threshold = numpy.quantile(historical_series, q=[EXTREME_PRECIP_QUANTILE])[0]
    LOGGER.info(
        f'extreme precip threshold: {threshold} '
        f'for quantile: {EXTREME_PRECIP_QUANTILE}')
    historical_sample = sample_from_gpd(historical_series, threshold, 10000)
    forecast_sample = sample_from_gpd(forecast_series, threshold, 10000)
    df = pandas.DataFrame({
        'historic_sample': numpy.sort(historical_sample),
        'forecast_sample': numpy.sort(forecast_sample),
    })
    df.to_csv(target_csv_path, index=False)
    return threshold


def downscale_precip(
        bootstrapped_dates_path, gridded_observed_precip,
        aoi_mask_path, target_netcdf_path, extreme_value_samples_path=None,
        extreme_precip_threshold=None):
    """Construct a simulated timeseries of gridded precipitation.

    Reshuffle the timeseries of ``gridded_observed_precip`` into the order
    determined by ``bootstrap_dates_precip``. If not a hindcast, also
    apply a magnitude correction based on extreme values analysis of GCM.

    Args:
        bootstrapped_dates_path (str): Path to CSV containing the series of
            bootstrapped dates from ``bootstrap_dates_precip``.
        gridded_observed_precip (str): Path to netCDF representing the gridded
            observed climate variable. This is the grid that will be downscaled
            onto.
        aoi_mask_path (str): Path to netCDF with a boolean `aoi` variable.
        target_netcdf_path (str): Path to netCDF - the downscaled product.
        extreme_value_samples_path (str, optional): Path to CSV with two
            columns: 'historic_sample', 'forecast_sample'.
            Not used for hindcasts.
        extreme_precip_threshold (float, optional): Value determined by extreme
            values analysis; units match ``gridded_observed_precip`` (mm).
            Not used for hindcasts.
    """
    LOGGER.info('Downscaling...')
    dates = pandas.read_csv(bootstrapped_dates_path)
    dates['date'] = pandas.to_datetime(dates.iloc[:, 0])
    if extreme_value_samples_path:
        extremes = pandas.read_csv(extreme_value_samples_path)
    with xarray.open_dataset(gridded_observed_precip) as dataset:
        with xarray.open_dataset(aoi_mask_path) as aoi_dataset:
            dataset['aoi'] = aoi_dataset.aoi
            dataset = dataset.where(dataset.aoi == 1, drop=True)
        precip_array = numpy.empty(
            (len(dates.index), *dataset.precipitation.shape[-2:]))
        for i, d in enumerate(dates.historic_date):
            array = dataset.precipitation.sel(time=d).to_numpy()

            if extreme_precip_threshold is not None:
                extreme_mask = array >= extreme_precip_threshold
                adjusted_values = []
                for x in array[extreme_mask].flatten():
                    idx = numpy.count_nonzero(extremes.historic_sample < x)
                    try:
                        delta = extremes.forecast_sample[idx] - \
                                extremes.historic_sample[idx]
                    except KeyError:
                        # Not sure this can happen, but if x is larger
                        # than all values sampled from the generalized-pareto
                        # distribution, then idx equals the sample size.
                        delta = extremes.forecast_sample[idx - 1] - \
                                extremes.historic_sample[idx - 1]
                    adjusted_values.append(x + delta)
                array[extreme_mask] = adjusted_values

            precip_array[i] = array
        da = xarray.DataArray(
            precip_array,
            coords={
                'time': dates.date.values,
                'lon': dataset.lon,
                'lat': dataset.lat
            },
            dims=dataset.dims
        )
    target_dataset = xarray.Dataset({'precipitation': da})
    target_dataset.to_netcdf(target_netcdf_path)
    LOGGER.info(f'Downscaling complete. Created file: {target_netcdf_path}')


def rasterize_aoi(aoi_path, netcdf_path, target_filepath, fill=0):
    LOGGER.info(f'rasterizing AOI onto {netcdf_path}')
    with xarray.open_dataset(netcdf_path) as dataset:
        coords = dataset.coords

    # netcdf will always have at least 2 coords in each dimension because
    # the bbox that extracted it was buffered by 1 netcdf pixel
    width = coords['lon'][1] - coords['lon'][0]
    height = coords['lat'][1] - coords['lat'][0]  # should be negative
    translation = Affine.translation(
        coords['lon'][0] - width / 2, coords['lat'][0] + numpy.abs(height) / 2)
    scale = Affine.scale(width, height)
    transform = translation * scale

    aoi_df = geopandas.read_file(aoi_path)
    out_shape = (len(coords['lat']), len(coords['lon']))
    raster = rasterize(
        aoi_df.geometry,
        out_shape=out_shape,
        fill=fill,
        transform=transform,
        dtype=int,
        all_touched=True)

    da = xarray.DataArray(
        raster,
        coords={
            'lon': coords['lon'],
            'lat': coords['lat']
        },
        dims=('lat', 'lon'))
    ds = xarray.Dataset({'aoi': da})
    ds.to_netcdf(target_filepath)


def reduce_netcdf(source_file_list, target_filepath, aoi_netcdf_path):
    LOGGER.info(f'averaging {source_file_list} values within AOI')
    with xarray.open_mfdataset(source_file_list,
                               combine='nested', concat_dim='time',
                               data_vars='minimal', coords='minimal',
                               compat='override') as dataset:
        with xarray.open_dataset(aoi_netcdf_path) as aoi_dataset:
            dataset['aoi'] = aoi_dataset.aoi
        dataset = dataset.where(dataset.aoi == 1, drop=True)
        dataset = dataset.mean(['lat', 'lon'])
        dataset.to_netcdf(target_filepath)


def extract_from_zarr(zarr_path, aoi_path, target_path, open_chunks=-1):
    LOGGER.info(f'extracting data from {zarr_path}')
    # TODO: validate aoi has geographic coords
    minx, miny, maxx, maxy = pygeoprocessing.get_vector_info(
        aoi_path)['bounding_box']
    with xarray.open_dataset(
            zarr_path,
            engine='zarr',
            chunks=open_chunks,
            backend_kwargs={"storage_options": {"token": 'anon'}}
            ) as dataset:
        # TODO: in order to accomodate an AOI that crosses 180 deg longitude,
        # it might be best to shift the AOI coordinates to 0:360, rather than
        # shifting the GCM to -180:180.
        dataset = shift_longitude_from_360(dataset)
        # expand the bbox by one pixel on all sides to accomodate
        # rasterization with all_touched=T
        width = numpy.abs(dataset.coords['lon'][1] - dataset.coords['lon'][0])
        height = numpy.abs(dataset.coords['lat'][1] - dataset.coords['lat'][0])
        dataset = dataset.where(
            (dataset.lon >= minx - width)
            & (dataset.lon <= maxx + width)
            & (dataset.lat >= miny - height)
            & (dataset.lat <= maxy + height), drop=True)

        # if the time datatype is not set as datetime64,
        # the time_bnds encoding may need to be set to match that of time
        if dataset.time.dtype == 'object' and ("time_bnds" in dataset and
                                               not dataset["time_bnds"].encoding):
            dataset['time_bnds'].encoding.update(dataset.time.encoding)

        dataset.to_netcdf(target_path)


def validate(dataset, start_date, end_date):
    """Validate that reference or prediction dates are in dataset's range

    Check that start date occurs before dataset's end date and end date
    occurs after dataset's start date (with offset)"""
    # TODO: Also validate for hindcasts,
    # that the prediction dates are within bounds of observed data.
    gcm_start_date = datetime.fromisoformat(dataset.time.min().item().isoformat())
    gcm_end_date = datetime.fromisoformat(dataset.time.max().item().isoformat())
    date_offset = timedelta(days=NEAR_WINDOW)
    if ((datetime.strptime(start_date, '%Y-%m-%d') + date_offset) > gcm_end_date or
            (datetime.strptime(end_date, '%Y-%m-%d') - date_offset) < gcm_start_date):
        raise ValueError(
            f'the requested time period {start_date} : '
            f'{end_date} is outside the time-range of the gcm'
            f'({gcm_start_date} : {gcm_end_date})')


def execute(args):
    """Create a downscaled precipitation product for an area of interest.

    Args:
        args['aoi_path'] (str): a path to a GDAL polygon vector. Coordinates
            represented by longitude, latitude decimal degrees (e.g. WGS84).
        args['workspace_dir'] (str): a path to the directory where this program
            writes output and other temporary files.
        args['reference_period_dates'] (sequence): ('YYYY-MM-DD', 'YYYY-MM-DD')
            first and last day in the reference period, which is used to
            calculate climate "normals".
        args['lower_precip_threshold'] (float): the lower boundary of the
            middle bin of precipitation states. Units: mm
        args['upper_precip_percentile'] (float): a percentile (from 0:100) with
            which to extract the absolute precipitation value that will be the
            upper boundary (inclusive) of the middle bin of precipitation states.
        args['hindcast'] (bool): If True, observed data (MSWEP) is substituted
            for GCM data and the prediction period is set to match the date
            range of the observed dataset (``MSWEP_DATE_RANGE``).
        args['prediction_dates'] (sequence, optional):
            ('YYYY-MM-DD', 'YYYY-MM-DD') first and last day in the simulation period.
            Required if `hindcast=False`.
        args['gcm_model_list'] (sequence, optional): a sequence of strings
            representing CMIP6 model codes. Each model will be used to generate
            a single downscaled product for each experiment in `gcm_experiment_list`.
            Available models are stored in ``GCM_MODEL_LIST``.
            Required if `hindcast=False`.
        args['gcm_experiment_list'] (sequence, optional): a sequence of strings
            representing CMIP6 SSP experiments. Available experiments are
            stored in ``GCM_EXPERIMENT_LIST``. If a CMIP model does not include
            a given experiment, that experiment will be skipped for that model.
            Required if `hindcast=False`.
        args['observed_dataset_path'] (string, optional): if provided, this
            dataset will be used instead of MSWEP as the source of observed,
            historical preciptation. The dataset should be a netCDF or other
            ``xarray.open_dataset`` readable format. It should contain
            coordinates and variables named & defined as,

                Coordinates:
                * ``lat``  - decimal degrees (-90 : 90)
                * ``lon``  - decimal degrees (-180 : 180) or (0 : 360)
                * ``time`` - daily timesteps in units that can be parsed to
                             ``numpy.datetime64``

                Variables:
                * ``precipitation`` - dimensions: (time, lat, lon)
                                      units: millimeter

        args['n_workers'] (int, optional): The number of worker processes to
            use. If omitted, computation will take place in the current process.
            If a positive number, tasks can be parallelized across this many
            processes, which can be useful if `gcm_model_list` or
            `gcm_experiement_list` contain multiple items.
    """
    LOGGER.info(pformat(args))
    taskgraph_working_dir = os.path.join(args['workspace_dir'], '.taskgraph')
    graph = taskgraph.TaskGraph(
        taskgraph_working_dir, args.get('n_workers', -1))
    intermediate_dir = os.path.join(args['workspace_dir'], 'intermediate')
    if not os.path.exists(intermediate_dir):
        os.mkdir(intermediate_dir)

    # Validate reference dates if using MSWEP data
    if 'observed_dataset_path' not in args or \
            args['observed_dataset_path'] is None:
        min_mswep_date = datetime.strptime(MSWEP_DATE_RANGE[0], '%Y-%m-%d')
        max_mswep_date = datetime.strptime(MSWEP_DATE_RANGE[1], '%Y-%m-%d')
        ref_start = datetime.strptime(args['reference_period_dates'][0], '%Y-%m-%d')
        ref_end = datetime.strptime(args['reference_period_dates'][1], '%Y-%m-%d')
        if (ref_start > max_mswep_date or ref_end < min_mswep_date):
            raise ValueError(
                f'the requested reference time period is outside the '
                f'time-range of MSWEP ({min_mswep_date} : {max_mswep_date})'
            )

    mswep_extract_path = os.path.join(intermediate_dir, 'extracted_mswep.nc')
    aoi_mask_mswep_path = os.path.join(intermediate_dir, 'aoi_mask_mswep.nc')
    mswep_netcdf_path = os.path.join(intermediate_dir, 'mswep_mean.nc')

    rasterize_dependent_task_list = []
    if 'observed_dataset_path' in args and \
            args['observed_dataset_path'] is not None:
        mswep_extract_path = args['observed_dataset_path']
    else:
        extract_mswep_task = graph.add_task(
            func=extract_from_zarr,
            kwargs={
                'zarr_path': MSWEP_STORE_PATH,
                'aoi_path': args['aoi_path'],
                'target_path': mswep_extract_path,
                'open_chunks': MSWEP_ZARR_CHUNKS,
            },
            task_name='Extract MSWEP data by bounding box',
            target_path_list=[mswep_extract_path],
            dependent_task_list=[]
        )
        rasterize_dependent_task_list.append(extract_mswep_task)
        graph.join()

    rasterize_aoi_mswep_task = graph.add_task(
        func=rasterize_aoi,
        kwargs={
            'aoi_path': args['aoi_path'],
            'netcdf_path': mswep_extract_path,
            'target_filepath': aoi_mask_mswep_path,
        },
        task_name='Rasterize AOI onto the MSWEP grid.',
        target_path_list=[aoi_mask_mswep_path],
        dependent_task_list=rasterize_dependent_task_list
    )

    reduce_mswep_task = graph.add_task(
        func=reduce_netcdf,
        kwargs={
            'source_file_list': [mswep_extract_path],
            'target_filepath': mswep_netcdf_path,
            'aoi_netcdf_path': aoi_mask_mswep_path
        },
        task_name='Reduce MSWEP to average value within AOI.',
        target_path_list=[mswep_netcdf_path],
        dependent_task_list=[rasterize_aoi_mswep_task]
    )

    if args['hindcast']:
        hindcast_target_csv_path = os.path.join(
            intermediate_dir, 'bootstrapped_dates_precip_hindcast.csv')
        hindcast_date_range = MSWEP_DATE_RANGE
        if 'observed_dataset_path' in args and \
                args['observed_dataset_path'] is not None:
            with xarray.open_dataset(mswep_netcdf_path) as dataset:
                min_date = str(dataset.time.min().values)[:10]
                max_date = str(dataset.time.max().values)[:10]
            hindcast_date_range = (min_date, max_date)
        hind_bootstrap_dates_task = graph.add_task(
            func=bootstrap_dates_precip,
            kwargs={
                'observed_data_path': mswep_netcdf_path,
                'prediction_dates': hindcast_date_range,
                'reference_period_dates': args['reference_period_dates'],
                'lower_precip_threshold': args['lower_precip_threshold'],
                'upper_precip_percentile': args['upper_precip_percentile'],
                'target_csv_path': hindcast_target_csv_path,
                'hindcast': True
            },
            task_name='Bootstrap dates for precipitation',
            target_path_list=[hindcast_target_csv_path],
            dependent_task_list=[reduce_mswep_task]
        )
        hindcast_target_netcdf_path = os.path.join(
            args['workspace_dir'], 'downscaled_precip_hindcast.nc')
        hind_downscale_precip_task = graph.add_task(
            func=downscale_precip,
            kwargs={
                'bootstrapped_dates_path': hindcast_target_csv_path,
                'gridded_observed_precip': mswep_extract_path,
                'aoi_mask_path': aoi_mask_mswep_path,
                'target_netcdf_path': hindcast_target_netcdf_path
            },
            task_name='Downscale Precipitation',
            target_path_list=[hindcast_target_netcdf_path],
            dependent_task_list=[hind_bootstrap_dates_task]
        )
        hindcast_target_pdf_path = os.path.splitext(
            hindcast_target_netcdf_path)[0] + '.pdf'
        report_task = graph.add_task(
            func=plot.plot,
            kwargs={
                'dates_filepath': hindcast_target_csv_path,
                'precip_filepath': hindcast_target_netcdf_path,
                'observed_mean_precip_filepath': mswep_netcdf_path,
                'observed_precip_filepath': mswep_extract_path,
                'aoi_netcdf_path': aoi_mask_mswep_path,
                'reference_period_dates': args['reference_period_dates'],
                'hindcast': True,
                'target_filename': hindcast_target_pdf_path
            },
            task_name='Report',
            target_path_list=[hindcast_target_pdf_path],
            dependent_task_list=[hind_downscale_precip_task]
        )

    try:
        gcm_model_list = args['gcm_model_list']
    except KeyError:
        gcm_model_list = []
    gcs_filesystem = gcsfs.GCSFileSystem(token='anon')
    for gcm_model in gcm_model_list:
        historical_gcm_files = gcs_filesystem.glob(
            f"{BUCKET}/{GCM_PREFIX}/{gcm_model}/{GCM_PRECIP_VAR}_day_{gcm_model}_historical_*.zarr")
        if len(historical_gcm_files) == 0:
            LOGGER.warning(
                f'No files found for model: {gcm_model}, experiment: historical'
                f'skipping model {gcm_model}.')
            continue
        if len(historical_gcm_files) > 1:
            LOGGER.warning(
                f'Ambiguous files found for model: {gcm_model}, experiment: historical'
                f'Found: {historical_gcm_files}'
                f'skipping model {gcm_model}.')
            continue

        # validate that reference dates fall within range of historical data
        with xarray.open_dataset(
                f'{GCS_PROTOCOL}{historical_gcm_files[0]}',
                decode_times=xarray.coders.CFDatetimeCoder(use_cftime=True),
                engine='zarr',
                backend_kwargs={"storage_options": {"token": 'anon'}}
                    ) as gcm_hist_dataset:
            validate(gcm_hist_dataset, *args['reference_period_dates'])
        # validate forecast dates fall within range of future data
        future_gcm_files = gcs_filesystem.glob(
                f"{BUCKET}/{GCM_PREFIX}/{gcm_model}/{GCM_PRECIP_VAR}_day_{gcm_model}_ssp*.zarr")
        if len(future_gcm_files) == 0:
            LOGGER.warning(
                f'No files found for model: {gcm_model}. Skipping.')
            continue
        with xarray.open_dataset(
                f'{GCS_PROTOCOL}{future_gcm_files[0]}',
                decode_times=xarray.coders.CFDatetimeCoder(use_cftime=True),
                engine='zarr',
                backend_kwargs={"storage_options": {"token": 'anon'}}
                ) as future_gcm_dataset:
            validate(future_gcm_dataset, *args['prediction_dates'])

        gcm_historical_extract_path = os.path.join(
            intermediate_dir, f'extracted_{gcm_model}_historical.nc')
        extract_historical_gcm_task = graph.add_task(
            func=extract_from_zarr,
            kwargs={
                'zarr_path': f'{GCS_PROTOCOL}{historical_gcm_files[0]}',
                'aoi_path': args['aoi_path'],
                'target_path': gcm_historical_extract_path,
                'open_chunks': CMIP_ZARR_CHUNKS
            },
            task_name='Extract GCM historical data by bounding box',
            target_path_list=[gcm_historical_extract_path],
            dependent_task_list=[]
        )

        aoi_mask_gcm_path = os.path.join(
            intermediate_dir, f'aoi_mask_{gcm_model}.nc')
        rasterize_aoi_gcm_task = graph.add_task(
            func=rasterize_aoi,
            kwargs={
                'aoi_path': args['aoi_path'],
                'netcdf_path': gcm_historical_extract_path,
                'target_filepath': aoi_mask_gcm_path,
            },
            task_name='Rasterize AOI onto the GCM grid.',
            target_path_list=[aoi_mask_gcm_path],
            dependent_task_list=[extract_historical_gcm_task]
        )
        for gcm_experiment in args['gcm_experiment_list']:
            future_gcm_files = gcs_filesystem.glob(
                f"{BUCKET}/{GCM_PREFIX}/{gcm_model}/{GCM_PRECIP_VAR}_day_{gcm_model}_{gcm_experiment}_*.zarr")

            if len(future_gcm_files) == 0:
                LOGGER.warning(
                    f'No files found for model: {gcm_model}, experiment: {gcm_experiment}'
                    f'skipping experment: {gcm_experiment} - {gcm_model}.')
                continue
            if len(future_gcm_files) > 1:
                LOGGER.warning(
                    f'Ambiguous files found for model: {gcm_model}, experiment: {gcm_experiment}'
                    f'Found: {future_gcm_files}'
                    f'skipping experiment: {gcm_experiment} - {gcm_model}.')
                continue
            LOGGER.info(f'Starting {gcm_model} {gcm_experiment}')

            target_csv_path = os.path.join(
                intermediate_dir,
                f'bootstrapped_dates_precip_{gcm_model}_{gcm_experiment}.csv')
            target_netcdf_path = os.path.join(
                args['workspace_dir'],
                f'downscaled_precip_{gcm_model}_{gcm_experiment}.nc')

            gcm_netcdf_path = os.path.join(
                intermediate_dir,
                f"{GCM_PRECIP_VAR}_day_{gcm_model}_{gcm_experiment}_mean.nc")
            gcm_future_extract_path = os.path.join(
                intermediate_dir, f'extracted_{gcm_model}_{gcm_experiment}.nc')
            extract_future_gcm_task = graph.add_task(
                func=extract_from_zarr,
                kwargs={
                    'zarr_path': f'{GCS_PROTOCOL}{future_gcm_files[0]}',
                    'aoi_path': args['aoi_path'],
                    'target_path': gcm_future_extract_path,
                    'open_chunks': CMIP_ZARR_CHUNKS
                },
                task_name='Extract GCM future data by bounding box',
                target_path_list=[gcm_future_extract_path],
                dependent_task_list=[]
            )

            target_extreme_values_path = os.path.join(
                intermediate_dir,
                f'sythesized_extreme_precip_{gcm_model}_{gcm_experiment}.csv')
            extreme_values_task = graph.add_task(
                func=synthesize_extreme_values,
                kwargs={
                    'historical_gcm_path': gcm_historical_extract_path,
                    'reference_period_dates': args['reference_period_dates'],
                    'forecast_gcm_path': gcm_future_extract_path,
                    'prediction_period_dates': args['prediction_dates'],
                    'target_csv_path': target_extreme_values_path,
                },
                task_name='Synthesize extreme values',
                target_path_list=[target_extreme_values_path],
                store_result=True,
                dependent_task_list=[extract_future_gcm_task,
                                     extract_historical_gcm_task]
            )

            reduce_gcm_task = graph.add_task(
                func=reduce_netcdf,
                kwargs={
                    'source_file_list': [
                        gcm_historical_extract_path, gcm_future_extract_path],
                    'aoi_netcdf_path': aoi_mask_gcm_path,
                    'target_filepath': gcm_netcdf_path,
                },
                task_name='Reduce GCM to average value within AOI.',
                target_path_list=[gcm_netcdf_path],
                dependent_task_list=[
                    rasterize_aoi_gcm_task,
                    extract_historical_gcm_task,
                    extract_future_gcm_task]
            )

            bootstrap_dates_task = graph.add_task(
                func=bootstrap_dates_precip,
                kwargs={
                    'observed_data_path': mswep_netcdf_path,
                    'prediction_dates': args['prediction_dates'],
                    'reference_period_dates': args['reference_period_dates'],
                    'gcm_netcdf_path': gcm_netcdf_path,
                    'lower_precip_threshold': args['lower_precip_threshold'],
                    'upper_precip_percentile': args['upper_precip_percentile'],
                    'target_csv_path': target_csv_path
                },
                task_name='Bootstrap dates for precipitation',
                target_path_list=[target_csv_path],
                dependent_task_list=[reduce_gcm_task, reduce_mswep_task]
            )

            # TODO: is it problematic for these task objects created in for-loop
            # to overwrite each other?
            extreme_precip_threshold = extreme_values_task.get()
            downscale_precip_task = graph.add_task(
                func=downscale_precip,
                kwargs={
                    'bootstrapped_dates_path': target_csv_path,
                    'gridded_observed_precip': mswep_extract_path,
                    'aoi_mask_path': aoi_mask_mswep_path,
                    'target_netcdf_path': target_netcdf_path,
                    'extreme_value_samples_path': target_extreme_values_path,
                    'extreme_precip_threshold': extreme_precip_threshold
                },
                task_name='Downscale Precipitation',
                target_path_list=[target_netcdf_path],
                dependent_task_list=[bootstrap_dates_task, extreme_values_task]
            )

            target_pdf_path = os.path.splitext(target_netcdf_path)[0] + '.pdf'
            report_task = graph.add_task(
                func=plot.plot,
                kwargs={
                    'dates_filepath': target_csv_path,
                    'precip_filepath': target_netcdf_path,
                    'observed_mean_precip_filepath': mswep_netcdf_path,
                    'observed_precip_filepath': mswep_extract_path,
                    'aoi_netcdf_path': aoi_mask_mswep_path,
                    'reference_period_dates': args['reference_period_dates'],
                    'hindcast': False,
                    'target_filename': target_pdf_path
                },
                task_name='Report',
                target_path_list=[target_pdf_path],
                dependent_task_list=[downscale_precip_task]
            )

    graph.join()
