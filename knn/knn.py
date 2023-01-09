from collections import defaultdict
from datetime import datetime, timedelta
import glob
import logging
import os
import re

from affine import Affine
import geopandas
import numpy
from numpy.lib.stride_tricks import sliding_window_view
import pandas
import rasterio
import taskgraph
import xarray

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
MODEL_LIST = [
    'CanESM5',
    'CESM2',
    'CESM2-WACCM',
    'CESM-FV2',
    'CMCC-CM2-HR4',
    'CMCC-CM2-SR5',
    'CMCC-CM2-ESM2',
    'FGOALS-g3',
    'GFDL-EMS4',
    'IPSL-CM6A-LR',
    'MIROC6',
    'MPI_ESM1-2-LR',
    # 'MPI_ESM1-2-HR'  # has a bad file with 0 bytes
]
GCM_EXPERIMENT_LIST = [
    'historical',
    'ssp119',
    'ssp126',
    'ssp245',
    'ssp370',
    'ssp460',
    'ssp585'
]
GCM_VAR_LIST = ['pr', 'tas']


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


def marginal_probability_of_transitions(projected_jp_matrix):
    """Add joint-probability matrices; calculate marginal probability matrix.

    Args:
        observed_matrix (numpy.array): 2-d array of shape (3,3)
        delta_matrix (numpy.array): 2-d array of shape (3,3)

    Returns:
        (numpy.array): 2-d array of shape (3,3)
    """
    # projected_jp_matrix = observed_matrix + delta_matrix
    projected_jp_matrix[projected_jp_matrix < 0] = 0
    # TODO: In R, it appears there's another normalization step here:
    # projected_jp_matrix / projected_jp_matrix.sum()

    def func(x):
        if x.sum() <= 0:
            # all joint probabilities <=0;
            # force marginal probabilities to be uniform
            # TODO: is there a better way to handle this?
            x = numpy.ones_like(x)
            LOGGER.info('all joint probabilities were <= 0')
        return x / x.sum()

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
        observed_data_path, prediction_start_date, prediction_end_date,
        reference_start_date, reference_end_date, gcm_netcdf_path,
        lower_precip_threshold, upper_precip_percentile, target_csv_path,
        hindcast=False):
    dates_lookup = {}

    with xarray.open_dataset(observed_data_path) as obs_dataset:
        LOGGER.info(
            f'computing observed JP matrices for reference period '
            f'{reference_start_date} : {reference_end_date}')
        obs_reference_period_ds = obs_dataset.sel(
                time=slice(reference_start_date, reference_end_date))
        lower_bound_obs = lower_precip_threshold
        upper_bound_obs = numpy.percentile(
            obs_reference_period_ds.pr.values, q=upper_precip_percentile)
        LOGGER.info(
            f'Threshold precip values used for observational record: '
            f'Lower bound: {lower_bound_obs}, '
            f'Upper Bound: {upper_bound_obs}')
        transitions_array_obs = state_transition_series(
            obs_reference_period_ds.pr.values,
            lower_bound_obs, upper_bound_obs)
        historic_obs_jp_matrix_lookup = compute_historical_jp_matrices(
            transitions_array_obs, obs_reference_period_ds.time)

    # If hindcast, mock the GCM with the observational record
    if hindcast:
        gcm_dataset = obs_dataset
        lower_bound_gcm = lower_bound_obs
        upper_bound_gcm = upper_bound_obs
        historic_gcm_jp_matrix_lookup = historic_obs_jp_matrix_lookup
    else:
        with xarray.open_dataset(gcm_netcdf_path, use_cftime=True) as gcm_dataset:
            validate(gcm_dataset, prediction_start_date, prediction_end_date)
            LOGGER.info(
                f'computing GCM JP matrices for reference period '
                f'{reference_start_date} : {reference_end_date}')
            gcm_reference_period_ds = gcm_dataset.sel(
                time=slice(reference_start_date, reference_end_date))
            lower_bound_gcm = lower_precip_threshold / KG_M2S_TO_MM
            upper_bound_gcm = numpy.percentile(
                gcm_reference_period_ds.pr.values, q=upper_precip_percentile)
            LOGGER.info(
                f'Threshold precip values used for GCM record: '
                f'Lower bound: {lower_bound_gcm}, '
                f'Upper Bound: {upper_bound_gcm}')
            transitions_array_gcm = state_transition_series(
                gcm_reference_period_ds.pr.values,
                lower_bound_gcm, upper_bound_gcm)
            historic_gcm_jp_matrix_lookup = compute_historical_jp_matrices(
                transitions_array_gcm, gcm_reference_period_ds.time)

    simulation_dates_index = date_range_no_leap(
        prediction_start_date, prediction_end_date)
    LOGGER.info(
        f'Simulating precip for period {simulation_dates_index.min()} : '
        f'{simulation_dates_index.max()}')
    day_one = simulation_dates_index[0]
    window_idx = slice_dates_around_dayofyear(
            obs_reference_period_ds.time, day_one.month, day_one.day, NEAR_WINDOW)
    a_date = obs_reference_period_ds.time[numpy.random.choice(window_idx)].values
    for sim_date in simulation_dates_index:
        # TODO: don't need to re-determine this here, we know the wetstate
        # because we chose it deliberately. Though having this revealed
        # a bug in our indexing....
        precip = obs_reference_period_ds.sel(time=a_date).pr.values
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
        array = gcm_dataset.sel(time=slice(window_start, window_end)).pr.values
        # TODO: Is it okay for the window to overlap the historical period?
        # This occurs in the first NEAR_WINDOW days of the prediction
        # period.

        jp_matrix = tri_state_joint_probability(
            array, lower_bound_gcm, upper_bound_gcm)
        gcm_ref_jp_matrix = historic_gcm_jp_matrix_lookup[sim_date.month][sim_date.day]
        delta_jp_doy_yr = jp_matrix - gcm_ref_jp_matrix
        observed_jp_doy = historic_obs_jp_matrix_lookup[sim_date.month][sim_date.day]
        projected_jp_matrix = observed_jp_doy + delta_jp_doy_yr
        margins_matrix = marginal_probability_of_transitions(
            projected_jp_matrix)
        next_wet_state = numpy.random.choice(
            [DRY, WET, VERY_WET], p=margins_matrix[current_wet_state])
        sim_dict['next_wet_state'] = next_wet_state

        valid_mask = numpy.array([False])
        search_window = NEAR_WINDOW
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
            numpy.abs(neighbors.pr.values - precip) + EPSILON)
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
    LOGGER.info(f'Simulation complete. Created file: {target_csv_path}')


def rasterize_aoi(aoi_path, netcdf_path, target_filepath, fill=0):
    LOGGER.info('rasterizing AOI')
    # TODO: in order to accomodate an AOI that crosses 180 deg longitude,
    # it might be best to shift the AOI coordinates to 0:360, rather than
    # shifting the GCM to -180:180.
    with xarray.open_dataset(netcdf_path) as dataset:
        dataset = shift_longitude_from_360(dataset)
        coords = dataset.coords

    width = coords['lon'][1] - coords['lon'][0]
    height = coords['lat'][1] - coords['lat'][0]  # should be negative
    translation = Affine.translation(
        coords['lon'][0] - width / 2, coords['lat'][0] + numpy.abs(height) / 2)
    scale = Affine.scale(width, height)
    transform = translation * scale

    aoi_df = geopandas.read_file(aoi_path)
    out_shape = (len(coords['lat']), len(coords['lon']))
    raster = rasterio.features.rasterize(
        aoi_df.geometry,
        out_shape=out_shape,
        fill=fill,
        transform=transform,
        dtype=int)

    da = xarray.DataArray(
        raster,
        coords={
            'lon': coords['lon'],
            'lat': coords['lat']
        },
        dims=('lat', 'lon'))
    ds = xarray.Dataset({'aoi': da})
    ds.to_netcdf(target_filepath)


def reduce_netcdf(source_files_list, aoi_netcdf_path, target_filepath):
    LOGGER.info('averaging GCM values within AOI')
    # TODO: validate that AOI has geographic coordinates
    # TODO: try smaller lat/lon chunks
    with xarray.open_mfdataset(source_files_list, parallel=True,
                               combine='nested', concat_dim='time',
                               data_vars='minimal', coords='minimal',
                               compat='override') as dataset:
        # TODO: this could be optimized by slicing out a date range
        # first. But since we always need to include the reference period
        # and the simulation period, it's convenient to just keep the whole
        # timeseries for now.
        with xarray.open_dataset(aoi_netcdf_path) as aoi_dataset:
            dataset['aoi'] = aoi_dataset.aoi
            dataset = dataset.where(dataset.aoi == 1, drop=True)
            # dataset = dataset.sel(aoi=1)
            dataset = dataset.mean(['lat', 'lon'])
            dataset.to_netcdf(target_filepath)


def validate(dataset, prediction_start_date, prediction_end_date):
    # TODO: Also validate for hindcasts,
    # that the prediction dates are within bounds of observed data.
    gcm_start_date = datetime.fromisoformat(dataset.time.min().item().isoformat())
    gcm_end_date = datetime.fromisoformat(dataset.time.max().item().isoformat())
    date_offset = timedelta(days=NEAR_WINDOW)
    if ((datetime.strptime(prediction_end_date, '%Y-%m-%d') + date_offset) > gcm_end_date or
            (datetime.strptime(prediction_start_date, '%Y-%m-%d') - date_offset) < gcm_start_date):
        raise ValueError(
            f'the requested prediction period {prediction_start_date} : '
            f'{prediction_end_date} is outside the time-range of the gcm'
            f'({gcm_start_date} : {gcm_end_date})')


def execute(args):
    """
    Args:
        args[data_store_path] (string): path to store of CMIP netcdf files,
            with subdirectories for each model (e.g. "CESM2"). Each NetCDF
            should have coordinates named 'lon', 'lat', 'time'. Filenames
            should follow this pattern:
                f"{args['gcm_var']}_day_{gcm_model}_{gcm_experiment}_*.nc"

    """
    LOGGER.info(args)
    taskgraph_working_dir = os.path.join(args['workspace_dir'], '.taskgraph')
    graph = taskgraph.TaskGraph(taskgraph_working_dir, -1)

    # MSWEP Files are named as '%Y%j.nc'; https://strftime.org/
    mswep_store = args['mswep_store_path']
    # mswep_store = os.path.join(
    #     args['data_store_path'], 'OBSERVATIONS', 'Global MSWEP2', 'Daily')
    # All files from years in the reference period. Don't worry about excluding
    # extra days if the ref period starts/ends in the middle of a year.
    mswep_files = [
        os.path.join(mswep_store, f) for f in os.listdir(mswep_store)
        if re.search(r'[0-9]{7}\.nc', f)
        and (int(f[:4]) >= int(args['ref_period_start_date'][:4]))
        and (int(f[:4]) <= int(args['ref_period_end_date'][:4]))]

    temp_aoi_mswep_path = os.path.join(
            args['workspace_dir'], 'aoi_mswep.nc')
    temp_mswep_netcdf_path = os.path.join(
            args['workspace_dir'], 'mswep_mean.nc')

    # Only the geotransform is needed from the netcdf, which is
    # the same for all mswep files, so always use the same file for
    # taskgraph benefits.
    representative_mswep_path = os.path.join(mswep_store, '1980.nc')
    rasterize_aoi_mswep_task = graph.add_task(
        func=rasterize_aoi,
        kwargs={
            'aoi_path': args['aoi_path'],
            'netcdf_path': representative_mswep_path,
            'target_filepath': temp_aoi_mswep_path,
        },
        task_name='Rasterize AOI onto the GCM grid.',
        target_path_list=[temp_aoi_mswep_path],
        dependent_task_list=[]
    )
    reduce_mswep_task = graph.add_task(
        func=reduce_netcdf,
        kwargs={
            'source_files_list': mswep_files,
            'aoi_netcdf_path': temp_aoi_mswep_path,
            'target_filepath': temp_mswep_netcdf_path,
        },
        task_name='Reduce GCM to average value within AOI.',
        target_path_list=[temp_mswep_netcdf_path],
        dependent_task_list=[rasterize_aoi_mswep_task]
    )

    if args['hindcast']:
        target_csv_path = os.path.join(
            args['workspace_dir'], 'downscaled_precip_hindcast.csv')
        temp_netcdf_path = None
        downscale_precipitation(
                args['observed_precip_path'],
                args['prediction_start_date'],
                args['prediction_end_date'],
                args['ref_period_start_date'],
                args['ref_period_end_date'],
                temp_netcdf_path,
                args['lower_precip_threshold'],
                args['upper_precip_percentile'],
                target_csv_path,
                hindcast=True)
        return None

    for gcm_model in args['gcm_model_list']:
        historical_gcm_files = glob.glob(
            os.path.join(args['data_store_path'], 'GCMs', 'CMIP6', gcm_model,
                         f"{args['gcm_var']}_day_{gcm_model}_historical_*.nc"))
        if len(historical_gcm_files) == 0:
            LOGGER.info(
                f'No files found for model: {gcm_model}, experiment: historical'
                f'skipping model {gcm_model}.')
            continue
        temp_aoi_path = os.path.join(
            args['workspace_dir'], f'aoi_{gcm_model}.nc')
        rasterize_aoi_gcm_task = graph.add_task(
            func=rasterize_aoi,
            kwargs={
                'aoi_path': args['aoi_path'],
                'netcdf_path': historical_gcm_files[0],
                'target_filepath': temp_aoi_path,
            },
            task_name='Rasterize AOI onto the GCM grid.',
            target_path_list=[temp_aoi_path],
            dependent_task_list=[]
        )
        for gcm_experiment in args['gcm_experiment_list']:
            future_gcm_files = glob.glob(
                os.path.join(args['data_store_path'], 'GCMs', 'CMIP6', gcm_model,
                             f"{args['gcm_var']}_day_{gcm_model}_{gcm_experiment}_*.nc"))

            if len(future_gcm_files) == 0:
                LOGGER.info(
                    f'No files found for model: {gcm_model}, experiment: {gcm_experiment}')
                LOGGER.info(f'skipping experment {gcm_experiment}.')
                continue
            LOGGER.info(f'Starting {gcm_model} {gcm_experiment}')

            target_csv_path = os.path.join(
                args['workspace_dir'],
                f'downscaled_precip_{gcm_model}_{gcm_experiment}.csv')
            temp_netcdf_path = os.path.join(args['workspace_dir'],
                f"{args['gcm_var']}_day_{gcm_model}_{gcm_experiment}_mean_in_aoi.nc")

            reduce_gcm_task = graph.add_task(
                func=reduce_netcdf,
                kwargs={
                    'source_files_list': [*historical_gcm_files, *future_gcm_files],
                    'aoi_netcdf_path': temp_aoi_path,
                    'target_filepath': temp_netcdf_path,
                },
                task_name='Reduce GCM to average value within AOI.',
                target_path_list=[temp_netcdf_path],
                dependent_task_list=[rasterize_aoi_gcm_task]
            )

            _ = graph.add_task(
                func=downscale_precipitation,
                kwargs={
                    'observed_data_path': temp_mswep_netcdf_path,
                    'prediction_start_date': args['prediction_start_date'],
                    'prediction_end_date': args['prediction_end_date'],
                    'reference_start_date': args['ref_period_start_date'],
                    'reference_end_date': args['ref_period_end_date'],
                    'gcm_netcdf_path': temp_netcdf_path,
                    'lower_precip_threshold': args['lower_precip_threshold'],
                    'upper_precip_percentile': args['upper_precip_percentile'],
                    'target_csv_path': target_csv_path
                },
                task_name='Downscale Precipitation',
                target_path_list=[target_csv_path],
                dependent_task_list=[reduce_gcm_task, reduce_mswep_task]
            )

    graph.join()
