import os

import pandas
import xarray

from .joint_probability import tri_state_joint_probability

LOWER_BOUND = 2.0e-06  # in units of the GCM pr variable
UPPER_BOUND = 4.0e-06
DAYS_IN_YEAR = 365
WINDOW_LENGTH = (DAYS_IN_YEAR - 1) / 2
REF_PERIOD_START_DATE = '1985-01-01'
REF_PERIOD_END_DATE = '2014-12-31'
PREDICTION_PERIOD_START_DATE = '2015-01-01'
PREDICTION_PERIOD_END_DATE = '2015-01-31'
DATA_STORE_PATH = 'H://Shared drives/GCM_Climate_Tool/required_files/GCMs'


def compute_delta_jp_matrices(
        gcm_dataset, reference_start_date, reference_end_date,
        prediction_start_date, prediction_end_date,
        lower_bound, upper_bound):

    jp_delta_matrix_lookup = {}

    precip_array = gcm_dataset.sel(time=slice(
        reference_start_date, reference_end_date)).pr.values  # pr is variable name
    gcm_ref_jp_matrix = tri_state_joint_probability(
        precip_array,
        lower_bound,
        upper_bound)
    print(gcm_ref_jp_matrix)
    print(gcm_ref_jp_matrix.sum())

    prediction_period_dates = pandas.date_range(
        prediction_start_date, prediction_end_date)

    for base_date in prediction_period_dates:
        # TODO: Is it okay for the window to overlap the historical period?
        # This occurs in the first WINDOW_LENGTH/2 days of the prediction
        # period.
        # Or should we not start predicting until July of the first year?
        date_offset = pandas.DateOffset(days=WINDOW_LENGTH)
        window_start = base_date - date_offset
        window_end = base_date + date_offset
        print(window_start)
        print(window_end)

        array = gcm_dataset.sel(time=slice(
            window_start, window_end)).pr.values

        jp_matrix = tri_state_joint_probability(
            array,
            lower_bound, upper_bound)

        # TODO: rescale this delta matrix from 0 - 1.
        jp_delta_matrix_lookup[base_date] = jp_matrix - gcm_ref_jp_matrix

    return jp_delta_matrix_lookup


if __name__ == "__main__":
    historical_ds_path = os.path.join(
        DATA_STORE_PATH,
        'Amazon__pr_day_CanESM5_historical_r1i1p1f1_gn_18500101-20141231.nc')
    future_ds_path = os.path.join(
        DATA_STORE_PATH,
        'Amazon__pr_day_CanESM5_ssp126_r1i1p1f1_gn_20150101-21001231.nc')

    with xarray.open_mfdataset([historical_ds_path, future_ds_path]) as mfds:
        jp_delta_matrix_lookup = compute_delta_jp_matrices(
            mfds.isel(lon=0, lat=0),  # for now, take first location
            REF_PERIOD_START_DATE, REF_PERIOD_END_DATE,
            PREDICTION_PERIOD_START_DATE, PREDICTION_PERIOD_END_DATE,
            LOWER_BOUND, UPPER_BOUND)

    print(jp_delta_matrix_lookup)



