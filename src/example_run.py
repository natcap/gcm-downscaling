from datetime import datetime
import logging
import os
import sys

import knn


data_store_path = 'H://Shared drives/GCM_Climate_Tool/required_files'

# TODO: If we expect that any of these arguments will *never* change
# we should define them in knn.py instead of here.
args = {
    'ref_period_start_date': '1985-01-01',
    'ref_period_end_date': '2014-12-31',
    'prediction_start_date': '2030-01-01',
    'prediction_end_date': '2060-01-01',
    'hindcast': False,
    'data_store_path': data_store_path,
    'mswep_store_path': 'C://Users/dmf/projects/gcm-project/mswep_annual',
    'gcm_var': 'pr',
    'gcm_experiment_list': knn.GCM_EXPERIMENT_LIST,
    'gcm_model_list': ['CanESM5'],
    'upper_precip_percentile': (75),
    'lower_precip_threshold': 1,  # millimeter
    'aoi_path': os.path.join(
        data_store_path, 'OBSERVATIONS/LLdM_AOI2/SHP/Basin_LldM.shp'),
    # 'observed_precip_path': os.path.join(
    #     data_store_path, 'OBSERVATIONS/LLdM_AOI2/series_pr_diario_regional_average.nc'),
    'workspace_dir': 'C://Users/dmf/projects/gcm-project/output_mswep'
}

if __name__ == '__main__':
    if not os.path.exists(args['workspace_dir']):
        os.mkdir(args['workspace_dir'])
    logfile = os.path.join(
        args['workspace_dir'],
        f'log_{datetime.now().strftime("%Y-%m-%d--%H_%M_%S")}.txt')
    formatter = logging.Formatter(knn.LOG_FMT)
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(logfile)
    stream_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)
    logging.basicConfig(
        level=logging.INFO,
        handlers=[stream_handler, file_handler])

    knn.execute(args)
