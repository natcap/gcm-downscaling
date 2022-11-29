import logging
import os
import sys

import knn


data_store_path = 'H://Shared drives/GCM_Climate_Tool/required_files'

args = {
    'ref_period_start_date': '1985-01-01',
    'ref_period_end_date': '2014-12-31',
    'prediction_start_date': '2030-01-01',
    'prediction_end_date': '2050-01-01',
    'hindcast': False,
    'data_store_path': data_store_path,
    'gcm_var': 'pr',
    'gcm_experiment_list': ['ssp126'],
    'gcm_model_list': ['CanESM5'],
    'upper_precip_percentile': (75),
    'lower_precip_threshold': 1,  # millimeter
    'aoi_path': os.path.join(
        data_store_path, 'OBSERVATIONS/LLdM_AOI2/SHP/Basin_LldM.shp'),
    'observed_precip_path': os.path.join(
        data_store_path, 'OBSERVATIONS/LLdM_AOI2/series_pr_diario_regional_average.nc'),
    'workspace_dir': 'C://Users/dmf/projects/gcm-project/output'
}

if __name__ == '__main__':
    if not os.path.exists(args['workspace_dir']):
        os.mkdir(args['workspace_dir'])
    logging.basicConfig(
        level=logging.INFO,
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(
                os.path.join(args['workspace_dir'], 'log.txt'))])

    knn.execute(args)
