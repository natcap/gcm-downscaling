from datetime import datetime
import logging
import os
import sys

from knn import knn

args = {
    'reference_period_dates': ('1985-01-01', '2014-12-31'),
    'prediction_dates': ('2030-01-01', '2060-01-01'),
    'hindcast': True,
    'gcm_experiment_list': [knn.GCM_EXPERIMENT_LIST[0]],
    'gcm_model_list': [knn.MODEL_LIST[0]],
    'upper_precip_percentile': 75,
    'lower_precip_threshold': 1,  # millimeter
    'aoi_path': 'H://Shared drives/GCM_Climate_Tool/required_files/OBSERVATIONS/LLdM_AOI2/SHP/Basin_LldM.shp',
    'workspace_dir': 'C://Users/dmf/projects/gcm-project/simple_test',
    'n_workers': -1
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
