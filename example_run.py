from datetime import datetime
import logging
import os
import sys

from knn import knn

args = {
    'aoi_path': '../lldm_aoi/Basin_LldM.shp',
    'workspace_dir': '../lldm_aoi/output',
    'reference_period_dates': ('1980-01-01', '2010-12-31'),
    'prediction_dates': ('2015-01-01', '2100-12-31'),
    'hindcast': True,
    'gcm_experiment_list': knn.GCM_EXPERIMENT_LIST,
    'gcm_model_list': knn.MODEL_LIST,
    'upper_precip_percentile': 75,
    'lower_precip_threshold': 1,  # millimeter
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
