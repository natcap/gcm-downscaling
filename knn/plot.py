import logging
import numpy
import pandas
import seaborn
import xarray

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

LOGGER = logging.getLogger(__name__)


def plot(dates_filepath, precip_filepath, observed_mean_precip_filepath,
         observed_precip_filepath, aoi_netcdf_path, reference_start_date,
         reference_end_date, target_filename):
    LOGGER.info(f'creating report for {precip_filepath}')
    data = pandas.read_csv(dates_filepath, parse_dates={'date': [0]})
    data.set_index('date', inplace=True)
    data.drop(columns=[
        'historic_date',
        'wet_state',
        'next_wet_state',
        'next_historic_date'], inplace=True)
    data.rename({'historic_precip': 'python_prediction'}, axis=1, inplace=True)

    with PdfPages(target_filename) as pdf_pages:
        plt.figure(figsize=(12, 3))
        seaborn.lineplot(
            data=data, x='date', y='python_prediction',
            linewidth=0.5
        )
        plt.title('bootstrapped values from observed regional average', fontsize=10)
        pdf_pages.savefig()
        plt.close()

        fig, axs = plt.subplots(ncols=2, figsize=(12, 6), sharey=True)
        seaborn.histplot(data, x='python_prediction', bins=100, ax=axs[0])
        axs[0].set_title('bootstrapped values from observed regional average', fontsize=10)
        with xarray.open_dataset(observed_mean_precip_filepath) as obs_mean_dataset:
            obs_mean_dataset = obs_mean_dataset.sortby('time')
            obs_mean_dataset = obs_mean_dataset.sel(time=slice(reference_start_date, reference_end_date))
            obs_mean_dataset.precipitation.plot.hist(bins=100, ax=axs[1])
        axs[1].set_title(
            f'observed regional average values ({reference_start_date} : {reference_end_date})')
        pdf_pages.savefig()
        plt.close()

        levels = numpy.arange(0, 20, 1)
        fig, axs = plt.subplots(ncols=2, figsize=(12, 5), sharey=True)
        with xarray.open_dataset(precip_filepath) as dataset:
            min_date = str(dataset.time.min().values)[:10]
            max_date = str(dataset.time.max().values)[:10]
            dataset.precipitation.mean(['time']).plot(ax=axs[0], levels=levels)
        axs[0].set_title(
            f'downscaled average daily precip ({min_date} : {max_date})', fontsize=10)
        with xarray.open_dataset(observed_precip_filepath) as obs_dataset:
            obs_dataset = obs_dataset.sortby('time')
            obs_dataset = obs_dataset.sel(time=slice(reference_start_date, reference_end_date))
            with xarray.open_dataset(aoi_netcdf_path) as aoi_dataset:
                obs_dataset['aoi'] = aoi_dataset.aoi
            obs_dataset = obs_dataset.where(obs_dataset.aoi == 1, drop=True)
            obs_dataset.precipitation.mean(['time']).plot(ax=axs[1], levels=levels)
        axs[1].set_title(
            f'observed average daily precip ({reference_start_date} : {reference_end_date})', fontsize=10)
        pdf_pages.savefig()
        plt.close()
