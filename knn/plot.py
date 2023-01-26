import logging
import math

import numpy
import pandas
import seaborn
import xarray

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

LOGGER = logging.getLogger(__name__)


def plot(dates_filepath, precip_filepath, observed_mean_precip_filepath,
         observed_precip_filepath, aoi_netcdf_path, reference_start_date,
         reference_end_date, hindcast, target_filename):
    LOGGER.info(f'creating report for {precip_filepath}')

    bootstrapped_data = pandas.read_csv(
        dates_filepath, parse_dates={'date': [0]})
    bootstrapped_data.set_index('date', inplace=True)
    bootstrapped_data.drop(columns=[
        'historic_date',
        'wet_state',
        'next_wet_state',
        'next_historic_date'], inplace=True)
    bootstrapped_data.rename(
        {'historic_precip': 'python_prediction'}, axis=1, inplace=True)

    mswep_mean = xarray.open_dataset(observed_mean_precip_filepath)
    mswep_series = mswep_mean.precipitation.to_pandas()
    mswep_df = pandas.DataFrame({'mswep_precip': mswep_series})

    if hindcast:
        # Forecast and Observed timeseries overlap, can be merged        
        data = pandas.merge(
            bootstrapped_data, mswep_df,
            left_index=True, right_index=True, how='left')
        data['residual'] = data['python_prediction'] - data['mswep_precip']
    else:
        data = bootstrapped_data
        # For forecasts, only the reference period of the observed data
        # is relevant
        mswep_df = mswep_df[reference_start_date: reference_end_date]

    with PdfPages(target_filename) as pdf_pages:
        fig, axs = plt.subplots(nrows=2, figsize=(12, 6), sharex=True)
        fig.frameon = False
        data.plot(
            y=['python_prediction', 'mswep_precip'],
            linewidth=0.5,
            alpha=0.7,
            ax=axs[0])
        axs[0].set_title(
            'bootstrapped values from observed regional average', fontsize=10)
        axs[0].grid(True, axis='y', linestyle='dotted', linewidth=0.1)
        if 'residual' in data.columns:
            data.plot(
                y='residual',
                linewidth=0.2,
                alpha=0.7,
                ax=axs[1])
            axs[1].grid(True, axis='y', linestyle='dotted', linewidth=0.1)
            pdf_pages.savefig()
            plt.close()

        fig, axs = plt.subplots(ncols=3, figsize=(12, 4), sharey=True)
        fig.frameon = False
        seaborn.histplot(
            data,
            x='python_prediction',
            bins=100,
            linewidth=None,
            ax=axs[0])
        seaborn.despine()
        axs[0].set_title(
            'bootstrapped values', fontsize=10)
        seaborn.histplot(
            mswep_df,
            x='mswep_precip',
            bins=100,
            linewidth=None,
            ax=axs[1])
        seaborn.despine()
        axs[1].set_title(
            f'observed regional avg \n'
            f'({str(mswep_df.index.min())[:10]} : '
            f'{str(mswep_df.index.max())[:10]})')
        if 'residual' in data.columns:
            seaborn.histplot(
                data,
                x='residual',
                bins=100,
                linewidth=None,
                ax=axs[2])
            seaborn.despine()
            axs[2].set_title('residual')
        pdf_pages.savefig()
        plt.close()

        data['month'] = data.index.month
        data_long = pandas.melt(
            data.drop(columns=['residual']), id_vars="month")
        nrows = 2
        ncols = 6
        fig, axs = plt.subplots(
            nrows=nrows, ncols=ncols, figsize=(12, 4), sharex=True, sharey=True)
        fig.frameon = False
        for r in range(nrows):
            for c in range(ncols):
                month = r*6 + c + 1
                seaborn.ecdfplot(
                    data_long[data_long['month'] == month],
                    x="value",
                    hue="variable",
                    linewidth=0.5,
                    ax=axs[r, c],
                    legend=False)
                seaborn.despine()
                axs[r, c].set_title(f'month {month}')
                axs[r, c].xaxis.label.set_visible(False)
        pdf_pages.savefig()
        plt.close()

        fig, axs = plt.subplots(ncols=2, figsize=(12, 5), sharey=True)
        fig.frameon = False
        with xarray.open_dataset(precip_filepath) as dataset:
            min_date = str(dataset.time.min().values)[:10]
            max_date = str(dataset.time.max().values)[:10]
            downscaled_mean = dataset.precipitation.mean(['time'])
            max_precip = downscaled_mean.max()
        with xarray.open_dataset(observed_precip_filepath) as obs_dataset:
            obs_dataset = obs_dataset.sortby('time')
            obs_dataset = obs_dataset.sel(time=slice(reference_start_date, reference_end_date))
            with xarray.open_dataset(aoi_netcdf_path) as aoi_dataset:
                obs_dataset['aoi'] = aoi_dataset.aoi
            obs_dataset = obs_dataset.where(obs_dataset.aoi == 1, drop=True)
            obs_mean = obs_dataset.precipitation.mean(['time'])
            obs_max = obs_mean.max()
            if obs_max > max_precip:
                max_precip = obs_max

        levels = numpy.arange(0, math.ceil(max_precip) + 1, 1)
        downscaled_mean.plot(ax=axs[0], levels=levels)
        axs[0].set_title(
            f'downscaled average daily precip\n'
            f'({min_date} : {max_date})',
            fontsize=10)
        obs_mean.plot(ax=axs[1], levels=levels)
        axs[1].set_title(
            f'observed average daily precip\n'
            f'({reference_start_date} : {reference_end_date})',
            fontsize=10)
        pdf_pages.savefig()
        plt.close()
