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
         observed_precip_filepath, aoi_netcdf_path, reference_period_dates,
         hindcast, target_filename):
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
        data_series_list = ['python_prediction', 'mswep_precip']
    else:
        data = bootstrapped_data
        data_series_list = ['python_prediction']
        # For forecasts, only the reference period of the observed data
        # is relevant
        mswep_df = mswep_df[reference_period_dates[0]: reference_period_dates[1]]

    data['month'] = data.index.month
    data_long = pandas.melt(
        data[['python_prediction', 'month']], id_vars='month')
    mswep_df['month'] = mswep_df.index.month
    mswep_long = pandas.melt(mswep_df, id_vars='month')
    long_df = pandas.concat([data_long, mswep_long])
    long_df.reset_index(drop=True, inplace=True)

    with PdfPages(target_filename) as pdf_pages:
        fig, axs = plt.subplots(nrows=2, figsize=(12, 6), sharex=hindcast)
        data.plot(
            y=data_series_list,
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
        else:
            axs[1].axis('off')
        pdf_pages.savefig()
        plt.close()

        fig, axs = plt.subplots(ncols=2, figsize=(8, 4))
        seaborn.ecdfplot(
            long_df,
            x='value',
            hue='variable',
            linewidth=0.7,
            ax=axs[0])
        seaborn.move_legend(axs[0], loc='lower center')
        seaborn.despine()
        axs[0].set_title(
            'bootstrapped values', fontsize=10)
        if 'residual' in data.columns:
            seaborn.histplot(
                data,
                x='residual',
                bins=100,
                linewidth=None,
                ax=axs[1])
            seaborn.despine()
            axs[1].set_title('residual')
        else:
            axs[1].axis('off')
        pdf_pages.savefig()
        plt.close()

        nrows = 2
        ncols = 6
        fig, axs = plt.subplots(
            nrows=nrows, ncols=ncols, figsize=(12, 4), sharex=True, sharey=True)
        for r in range(nrows):
            for c in range(ncols):
                month = r*6 + c + 1
                seaborn.ecdfplot(
                    long_df[long_df['month'] == month],
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
            obs_dataset = obs_dataset.sel(
                time=slice(*reference_period_dates))
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
            f'({reference_period_dates[0]} : {reference_period_dates[1]})',
            fontsize=10)
        pdf_pages.savefig()
        plt.close()
