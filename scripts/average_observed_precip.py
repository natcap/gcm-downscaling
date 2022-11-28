#!/usr/bin/env python
# coding: utf-8

""" This script pre-processes daily observed precip values at a
set of sites in an AOI. It averages the daily values across all
sites. It exports the result in CSV and NetCDF formats.
"""

import os

import pandas
import xarray


data_store_path = 'H://Shared drives/GCM_Climate_Tool/required_files'
observed_precip_path = os.path.join(
    data_store_path, 'OBSERVATIONS/LLdM_AOI2/series_pr_diario.csv')

target_csv_path = os.path.join(
        data_store_path,
        'OBSERVATIONS/LLdM_AOI2/series_pr_diario_regional_average.csv')
target_netcdf_path = os.path.join(
    data_store_path,
    'OBSERVATIONS/LLdM_AOI2/series_pr_diario_regional_average.nc')
var_name = 'regional_pr'

# Make a CSV with the same formatting as original, but with a column
# for the average precip across all sites:
data = pandas.read_csv(observed_precip_path)
sites = [col for col in data.columns if col.startswith('pr_')]
data[var_name] = data[sites].mean(axis=1)
data[['Year', 'Month', 'Day', var_name]].to_csv(target_csv_path)

# Re-format as a netcdf matching CMIP conventions:
data = pandas.read_csv(
        target_csv_path,
        usecols=['Year', 'Month', 'Day', var_name],
        parse_dates={'time': ['Year', 'Month', 'Day']})
data.set_index('time', inplace=True)
data.rename({var_name: 'pr'}, axis=1, inplace=True)
ds = xarray.Dataset.from_dataframe(data)
ds.to_netcdf(target_netcdf_path)
