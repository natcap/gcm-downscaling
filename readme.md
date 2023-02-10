## Setup

1. Request access to climate data
    * Send Dave a google account email address you wish to use for Google Cloud.  
    (gmail or @stanford.edu addresses should work)
    * _(Users will be granted Role `Storage Object Viewer` for bucket `natcap-climate-data`)_

2. Authenticate with Google Cloud
    * install `gcloud` if needed (https://cloud.google.com/sdk/docs/install)
    * `gcloud config set project natcap-server`
    * `gcloud auth application-default login`

## Usage option 1: with a conda environment

1. Setup the python environment
`git clone https://github.com/natcap/gcm-downscaling.git`
`cd gcm-downscaling`
`conda env create -p ./env --file requirements.yml`  
`conda activate ./env`
`pip install .`

2. Make a copy of `example_run.py` and modify the `args` dictionary

3. `python copy_of_example_run.py`

## Usage option 2: with a docker container

1. Setup a directory containing your AOI vector and a copy of `example_run.py`
2. modify the `args` dictionary in your copy of `example_run.py`. 
    * set `args['workspace_dir']` to be a relative path within this directory
3. 
`docker run --rm -ti 
-v %CD%:/workspace 
-v %appdata%/gcloud:/home/mambauser/.config/gcloud 
-w /workspace 
-e GOOGLE_CLOUD_PROJECT='natcap-servers'
ghcr.io/natcap/gcm-downscaling:latest python copy_of_example_run.py`

Note: You may need to give Docker more RAM and CPUs than it is allowed to use by default.  
Adjust in `Docker Desktop > Settings > Resources`. 6GB of RAM should do it.

### The `args` dictionary
**'aoi_path**' (str): a path to a GDAL polygon vector. Coordinates
    represented by longitude, latitude decimal degrees (WGS84).
**'workspace_dir'** (str): a path to the directory where this program
    writes output and other temporary files.
**'reference_period_dates'** (sequence): ('YYYY-MM-DD', 'YYYY-MM-DD')
    first and last day in the reference period, which is used to
    calculate climate "normals".
**'prediction_dates'** (sequence, optional): ('YYYY-MM-DD', 'YYYY-MM-DD')
    first and last day in the simulation period.
    Required if `hindcast=False`.
**'lower_precip_threshold'** (float): the lower boundary of the
    middle bin of precipitation states. Units: mm
**'upper_precip_percentile'** (float): a percentile (from 0:100) with
    which to extract the absolute precipitation value that will be the
    upper boundary (inclusive) of the middle bin of precipitation states.
**'hindcast'** (bool): If True, observed data (MSWEP) is substituted
    for GCM data and the prediction period is set to match the date
    range of the observed dataset (``knn.MSWEP_DATE_RANGE``).
**'gcm_model_list'** (sequence, optional): a sequence of strings
    representing CMIP6 model codes. Each model will be used to generate
    a single downscaled product for each experiment in `gcm_experiment_list`.
    Available models are stored in ``knn.GCM_MODEL_LIST``.
    Required if `hindcast=False`.
**'gcm_experiment_list'** (sequence, optional): a sequence of strings
    representing CMIP6 SSP experiments. Available experiments are
    stored in ``GCM_EXPERIMENT_LIST``. If a CMIP model does not include
    a given experiment, that experiment will be skipped for that model.
    Required if `hindcast=False`.
**'n_workers'** (int, optional): The number of worker processes to
    use. If omitted, computation will take place in the current process.
    If a positive number, tasks can be parallelized across this many
    processes, which can be useful if `gcm_model_list` or
    `gcm_experiement_list` contain multiple items.

# About Global Climate Data
## Data Availablity

This workflow derives downscaled climate data from,
* CMIP6 General Circulation Models (ee `knn.MODEL_LIST` for list of available models)
* MSWEP historical precipitation data.

## Data Storage
Analysis-ready data are stored in `zarr` format in a private google cloud bucket
(`natcap-climate-data`) in the `NatCap Servers` cloud project.  

Raw netCDF data are stored on Stanford's Oak Storage Service at
`/oak/stanford/groups/gdaily/`.  
See `scripts/preprocessing/` for workflows to create `zarr` from `netCDF`.

## Adding new data sources
### global data:
New CMIP6 models or other new global data source can be made available
by following the examples in `scripts/preprocessing/` to create `zarr` stores
and move them to the `natcap-climate-data` bucket. 

### local data:
To use this workflow with local observational data instead of MSWEP data,
the data should be formatted as a netCDF (or other `xarray.open_dataset` readable format).

The netCDF should contain coordinates and variables named & defined as,

`Coordinates:`
* `lat`  - decimal degrees (-90 : 90)
* `lon`  - decimal degrees (-180 : 180) or (0 : 360)
* `time` - daily timesteps in units that can be parsed to `numpy.datetime64`

`Variables:`
* `precipitation` - dimensions: `(time, lat, lon)`; units: millimeter

The downscaled product will have the same spatial resolution as the observation data.
