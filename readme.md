## Setup

1. Request access to climate data
* Send Dave a google account email address you wish to use for Google Cloud.  
(gmail or @stanford.edu addresses should work)
2. Authenticate with Google Cloud
* install `gcloud` if needed (https://cloud.google.com/sdk/docs/install)
* `gcloud config set project natcap-server`
* `gcloud auth application-default login`

## Usage, Option 1: with a conda environment

1. Create and activate a python environment
`conda env create -p ./env --file requirements.yml`  
`conda activate ./env`

2. Make a copy of `example_run.py` and modify the `args` dictionary

3. `python example_run.py`

## Usage, Option 2: with a docker container

1. Setup a directory containing your AOI vector and a copy of `example_run.py`
2. modify the `args` dictionary in `example_run.py`. 
* set `workspace_dir` to be a relative path within this directory
3. `docker run --rm -ti -v %CD%:/workspace -w /workspace ghcr.io/natcap/gcm-downscaling:latest python copy_of_example_run.py`

# About Global Climate Data
## Data Availablity

This workflow derives downscaled climate data from,
* CMIP6 General Circulation Models
-- See `knn.MODEL_LIST` for list of available models
* MSWEP historical precipitation data.

## Data Storage
Analysis-ready data are stored in `zarr` format in a private google cloud bucket
(`natcap-climate-data`) in the `NatCap Servers` cloud project.  

Raw netCDF data are stored on Stanford's Oak Storage Service at
`/oak/stanford/groups/gdaily/`.  
See `scripts/preprocessing/` for workflows to create `zarr` from `netCDF`.

## Adding new data sources
### global data:
New CMIP6 models, or any other new global data source can be made available
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

