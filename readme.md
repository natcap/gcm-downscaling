## Usage with a conda environment

1. Create and activate a python environment
`conda env create -p ./env --file requirements.yml`  
`conda activate ./env`

2. Make a copy of `example_run.py` and modify the `args` dictionary

3. `python example_run.py`

## Usage with a docker container

1. Setup a directory containing your AOI vector and a copy of `example_run.py`
2. modify the `args` dictionary in `example_run.py`. 
* set `workspace_dir` to be a relative path within this directory
3. docker run --rm -ti -v %CD%:/workspace -w /workspace ghcr.io/natcap/gcm-downscaling:latest python copy_of_example_run.py

## Run Tests
`pytest src/tests`
