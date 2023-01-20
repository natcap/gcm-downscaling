import argparse
import os
import re
import shutil

from dask.distributed import Client, progress
import rechunker
import xarray
from zarr.convenience import consolidate_metadata


mswep_store = '/oak/stanford/groups/gdaily/mswep2/annual'
zarr_store = os.path.join(mswep_store, '..', 'zarr')
netcdf_files = [
    os.path.join(mswep_store, f)
    for f in os.listdir(mswep_store)
    if re.search(r'[0-9]{4}.nc', f)]

if not os.path.exists(zarr_store):
    os.mkdir(zarr_store)


def preprocessor(ds):
    return ds.sortby('time')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--n_workers', type=int, default=-1,
        help='number of workers for Taskgraph.')
    parser.add_argument(
        '--max_mem', type=int, default=4,
        help='max memory in GB for a dask worker')
    args = parser.parse_args()

    client = Client(n_workers=args.n_workers,
                    processes=True,
                    threads_per_worker=1,
                    memory_limit=str(args.max_mem) + 'GB')

    with xarray.open_mfdataset(
            netcdf_files, parallel=True, combine='nested', concat_dim='time',
            data_vars='minimal', coords='minimal', compat='override',
            preprocess=preprocessor, autoclose=True,
            chunks={'lon': 360, 'lat': 360}) as dataset:
        # Cannot chunk the concat dim on opening
        dataset = dataset.chunk({'time': dataset.time.size})

        target_store = os.path.join(zarr_store, 'mswep_1980_2020.zarr')
        temp_store = os.path.join(zarr_store, 'temp', 'chunking.zarr')
        target_chunks = {
            'precipitation': {
                'time': len(dataset.time),
                'lon': 90,
                'lat': 90
            },
            'time': None,  # don't rechunk these
            'lon': None,
            'lat': None,
        }

        if os.path.exists(target_store):
            shutil.rmtree(target_store)
        if os.path.exists(temp_store):
            shutil.rmtree(temp_store)
        array_plan = rechunker.rechunk(
            dataset,
            target_chunks,
            str(args.max_mem / 2) + 'GB',  # observed dask actually using 2x max_mem
            target_store,
            temp_store=temp_store)

        print(array_plan)
        future = array_plan.execute()
        progress(future)
        shutil.rmtree(temp_store)
        consolidate_metadata(target_store)


if __name__ == '__main__':
    main()
