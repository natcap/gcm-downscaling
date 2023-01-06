import argparse
import glob
import multiprocessing
import os
import re
import shutil

from dask.distributed import Client, progress
import rechunker
import xarray

import knn


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--n_workers', type=int, default=multiprocessing.cpu_count(),
        help='number of workers for Taskgraph.')
    parser.add_argument(
        '--max_mem', type=int, default=4,
        help='max memory in GB for a dask worker')
    args = parser.parse_args()

    client = Client(n_workers=args.n_workers,
                    processes=True,
                    threads_per_worker=1,
                    memory_limit=str(args.max_mem) + 'GB')

    # cmip_store = '/oak/stanford/groups/gdaily/cmip6'
    cmip_store = 'H://Shared drives/GCM_Climate_Tool/required_files/GCMs/CMIP6'
    zarr_store = os.path.join(cmip_store, 'zarr')

    filemap = {}
    for model in knn.MODEL_LIST:
        for experiment in knn.GCM_EXPERIMENT_LIST:
            for var in knn.GCM_VAR_LIST:
                nc_files = glob.glob(
                    os.path.join(cmip_store),
                    model, f'{var}_day_{model}_{experiment}_*.nc')
                begin_dates = []
                end_dates = []
                for ncf in nc_files:
                    variant, grid, date_range = os.path.splitext(
                        os.path.basename(ncf))[0].split('_')[-3:]
                    begin, end = [int(x) for x in date_range.split('-')]
                    begin_dates.append(begin)
                    end_dates.append(end)
                zarr_filename = \
                    f'{var}_day_{model}_{experiment}_{variant}_{grid}_{min(begin_dates)}-{max(end_dates)}.zarr'
                filemap[zarr_filename] = nc_files

    # if not os.path.exists(zarr_store):
    #     os.mkdir(zarr_store)
    print(filemap)

    def make_zarr(nc_file_list, target_path):
        with xarray.open_mfdataset(
                nc_file_list, parallel=True, combine='nested', concat_dim='time',
                data_vars='minimal', coords='minimal', compat='override',
                autoclose=True, chunks=-1) as dataset:
            # Cannot chunk the concat dim on opening
            dataset = dataset.chunk({'time': dataset.time.size})

        temp_store = os.path.join(
            zarr_store, 'temp', os.path.basename(target_path))
        # I think some GCM grids are 1deg resolution and some 2deg.
        # size 10 chunks are then 10deg or 20deg chunks, which can
        # typically contain an entire downscaling AOI. So a likely
        # worst case is needing to open 4 chunks to extract for an AOI.
        # And chunks are still only ~50MB for a 150yr timeseries
        target_chunks = {
            'precipitation': {
                'time': len(dataset.time),
                'lon': 10,
                'lat': 10
            },
            'time': None,  # don't rechunk these
            'lon': None,
            'lat': None,
        }

        if os.path.exists(target_path):
            shutil.rmtree(target_path)
        if os.path.exists(temp_store):
            shutil.rmtree(temp_store)
        array_plan = rechunker.rechunk(
            dataset,
            target_chunks,
            str(args.max_mem / 2) + 'GB',  # observed dask actually using 2x max_mem
            target_path,
            temp_store=temp_store)

        print(array_plan)
        future = array_plan.execute()
        progress(future)
        shutil.rmtree(temp_store)

    # taskgraph_working_dir = os.path.join(zarr_store, '.taskgraph')
    # graph = taskgraph.TaskGraph(taskgraph_working_dir, args.n_workers)
    # for zarr in filemap:
    #     target_path = os.path.join(zarr_store, zarr)
    #     graph.add_task(
    #         func=make_zarr,
    #         kwargs={
    #             'nc_file_list': filemap[zarr],
    #             'target_path': target_path
    #         },
    #         target_path_list=[target_path],
    #         dependent_task_list=[]
    #     )


if __name__ == '__main__':
    main()
