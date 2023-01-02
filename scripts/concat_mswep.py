import argparse
from collections import defaultdict
import multiprocessing
import os
import re

import taskgraph
import xarray

mswep_store_path = '/oak/stanford/groups/gdaily/mswep2'


def concat_netcdfs(year, filemap, target_path):
    print(year)
    file_list = [os.path.join(mswep_store_path, f) for f in filemap[year]]
    with xarray.open_mfdataset(
            file_list, parallel=False, combine='nested', concat_dim='time',
            data_vars='minimal', coords='minimal', compat='override',
            autoclose=True) as dataset:
        dataset.to_netcdf(target_path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--n_workers', type=int, default=multiprocessing.cpu_count(),
        help='number of workers for Taskgraph.')
    args = parser.parse_args()

    workspace_dir = os.path.join(mswep_store_path, 'annual')
    taskgraph_working_dir = os.path.join(mswep_store_path, '.taskgraph')
    graph = taskgraph.TaskGraph(taskgraph_working_dir, args.n_workers)
    if not os.path.exists(workspace_dir):
        os.mkdir(workspace_dir)

    filemap = defaultdict(list)
    for file in os.listdir(mswep_store_path):
        if re.search(r'[0-9]{7}\.nc', file):
            year = file[:4]
            filemap[year].append(file)
    del filemap['1979']  # an incomplete year

    for year in filemap:
        target_path = os.path.join(workspace_dir, f'{year}.nc')
        graph.add_task(
            func=concat_netcdfs,
            kwargs={
                'year': year,
                'filemap': filemap,
                'target_path': target_path
            },
            target_path_list=[target_path],
            dependent_task_list=[]
        )

    graph.join()


if __name__ == '__main__':
    main()
