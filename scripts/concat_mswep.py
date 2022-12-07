from collections import defaultdict
import os
import re

import taskgraph
import xarray

data_store_path = 'H://Shared drives/GCM_Climate_Tool/required_files'
mswep_store = os.path.join(data_store_path, 'OBSERVATIONS', 'Global MSWEP2', 'Daily')


def concat_netcdfs(year, target_path):
    print(year)
    file_list = [os.path.join(mswep_store, f) for f in filemap[year]]
    with xarray.open_mfdataset(
            file_list, parallel=True, combine='nested', concat_dim='time',
            data_vars='minimal', coords='minimal', compat='override') as dataset:
        dataset.to_netcdf(target_path)


workspace_dir = 'C:/Users/dmf/projects/gcm-project/mswep_annual'
taskgraph_working_dir = os.path.join(workspace_dir, '..', '.taskgraph')
graph = taskgraph.TaskGraph(taskgraph_working_dir, -1)
if not os.path.exists(workspace_dir):
    os.mkdir(workspace_dir)

filemap = defaultdict(list)
for file in os.listdir(mswep_store):
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
            'target_path': target_path
        },
        target_path_list=[target_path],
        dependent_task_list=[]
    )
