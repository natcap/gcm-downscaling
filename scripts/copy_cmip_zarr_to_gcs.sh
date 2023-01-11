#!/bin/bash
#
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4000M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davefisher@stanford.edu
#SBATCH --partition=hns,normal
#SBATCH --output=/scratch/users/dfisher5/slurm-logfiles/slurm-%j.%x.out
#
# --partition=hns,normal means that this will be submitted to both queues, whichever gets to it first will be used.

ml load system rclone
rclone copy --progress \
    /oak/stanford/groups/gdaily/cmip6/zarr/ \
    gcs-remote:natcap-climate-data/cmip6/
