#!/bin/bash
#
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=7000M
#SBATCH --mail-type=ALL
#SBATCH --mail-user=davefisher@stanford.edu
#SBATCH --partition=hns,normal
#SBATCH --output=/scratch/users/dfisher5/slurm-logfiles/slurm-%j.%x.out
#
# --partition=hns,normal means that this will be submitted to both queues, whichever gets to it first will be used.

CONTAINER=ghcr.io/natcap/gcm-downscaling:latest

WORKSPACE_DIR="$L_SCRATCH/$WORKSPACE_NAME"

set -x  # Be eXplicit about what's happening.
FAILED=0
singularity run \
    docker://$CONTAINER python scripts/rechunk_to_zarr_cmip6.py \
    --n_workers=32 \
    --max_mem=7
