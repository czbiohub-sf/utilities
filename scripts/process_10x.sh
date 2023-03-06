#!/bin/bash

set -e

group=group.data.science

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# check whether access token is available
if [ -z $TOWER_ACCESS_TOKEN ]; then
  echo Could not find a tower access token!
  exit 1
fi

export TOWER_WORKSPACE_ID=124983787544228
export NXF_VER=22.04.5
export NXF_SINGULARITY_CACHEDIR="/hpc/scratch/$group/singularity_images"
export SCRATCH_DIR="/hpc/scratch/$group/nextflow_$USER"

# temp settings
export APPTAINER_CACHEDIR="$SCRATCH_DIR/temp"
export NXF_TEMP="$SCRATCH_DIR/temp"
export APPTAINERENV_TMPDIR="$SCRATCH_DIR/temp"

# dry run:
(cd "$SCRATCH_DIR" && nextflow run czbiohub/utilities \
  -r main_build \
  -main-script src/mapping/process_smartseq2/main.nf \
  -entry auto \
  -resume \
  -latest \
  -profile singularity \
  -c $SCRIPT_DIR/czbhpc.config \
  $@)