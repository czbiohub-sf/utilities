#!/bin/bash

set -e

# pipeline settings
WF_REPOSITORY="czbiohub/utilities"
WF_MAIN_SCRIPT="src/mapping/process_10x/main.nf"
WF_REVISION="main_build"
WF_ENTRY="auto"

# data science settings
group=group.data.science
export TOWER_WORKSPACE_ID=124983787544228

# used for finding the czbhpc.config
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# check whether access token is available
if [ -z $TOWER_ACCESS_TOKEN ]; then
  echo Could not find a tower access token!
  exit 1
fi

# workdir settings
export NXF_VER=22.04.5
export NXF_SINGULARITY_CACHEDIR="/hpc/scratch/$group/singularity_images"
export SCRATCH_DIR="/hpc/scratch/$group/nextflow_$USER"

# temp settings
export APPTAINER_CACHEDIR="$SCRATCH_DIR/temp"
export APPTAINERENV_TMPDIR="$SCRATCH_DIR/temp"
export NXF_TEMP="$SCRATCH_DIR/temp"
export NXF_WORK="$SCRATCH_DIR/work"

# create scratch ddir if it does not exist
[ -d "$SCRATCH_DIR" ] || mkdir -p "$SCRATCH_DIR"

(cd "$SCRATCH_DIR" && nextflow run  \
  "$WF_REPOSITORY" \
  -r "$WF_REVISION" \
  -main-script "$WF_MAIN_SCRIPT" \
  -entry "$WF_ENTRY" \
  -resume \
  -latest \
  -profile singularity \
  -c $SCRIPT_DIR/czbhpc.config \
  $@)