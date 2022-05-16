#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR="$HOME/.cache/singularity"

bin/nextflow \
  run . \
  -main-script workflows/demux_cellranger/main.nf \
  -resume \
  -with-singularity \
  --id tiny_bcl \
  --input resources_test/cellranger_tiny_bcl/bcl \
  --sample_sheet resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv \
  --publishDir temp