#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

bin/nextflow \
  run . \
  -main-script workflows/mapping_cellranger/main.nf \
  -resume \
  -with-docker \
  --id tiny_fastq \
  --input resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq \
  --reference resources_test/cellranger_tiny_fastq/cellranger_tiny_ref \
  --publishDir temp