#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.10.4

nextflow run . \
  -main-script src/mapping/process_smartseq2/main.nf \
  -entry auto \
  -resume \
  --input_dir resources_test/tsp10_fat_mat_a \
  --reference_index resources_test/reference/gencode_v41_ercc_star \
  --reference_gtf resources_test/reference/gencode_v41_ercc.gtf.gz \
  --publish_dir output2 \
  -stub
