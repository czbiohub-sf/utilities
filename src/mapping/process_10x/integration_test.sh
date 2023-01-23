#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.04.5

nextflow \
  run . \
  -main-script workflows/process/main.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/wf_utils/labels_ci.config