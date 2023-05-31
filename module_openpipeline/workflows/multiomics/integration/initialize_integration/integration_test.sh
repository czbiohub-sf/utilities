#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow run . \
  -main-script workflows/multiomics/integration/initialize_integration/main.nf \
  -profile docker,no_publish \
  -entry test_wf \
  -with-trace work/trace.txt \
  -with-dag workflows/multiomics/integration/initialize_integration/graph.dot