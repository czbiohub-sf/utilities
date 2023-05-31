#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

# Same as above but with remote yaml file.
nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.yaml

# Same as above but with remote json file.
nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.json

# Same as above but with remote csv file.
nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.csv

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf2 \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf3 \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

