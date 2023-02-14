#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.10.4

# nextflow \
#   run . \
#   -main-script src/mapping/process_10x/main.nf \
#   -entry test_wf \
#   -resume \
#   -profile docker,no_publish \
#   -c src/wf_utils/labels_ci.config


# nextflow \
#   run . \
#   -main-script src/mapping/process_10x/main.nf \
#   -resume \
#   -profile docker \
#   --id pbmc_1k_v3 \
#   --input $RES_DIR/ \
#   --reference resources_test/reference/gencode_v41_ercc_cellranger

nextflow run . \
  -main-script src/mapping/process_10x/main.nf \
  -entry auto \
  -resume \
  -profile docker \
  --input_dir resources_test/pbmc_1k_v3/fastqs/pbmc_1k_v3_fastqs \
  --reference resources_test/reference/gencode_v41_ercc_cellranger \
  --publish_dir output \
  -stub



nextflow run czbiohub/utilities \
  -r implement_auto_workflows \
  -main-script src/mapping/process_10x/main.nf \
  -entry auto \
  -resume \
  -profile docker \
  --input_dir resources_test/pbmc_1k_v3/fastqs/pbmc_1k_v3_fastqs \
  --reference resources_test/reference/gencode_v41_ercc_cellranger \
  --publish_dir output \
  --dry_run