#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=pbmc_1k_protein_v3
OUT="resources_test/$ID/"
DIR="$OUT"
S3DIR="s3://czbiohub-pipelines/$DIR"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -O resources_test/reference/refdata-gex-GRCh38-2020-A.tar.gz

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_fastqs.tar -O resources_test/reference/pbmc_1k_protein_v3_fastqs.tar

# todo: untar

if [ ! -f "${OUT}/bam" ]; then
  mkdir -p "$OUT/bam"

  target/docker/alignment/cellranger_count/cellranger_count \
    --input "${OUT}/fastqs" \
    --reference "resources_test/reference/refdata-gex-GRCh38-2020-A" \
    --output "${OUT}/bam"
fi

H5AD="${OUT}/filtered_feature_bc_matrix.h5ad"
if [ ! -f "$H5AD" ]; then
  target/docker/convert/convert_10x_h5_to_h5ad/convert_10x_h5_to_h5ad \
    --input "${OUT}/bam/filtered_feature_bc_matrix.h5" \
    --output "$H5AD"
fi

H5MU="${OUT}/filtered_feature_bc_matrix.h5mu"
if [ ! -f "$H5MU" ]; then
  target/docker/convert/convert_10x_h5_to_h5mu/convert_10x_h5_to_h5mu \
    --input "${OUT}/bam/filtered_feature_bc_matrix.h5" \
    --output "$H5MU"
fi