#!/bin/bash

set -e

# settings
ID=pbmc_1k_v3
DIR=resources_test/$ID

# get fastqs
mkdir -p "$DIR/fastqs"
wget "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar" \
  -O "$DIR/temp_pbmc_1k_v3_fastqs.tar"
tar -xf "$DIR/temp_pbmc_1k_v3_fastqs.tar" -C "$DIR/fastqs"

