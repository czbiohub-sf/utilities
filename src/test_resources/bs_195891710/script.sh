#!/bin/bash

DIR=resources_test/bs_195891710
ID=`basename "$OUT"`
S3DIR=`echo "$DIR" | sed 's#test_resources#s3://czbiohub-pipelines#'`

# download bcl files from basespace
target/docker/download/download_basespace/download_basespace \
  --id 195891710 \
  --output "$DIR/bcl_data/"

# convert to fastq
target/docker/demux/bcl2fastq/bcl2fastq \
  --input "$DIR/bcl_data/" \
  --output "$DIR/fastqs/"