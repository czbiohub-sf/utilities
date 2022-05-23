#!/bin/bash

DIR=resources_test/bs_195891710
ID=`basename "$OUT"`

# download bcl files from basespace
target/docker/download/download_basespace/download_basespace \
  --id 195891710 \
  --output "$DIR/bcl/"

# convert to fastq
target/docker/demux/bcl2fastq/bcl2fastq \
  --input "$DIR/bcl/" \
  --output "$DIR/fastqs/"
