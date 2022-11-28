#!/usr/bin/env bash

### 1 ###

## VIASH START

input="components/test_resources/demux/input"
output="components/test_resources/demux/output"
samplesheet=""

## VIASH END
if [ -z "$samplesheet" ]
then
      spaceranger mkfastq --run $input --output-dir $output
else
      spaceranger mkfastq --run $input --samplesheet $samplesheet --output-dir $output
fi