#!/bin/bash

echo Running $meta_functionality_name with standard options
./$meta_functionality_name \
  --id 195891710 --output test_output/

[ ! -f test_output/SampleSheet.csv ] && echo "Output 'SampleSheet.csv' could not be found" && exit 1
[ ! -f test_output/RunParameters.xml ] && echo "Output 'RunParameters.xml' could not be found" && exit 1

if ! grep -q 'RNAEnrichment_RVOPv2_iSeq_2' test_output/SampleSheet.csv; then
  echo Could not find content
  exit 1
fi

echo Running $meta_functionality_name only for csv files
./$meta_functionality_name \
  --id 195891710 --output test_output2/ --extension csv

[ ! -f test_output2/SampleSheet.csv ] && echo "Output 'SampleSheet.csv' could not be found" && exit 1
[ -f test_output2/RunParameters.xml ] && echo "Output 'RunParameters.xml' should not be found" && exit 1

echo Test succeeded!