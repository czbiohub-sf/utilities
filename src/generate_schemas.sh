#!/bin/bash

generate_schema \
  --input src/mapping/process_10x/auto.vsh.yaml \
  --output src/mapping/process_10x/nextflow_auto_schema.json
generate_schema \
  --input src/mapping/process_smartseq2/auto.vsh.yaml \
  --output src/mapping/process_smartseq2/nextflow_auto_schema.json