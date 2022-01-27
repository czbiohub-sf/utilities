#!/bin/bash

# process parameters

# TODO: define other variables

cellranger \
  count \
  --localmem=240 \
  --nosecondary \
  --disable-ui \
  --expect-cells $par_cell_count \
  --id $SAMPLE_ID \
  --fastqs $FASTQ_PATH \
  --transcriptome $TRANSCRIPTIME_DIR \
  --sample $SAMPLE_NAME

# output checks