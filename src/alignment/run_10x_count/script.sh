#!/bin/bash

# process parameters

SAMPLE_ID=`basename $par_input_path`
FASTQ_PATH=
TRANSCRIPTOME_DIR
SAMPLE_NAME
# TODO: define other variables

extra_params=( )

if [ -z "$par_cores" ]; then 
  extra_params+=( "--localcores" "$par_cores" )
fi
if [ -z "$par_memory" ]; then 
  extra_params+=( "--localmem" "$par_memory" )
fi

cellranger \
  count \
  --nosecondary \
  --disable-ui \
  --expect-cells $par_cell_count \
  --id $SAMPLE_ID \
  --fastqs $FASTQ_PATH \
  --transcriptome $TRANSCRIPTIME_DIR \
  --sample $SAMPLE_NAME \
  "${extra_params[@]}"

# output checks