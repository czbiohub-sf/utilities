#!/bin/bash

## VIASH START
meta_resources_name=foo
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_resources_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# if par_input not is a folder, untar first
if [ ! -d "$par_input" ]; then  
  echo "Assuming input is a tar.gz, untarring"
  input_dir="$tmpdir/bcl"
  mkdir -p "$input_dir"
  tar -xzf "$par_input" -C "$input_dir" --strip-components=1
else
  input_dir="$par_input"
fi

# create output dir, if need be
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

# add additional params
extra_params=( )

if [ ! -z "$par_cores" ]; then 
  extra_params+=( "--localcores" "$par_cores" )
fi
if [ ! -z "$par_memory" ]; then 
  extra_params+=( "--localmem" "$par_memory" )
fi


echo "Running cellranger demux"
work_dir="$tmpdir/out"
mkdir -p "$work_dir"

id=myoutput

cellranger mkfastq \
  --id "$id" \
  --sample-sheet "$par_sample_sheet" \
  --run "$par_input" \
  --output-dir "$work_dir" \
  "${extra_params[@]}" \
  --disable-ui

if [ -d "$work_dir/$id/outs/" ]; then
  mv "$work_dir/$id/outs/fastq_path" "$par_output/"
fi
