#!/bin/bash


# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_resources_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT


# add additional params
extra_params=( )

if [ ! -z "$par_cores" ]; then 
  extra_params+=( "--localcores" "$par_cores" )
fi
if [ ! -z "$par_memory" ]; then 
  extra_params+=( "--localmem" "$par_memory" )
fi
if [ ! -z "$par_cell_count" ]; then 
  extra_params+=( "--expect-cells" "$par_cell_count" )
fi
if [ ! -z "$par_libraries" ]; then 
  extra_params+=( "--libraries" "$par_libraries" )
fi
if [ ! -z "$par_feature_ref" ]; then 
  extra_params+=( "--feature-ref" "$par_feature_ref" )
fi

echo Running cellranger count
id=myoutput

cellranger count \
  --id "$id" \
  --fastqs "$par_input" \
  --transcriptome "$par_transcriptome" \
  "${extra_params[@]}" \
  --disable-ui



# create output dir
if [ ! -d "$par_output" ]; then
  mkdir -p "$par_output"
fi

# copy output
if [ -d "$id/outs/" ]; then
  mv "$id"/outs/* "$par_output/"
fi