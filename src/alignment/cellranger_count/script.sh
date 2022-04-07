#!/bin/bash

## VIASH START
par_input='resources_test/pbmc_1k_protein_v3/fastqs/pbmc_1k_protein_v3_gex_fastqs/'
par_transcriptome='resources_test/reference/refdata-gex-GRCh38-2020-A'
par_output='resources_test/pbmc_1k_protein_v3/bam'
par_input=`realpath $par_input`
par_transcriptome=`realpath $par_transcriptome`
par_output=`realpath $par_output`
## VIASH END


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
if [ ! -z "$par_expect_cells" ]; then 
  extra_params+=( "--expect-cells" "$par_expect_cells" )
fi
# if [ ! -z "$par_libraries" ]; then 
#   extra_params+=( "--libraries" "$par_libraries" )
# fi
# if [ ! -z "$par_feature_ref" ]; then 
#   extra_params+=( "--feature-ref" "$par_feature_ref" )
# fi

echo "Running cellranger count"
orig_pwd=`pwd`
work_dir="$tmpdir/out"
mkdir -p "$work_dir"

id=myoutput

cd "$work_dir"
cellranger count \
  --id "$id" \
  --fastqs "$par_input" \
  --transcriptome "$par_transcriptome" \
  "${extra_params[@]}" \
  --disable-ui \
  --nosecondary \
  --output-dir "$par_output"
cd "$orig_pwd"

# echo "Copying output"
# if [ -d "$work_dir/$id/outs/" ]; then
#   if [ ! -d "$par_output" ]; then
#     mkdir -p "$par_output"
#   fi
#   mv "$work_dir/$id/outs/"* "$par_output"
# fi
