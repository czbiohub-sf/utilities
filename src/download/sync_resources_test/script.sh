#!/bin/env bash

## VIASH START
par_input='s3://czbiohub-pipelines/resources_test'
par_output='resources_test'
par_dryrun='false'
## VIASH END

extra_params=( )

if [ "$par_dryrun" == "true" ]; then 
  extra_params+=( "--dryrun" )
fi

aws s3 sync "$par_input" "$par_output" --no-sign-request "${extra_params[@]}"
