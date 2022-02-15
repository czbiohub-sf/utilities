#!/bin/bash

## VIASH START
par_id='195891710'
par_output='output/'
par_verbose='false'
par_compress='false'
par_overwrite='false'
#par_extension=''
par_api_server=https://api.basespace.illumina.com
par_access_token=610b74406284445fa38691d551c30e4a
## VIASH END

extra_params=( )

if [ "$par_verbose" == "true" ]; then 
  extra_params+=( "--verbose" )
fi
if [ "$par_compress" == "true" ]; then 
  extra_params+=( "--compress" )
fi
if [ "$par_overwrite" == "true" ]; then 
  extra_params+=( "--overwrite" )
fi
if [ ! -z "$par_extension" ]; then
  IFS=":"
  for var in $par_extension; do
    unset IFS
    extra_params+=( "--extension" "$var" )
  done
fi

bs download run \
  --id "$par_id" \
  --output "$par_output" \
  --api-server "$par_api_server" \
  --access-token "$par_access_token" \
  "${extra_params[@]}"