import os
import stat

## VIASH START
par = {
  'repository': 'czbiohub/utilities',
  'main_script': 'src/mapping/process_10x/main.nf',
  'entry': 'auto',
  'revision': 'main_build',
  'nextflow_version': '22.04.5',
  'user_group': 'group.data.science',
  'tower_workspace_id': '124983787544228'
}
## VIASH END

# fill in template
str = f"""
#!/bin/bash

UTILITIES_TAG="${{UTILITIES_TAG:-{par['revision']}}}"

set -e

# load nextflow
module load nextflow

# nextflow tower settings
export TOWER_WORKSPACE_ID="{par['tower_workspace_id']}"

# check whether access token is available
if [ -z $TOWER_ACCESS_TOKEN ]; then
  echo Could not find a tower access token!
  exit 1
fi

# nextflow settings
export NXF_VER="{par['nextflow_version']}"

# singularity settings
export NXF_SINGULARITY_CACHEDIR="{par['singularity_cache_dir']}"

# temp settings
export APPTAINER_CACHEDIR="`pwd`/temp"
export APPTAINERENV_TMPDIR="`pwd`/temp"
export NXF_TEMP="`pwd`/temp"
export NXF_WORK="`pwd`/work"

# workaround for submodule being unintentionally removed
# when 'nextflow pull' is not run beforehand.
if [ ! -d "$HOME/.nextflow/assets/czbiohub/utilities/module_openpipeline/" ]; then
  echo "Couldn't find OpenPipelines modules. Recloning the utilities Nextflow pipelines."
  rm -rf "$HOME/.nextflow/assets/czbiohub/utilities/"
  nextflow pull "czbiohub/utilities" -r "$UTILITIES_TAG"
fi

# do run
nextflow run \\
  "{par['repository']}" \\
  -r "$UTILITIES_TAG" \\
  -main-script "{par['main_script']}" \\
  -entry "{par['entry']}" \\
  -resume \\
  -latest \\
  -with-tower \\
  -profile singularity \\
  -c "{par['nextflow_config']}" \\
  $@
"""

# write to file
f = open(par["output"], "w")
f.write(str)
f.close()

# make output file executable
os.chmod(par["output"], 0o775)