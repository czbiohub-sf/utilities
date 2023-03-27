## utilities scripts

These are helper scripts used to run pipelines on the CZ Biohub HPC.

## Procedure

1. **Create a new tmux session**. Example: `tmux new -s mysession`.

2. **Go into a working directory**. Example: `cd /hpc/scratch/group.data.science/nextflow_$USER/mysession`. First create this directory if it does not yet exist.

3. **Launch a pipeline**. The `/hpc/scratch/group.data.science/utilities/scripts/` directory contains helper executables for launching pipeline scripts. 
  Run `/hpc/scratch/group.data.science/utilities/scripts/<pipeline> --help` for more information on the pipeline's parameters.

4. (Optional) **Leave the tmux session**. Press `Control+B` and then `D` to leave the tmux session without killing it.

5. (Optional) **Rejoin the tmux session**. Run `tmux attach -t mysession` to rejoin the session.

6. **Clean up**. After everything finished successfully, remove the original working directory, as the `work/` directory is quite large.

## Example commands

### Process 10x data

```bash
# create a new tmux session
tmux new -s tsp27_10x

# change into WDIR, create it if it does not exist yet
WDIR=/hpc/scratch/group.data.science/nextflow_$USER/tsp27_10x
[ -d "$WDIR" ] || mkdir -p "$WDIR"
cd "$WDIR"

# launch pipeline
/hpc/scratch/group.data.science/utilities/scripts/process_10x \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot27/fastqs/10X/" \
  --reference "/hpc/reference/sequencing_alignment/alignment_references/human_gencode_v41_ercc_cellranger.tgz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP27/mapping/10X/"
```

### Process Smartseq2 data

```bash
# create a new tmux session
tmux new -s tsp2_batch1_ss2

# change into WDIR, create it if it does not exist yet
WDIR=/hpc/scratch/group.data.science/nextflow_$USER/tsp2_batch1_ss2
[ -d "$WDIR" ] || mkdir -p "$WDIR"
cd "$WDIR"

# launch pipeline
/hpc/scratch/group.data.science/utilities/scripts/process_smartseq2 \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot2/fastqs/smartseq2/batch1" \
  --reference_index "/hpc/reference/sequencing_alignment/alignment_references/human_gencode_v41_ercc_star.tgz" \
  --reference_gtf "/hpc/reference/sequencing_alignment/gff_files/human_gencode_v41.gtf.gz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP2/mapping/smartseq2/batch1"
```

## Help message

Run `/hpc/scratch/group.data.science/utilities/scripts/process_smartseq2 --help` to view a pipeline's help message.

## Stub run

Add `-stub` to one of the commands to check whether a pipeline is likely to run before actually running it.

## Additional information

The helper scripts were generated using:

```bash
viash run src/operations/create_runner_script/config.vsh.yaml -- \
  --output scripts/process_10x \
  --repository czbiohub/utilities \
  --main_script src/mapping/process_10x/main.nf \
  --entry auto \
  --revision main_build \
  --singularity_cache_dir /hpc/scratch/group.data.science/singularity_images \
  --tower_workspace_id 124983787544228 \
  --nextflow_config /hpc/scratch/group.data.science/utilities/scripts/czbhpc.config

viash run src/operations/create_runner_script/config.vsh.yaml -- \
  --output scripts/process_smartseq2 \
  --repository czbiohub/utilities \
  --main_script src/mapping/process_smartseq2/main.nf \
  --entry auto \
  --revision main_build \
  --singularity_cache_dir /hpc/scratch/group.data.science/singularity_images \
  --tower_workspace_id 124983787544228 \
  --nextflow_config /hpc/scratch/group.data.science/utilities/scripts/czbhpc.config
```