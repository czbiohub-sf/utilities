## utilities scripts

These are helper scripts used to run pipelines on the CZ Biohub HPC.

Example command to run 10x run:

```bash
tmux new -s nxf_run
cd /hpc/scratch/group.data.science/nextflow_$USER/

/hpc/scratch/group.data.science/utilities/scripts/process_10x \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot27/fastqs/10X/" \
  --reference "/hpc/reference/sequencing_alignment/alignment_references/human_gencode_v41_ercc_cellranger.tgz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP27/mapping/10X/"
```

Example command to perform ss2 run:

```bash
tmux new -s nxf_run
cd /hpc/scratch/group.data.science/nextflow_$USER/

/hpc/scratch/group.data.science/utilities/scripts/process_smartseq2 \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot2/fastqs/smartseq2/" \
  --reference_index "/hpc/reference/sequencing_alignment/alignment_references/human_gencode_v41_ercc_star.tgz" \
  --reference_gtf "/hpc/reference/sequencing_alignment/gff_files/human_gencode_v41.gtf.gz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP2/mapping/smartseq2/"
```

Press `Control+B` and then `D` to leave the tmux session without killing it.

Run `tmux attach -t nxf_run` to rejoin the session.



## Set up

These scripts were generated using:

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