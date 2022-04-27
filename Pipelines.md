# Instructions on Running Pipelines

## Run Published Pipeline

These are instructions on how to launch a pipeline using the published code without cloning and using this repository.

## Run Local Pipeline

These are instructions on running a utitlities pipeline from within the github reposository. 

### CellRanger demultiplexing

Running `cellranger` demultiplexing.

```
NXF_VER=21.10.6 nextflow run . \
    -main-script workflows/demux_cellranger/main.nf \
    --id example_run \
    --input ./resources_test/cellranger_tiny_bcl/bcl \
    --sample_sheet ./resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv \
    --output ./output/demux
```

### CellRanger alignment

Running `cellranger` alignment.

```
NXF_VER=21.10.6 nextflow run . \
    -main-script workflows/mapping_cellranger/main.nf \
    --id example_run \
    --input resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq \
    --reference ./resources_test/cellranger_tiny_fastq/cellranger_tiny_ref \
    --output ./output/alignment
```
