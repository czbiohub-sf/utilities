# Instructions on Running Pipelines

## Run Published Pipeline

These are instructions on how to launch a pipeline using the published code without cloning and using this repository.

## Run Local Pipeline

These are instructions on running a utitlities pipeline from within the github reposository. 

### CellRanger Preprocessing

Running `cellranger` demultiplexing and alignment.

```
bin/nextflow run ./cellranger.nf \
    --id example_run \
    --input ./resources_test/cellranger_tiny_bcl/bcl \
    --sample_sheet ./resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv \
    --transcriptome ./resources_test/cellranger_tiny_fastq/cellranger_tiny_ref \
    --publishDir ./example_run_output
```