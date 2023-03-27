# czb-util

A collection of scripts for common data management and processing tasks.

The provided pipelines are built using the [Viash framework](http://www.viash.io) on top of the 
Nextflow workflow system. For more information on Nextflow please visit the [Nextflow github page](https://github.com/nextflow-io/nextflow) 
and the [Nextflow read the docs page](https://www.nextflow.io/docs/latest/index.html).

## Requirements

- Install `Docker`, `Podman` or `Singularity`.
- Install `Nextflow`.

## Run versioned pipeline

These are instructions on how to launch a pipeline using the published code without cloning and using this repository.

### Running Cell Ranger Demux

```sh
nextflow \
  run https://github.com/czbiohub/utilities \
  -r 1.0.0 \
  -main-script workflows/1_ingestion/cellranger_demux/main.nf \
  -resume \
  -latest \
  -with-docker \
  --id tiny_bcl \
  --input resources_test/cellranger_tiny_bcl/bcl \
  --sample_sheet resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv \
  --publishDir temp
```

### Running Cell Ranger Mapping

```sh
nextflow \
  run https://github.com/czbiohub/utilities \
  -r 1.0.0 \
  -main-script workflows/1_ingestion/cellranger_mapping/main.nf \
  -resume \
  -latest \
  -with-docker \
  --id tiny_fastq \
  --input resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq \
  --reference resources_test/cellranger_tiny_fastq/cellranger_tiny_ref \
  --publishDir temp
```

## Setting up development environment

Install Viash and Nextflow. 

Run `viash ns build --parallel --setup cb` to build all components and Docker containers locally.

Run `viash ns test --parallel` to unit test all components.
