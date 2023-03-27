#!/bin/bash

set -e

# settings
ID=reference_gencodev41_chr1
OUT=resources_test/$ID

mkdir -p $OUT

NXF_VER=21.10.6 nextflow \
  run openpipelines-bio/openpipeline \
  -r "0.6.2" \
  -main-script workflows/ingestion/make_reference/main.nf \
  -profile docker \
  --id "$ID" \
  --genome_fasta "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz" \
  --transcriptome_gtf "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" \
  --target "cellranger:star" \
  --output_fasta "reference.fa.gz" \
  --output_gtf "reference.gtf.gz" \
  --output_cellranger "reference_cellranger.tar.gz" \
  --output_star "reference_star" \
  --subset_regex "chr1" \
  --publish_dir $OUT