## utilities scripts

These are helper scripts used to run pipelines on the CZ Biohub HPC.

Example commands: 

/hpc/scratch/group.data.science/utilities/scripts/process_smartseq2.sh \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot2/fastqs/smartseq2/" \
  --reference_index "/hpc/reference/sequencing_alignment/alignment_references/human_gencode_v41_ercc_star.tgz" \
  --reference_gtf "/hpc/reference/sequencing_alignment/gff_files/human_gencode_v41.gtf.gz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP2/mapping/smartseq2/"
  
  
/hpc/scratch/group.data.science/utilities/scripts/process_10x.sh \
  --input_dir "/hpc/archives/AWS/buckets/tabula-sapiens/Pilot27/fastqs/10X/" \
  --reference_index "/hpc/reference/sequencing_alignment/alignment_references/gencode_v41_ercc_cellranger.tgz" \
  --publish_dir "/hpc/scratch/group.data.science/robrecht_temp/TSP27/mapping/10X/"
  
  
 
