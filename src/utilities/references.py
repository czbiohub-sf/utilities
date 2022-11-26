import tarfile

import boto3


# reference genome bucket name for different regions
S3_REFERENCE = "czbiohub-reference"

# valid and deprecated reference genomes
reference_genomes = {
    "hg38-plus": "HG38-PLUS",
    "homo.gencode.v30.ERCC.chrM": "homo.gencode.v30.annotation.ERCC92",
    "mm10-plus": "MM10-PLUS",
    "mm10-1.2.0": "mm10-1.2.0",
    "mm10-1.2.0-premrna": "mm10-1.2.0-premrna",
    "hg19-mm10-3.0.0": "hg19-mm10-3.0.0",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
    "GRCh38_premrna": "GRCh38_premrna",
    "zebrafish-plus": "danio_rerio_plus_STAR2.6.1d",
    "botryllus": "botryllus",
    "zebrabow": "Danio.rerio_ZebraBow_genome",
    "GRCh38premrna_and_SARSCoV2": "GRCh38premrna_and_SARSCoV2",
    "arc-GRCh38": "refdata-cellranger-arc-GRCh38-2020-A-2.0.0",
    "arc-mm10": "refdata-cellranger-arc-mm10-2020-A-2.0.0",
    "arc-mm10-mcherry" : "arc-mm10-mcherry",
    "gex-GRCh38": "refdata-gex-GRCh38-2020-A",
    "gex-mm10": "refdata-gex-mm10-2020-A",
    "gex-GRCh38-and-mm10": "refdata-gex-GRCh38-and-mm10-2020-A",
    "mm10_genome_SNseq_Harwell": "mm10_genome_SNseq_Harwell",
    "gencode_mouse_MTB": "gencode_mouse_MTB",
    "zebrafishChromacode" : "Danio.rerio_Chromacode",
    "mouse_genome_mcherry": "mouse_genome_mcherry",
    "SARS.GRCh38_genome" : "SARS.GRCh38_genome",
    "rhesus_human" : "rhesus_human_genome",
    "rhesus_genome" : "rhesus_genome"
}


def validate_taxon(taxon):
    if taxon not in reference_genomes:
        raise ValueError(f"unknown taxon {taxon}")


def download_cellranger_reference(taxon, genome_base_dir, logger):
    """
    Downloads the reference and returns the local path, which can be used in
    calls to cellranger.
    """
    validate_taxon(taxon)
    genome_name = reference_genomes[taxon]

    ref_genome_10x_file = f"cellranger/{genome_name}.tgz"

    s3 = boto3.resource("s3")

    logger.info(f"Downloading and extracting genome data {genome_name}")

    s3_object = s3.Object(S3_REFERENCE, ref_genome_10x_file)

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        
        import os
        
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tf, path=genome_base_dir)

    return genome_base_dir / genome_name
