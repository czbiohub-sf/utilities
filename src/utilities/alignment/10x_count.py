#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./10x_count.py
import argparse
import os
import sys
import subprocess
import tarfile

from utilities.log_util import get_logger, log_command

import boto3


CELLRANGER = "cellranger"

S3_RETRY = 5
S3_LOG_DIR = "s3://jamestwebber-logs/10xcount_logs/"
S3_REFERENCE = {"east": "czi-hca", "west": "czbiohub-reference"}

reference_genomes = {
    "homo": "HG38-PLUS",
    "hg38-plus": "HG38-PLUS",
    "mus": "MM10-PLUS",
    "mm10-plus": "MM10-PLUS",
    "mm10-1.2.0": "mm10-1.2.0",
    "mus-premrna": "mm10-1.2.0-premrna",
    "mm10-1.2.0-premrna": "mm10-1.2.0-premrna",
    "hg19-mm10-3.0.0": "hg19-mm10-3.0.0",
    "microcebus": "MicMur3-PLUS",
    "gencode.vM19": "gencode.vM19",
}

deprecated = {
    "homo": "hg38-plus",
    "mus": "mm10-plus",
    "mus-premrna": "mm10-1.2.0-premrna",
}


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=256000, storage=2000, ecr_image="demuxer"
    )


def get_parser():
    parser = argparse.ArgumentParser(
        prog="10x_count.py", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--s3_input_dir", required=True)
    parser.add_argument("--s3_output_dir", required=True)
    parser.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
    )
    parser.add_argument("--cell_count", type=int, default=3000)

    parser.add_argument(
        "--region",
        default="east",
        choices=("east", "west"),
        help=(
            "Region you're running jobs in."
            " Should match the location of"
            " the fastq.gz files"
        ),
    )

    parser.add_argument("--glacier", action="store_true")
    parser.add_argument("--root_dir", default="/mnt")

    return parser


def main(logger):
    parser = get_parser()

    args = parser.parse_args()

    if os.environ.get("AWS_BATCH_JOB_ID"):
        args.root_dir = os.path.join(args.root_dir, os.environ["AWS_BATCH_JOB_ID"])

    # local directories
    sample_id = os.path.basename(args.s3_input_dir)
    result_path = os.path.join(args.root_dir, "data", "hca", sample_id)
    fastq_path = os.path.join(result_path, "fastqs")
    os.makedirs(fastq_path)

    genome_base_dir = os.path.join(args.root_dir, "genome", "cellranger")
    os.makedirs(genome_base_dir)

    if args.taxon in reference_genomes:
        if args.taxon in deprecated:
            logger.warn(
                f"'{args.taxon}' will be removed in the future,"
                f" use '{reference_genomes[args.taxon]}'"
            )

        genome_name = reference_genomes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    if args.region != "west" and genome_name not in ("HG38-PLUS", "MM10-PLUS"):
        raise ValueError(f"you must use --region west for {genome_name}")

    # files that should be uploaded outside of the massive tgz
    # path should be relative to the run folder
    files_to_upload = {
        "outs/raw_gene_bc_matrices_h5.h5": "raw_gene_bc_matrices_h5.h5",
        "outs/web_summary.html": "web_summary.html",
        "outs/metrics_summary.csv": "metrics_summary.csv",
    }

    if args.taxon == "hg19-mm10-3.0.0":
        genome_list = ("hg19", "mm10")
    else:
        genome_list = (genome_name,)

    for gn in genome_list:
        files_to_upload.update(
            {
                f"outs/raw_gene_bc_matrices/{gn}/genes.tsv": f"genes.{gn}.tsv",
                f"outs/raw_gene_bc_matrices/{gn}/barcodes.tsv": f"barcodes.{gn}.tsv",
                f"outs/raw_gene_bc_matrices/{gn}/matrix.mtx": f"matrix.{gn}.mtx",
            }
        )

    s3 = boto3.resource("s3")

    # download the ref genome data
    logger.info(f"Downloading and extracting genome data {genome_name}")

    s3_object = s3.Object(S3_REFERENCE[args.region], f"cellranger/{genome_name}.tgz")

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=genome_base_dir)

    genome_dir = os.path.join(genome_base_dir, genome_name)

    sys.stdout.flush()

    # download the fastq files
    command = [
        "aws",
        "s3",
        "cp",
        "--no-progress",
        "--recursive",
        "--force-glacier-transfer" if args.glacier else "",
        args.s3_input_dir,
        fastq_path,
    ]
    log_command(logger, command, shell=True)

    # Run cellranger
    os.chdir(result_path)
    command = [
        CELLRANGER,
        "count",
        "--localmem=240",
        "--nosecondary",
        "--disable-ui",
        f"--expect-cells={args.cell_count}",
        f"--id={sample_id}",
        f"--fastqs={fastq_path}",
        f"--transcriptome={genome_dir}",
    ]
    log_command(
        logger, command, shell=True, stderr=subprocess.STDOUT, universal_newlines=True
    )

    # Move results(websummary, cell-gene table, tarball) data back to S3
    for file_name, dest_name in files_to_upload.items():
        command = [
            "aws",
            "s3",
            "cp",
            "--quiet",
            os.path.join(result_path, sample_id, file_name),
            os.path.join(args.s3_output_dir, dest_name),
        ]
        for i in range(S3_RETRY):
            try:
                log_command(logger, command, shell=True)
                break
            except subprocess.CalledProcessError:
                logger.info(f"retrying cp {file_name}")
        else:
            raise RuntimeError(f"couldn't sync {file_name}")

    tarball_file = f"{os.path.join(result_path, sample_id)}.tgz"

    command = ["tar", "czf", tarball_file, sample_id]
    log_command(logger, command, shell=True)

    command = ["aws", "s3", "cp", "--quiet", tarball_file, f"{args.s3_output_dir}/"]
    for i in range(S3_RETRY):
        try:
            log_command(logger, command, shell=True)
            break
        except subprocess.CalledProcessError:
            logger.info(f"retrying cp {sample_id}.tgz")
    else:
        raise RuntimeError(f"couldn't sync {sample_id}.tgz")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)

    try:
        main(mainlogger)
    except:
        mainlogger.info("An exception occurred", exc_info=True)
        raise
    finally:
        if log_file:
            log_cmd = f"aws s3 cp --quiet {log_file} {S3_LOG_DIR}"
            mainlogger.info(log_cmd)

            file_handler.close()
            subprocess.check_output(log_cmd, shell=True)
