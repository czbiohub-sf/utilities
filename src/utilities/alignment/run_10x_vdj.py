#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./10x_count.py
import argparse
import os
import pathlib
import sys
import subprocess
import tarfile
import posixpath

from utilities.log_util import get_logger, log_command

import boto3


# reference genome bucket name for different regions
S3_REFERENCE = {"east": "czbiohub-reference-east", "west": "czbiohub-reference"}

# valid and deprecated reference genomes
reference_genomes = {
    "GRCh38-VDJ":"refdata-cellranger-vdj-GRCh38-alts-ensembl-4.0.0",
}

# other helpful constants
CELLRANGER = "cellranger"
S3_RETRY = 5


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=256000, storage=2000, ecr_image="demuxer"
    )


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = argparse.ArgumentParser(
        prog="run_10x_vdj.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run alignment jobs using 10x",
    )

    # required arguments
    requiredNamed = parser.add_argument_group("required arguments")

    requiredNamed.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
        help="Reference genome for the alignment run",
    )

    requiredNamed.add_argument(
        "--s3_input_path", required=True, help="The folder with fastq.gz files to align"
    )

    requiredNamed.add_argument(
        "--s3_output_path",
        required=True,
        help="The folder to store the alignment results",
    )

    requiredNamed.add_argument(
        "--num_partitions",
        type=int,
        required=True,
        default=10,
        help="Number of groups to divide samples "
        "into for the alignment run. Enter 10 as the default "
        "value here since we don't divide a single sample",
    )

    requiredNamed.add_argument(
        "--partition_id",
        type=int,
        required=True,
        help="Index of sample group. Enter 0 as "
        "the default value here since we only have one sample",
    )

    # optional arguments
    parser.add_argument("--cell_count", type=int, default=3000)

    parser.add_argument(
        "--dobby",
        action="store_true",
        help="Use if 10x run was demuxed locally (post November 2019)",
    )

    parser.add_argument(
        "--region",
        default="west",
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
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()

    args = parser.parse_args()

    args.root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        args.root_dir = args.root_dir / os.environ["AWS_BATCH_JOB_ID"]

    # local directories
    if args.s3_input_path.endswith("/"):
        args.s3_input_path = args.s3_input_path[:-1]

    sample_id = os.path.basename(args.s3_input_path)
    result_path = args.root_dir / "data" / sample_id
    if args.dobby:
        fastq_path = result_path
    else:
        fastq_path = result_path / "fastqs"
    fastq_path.mkdir(parents=True)

    genome_base_dir = args.root_dir / "genome" / "cellranger"
    genome_base_dir.mkdir(parents=True)

    # check if the input genome and region are valid
    if args.taxon in reference_genomes:
        genome_name = reference_genomes[args.taxon]
    else:
        raise ValueError(f"unknown taxon {args.taxon}")

    genome_dir = genome_base_dir / genome_name
    ref_genome_10x_file = f"cellranger/{genome_name}.tgz"

    if args.region != "west" and genome_name not in ("HG38-PLUS", "MM10-PLUS"):
        raise ValueError(f"you must use --region west for {genome_name}")

    if args.region == "east":
        ref_genome_10x_file = f"ref-genome/{ref_genome_10x_file}"

    logger.info(
        f"""Run Info: partition {args.partition_id} out of {args.num_partitions}
                   genome_dir:\t{genome_dir}
         ref_genome_10x_file:\t{ref_genome_10x_file}
                        taxon:\t{args.taxon}
                s3_input_path:\t{args.s3_input_path}"""
    )

    s3 = boto3.resource("s3")

    # download the reference genome data
    logger.info(f"Downloading and extracting genome data {genome_name}")

    s3_object = s3.Object(S3_REFERENCE[args.region], ref_genome_10x_file)

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=genome_base_dir)

    sys.stdout.flush()

    # download the fastq files
    command = [
        "aws",
        "s3",
        "cp",
        "--no-progress",
        "--recursive",
        "--force-glacier-transfer" if args.glacier else "",
        args.s3_input_path,
        f"{fastq_path}",
    ]
    log_command(logger, command, shell=True)

    logger.info(f"Running partition {args.partition_id} of {args.num_partitions}")

    # check the input folder for existing runs
    sample_name = {
        os.path.basename(fn).rsplit("_", 4)[0] for fn in fastq_path.glob("*fastq.gz")
    }
    assert len(sample_name) == 1, "Should only have one sample name to process"
    sample_name = sample_name.pop()

    # Run cellranger
    os.chdir(result_path)
    command = [
        CELLRANGER,
        "vdj",
        "--localmem=240",
        f"--id={sample_id}",
        f"--fastqs={fastq_path}",
        f"--reference={genome_dir}",
    ]
    if args.dobby:
        command.append(f"--sample={sample_name}")

    failed = log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    if failed:
        raise RuntimeError("cellranger vdj failed")

    # Move outs folder to S3
    command = [
        "aws",
        "s3",
        "sync",
        "--no-progress",
        os.path.join(result_path, sample_id, "outs"),
        posixpath.join(args.s3_output_path, sample_name.rsplit("_", 2)[0]),
    ]
    for i in range(S3_RETRY):
        if not log_command(logger, command, shell=True):
            break
        logger.info(f"retrying sync")
    else:
        raise RuntimeError(f"couldn't sync output")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)

    main(mainlogger)
