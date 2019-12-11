#!/usr/bin/env python

# Example: TAXON=homo CELL_COUNT=3000 S3_DIR=s3://biohub-spyros/data/10X_data/CK_Healthy/ ./10x_count.py
import argparse
import os
import pathlib
import sys
import subprocess
import tarfile

from utilities.log_util import get_logger, log_command

import boto3


CELLRANGER = "cellranger"

S3_RETRY = 5

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
    "zebrafish-plus": "danio_rerio_plus_STAR2.6.1d"
}

deprecated = {"homo", "mus", "mus-premrna"}


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
        "--taxon", required=True, choices=list(reference_genomes.keys())
    )
    parser.add_argument("--cell_count", type=int, default=3000)

    parser.add_argument("--dobby", action="store_true",
                        help="Use if 10x run was demuxed locally (post November 2019)")

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
    parser = get_parser()

    args = parser.parse_args()

    args.root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        args.root_dir = args.root_dir / os.environ["AWS_BATCH_JOB_ID"]

    # local directories
    if args.s3_input_dir.endswith("/"):
        args.s3_input_dir = args.s3_input_dir[:-1]

    sample_id = os.path.basename(args.s3_input_dir)
    result_path = args.root_dir / "data" / sample_id
    if args.dobby:
        fastq_path = result_path
    else:
        fastq_path = result_path / "fastqs"
    fastq_path.mkdir(parents=True)

    genome_base_dir = args.root_dir / "genome" / "cellranger"
    genome_base_dir.mkdir()

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

    s3 = boto3.resource("s3")

    # download the ref genome data
    logger.info(f"Downloading and extracting genome data {genome_name}")

    if args.region == "east":
        s3_object = s3.Object(
            S3_REFERENCE[args.region], f"ref-genome/cellranger/{genome_name}.tgz"
        )
    else:
        s3_object = s3.Object(
            S3_REFERENCE[args.region], f"cellranger/{genome_name}.tgz"
        )

    with tarfile.open(fileobj=s3_object.get()["Body"], mode="r|gz") as tf:
        tf.extractall(path=genome_base_dir)

    genome_dir = genome_base_dir / genome_name

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

    sample_name = {
        os.path.basename(fn).rsplit("_", 3)[0] for fn in fastq_path.glob("*fastq.gz")
    }
    assert len(sample_name) == 1
    sample_name = sample_name.pop()

    # Run cellranger
    os.chdir(result_path)
    command = [
        CELLRANGER,
        "count",
        "--localmem=240",
        "--nosecondary",
        "--disable-ui",
        f"--expect-cells={args.cell_count}",
        f"--fastqs={fastq_path}",
        f"--transcriptome={genome_dir}",
    ]
    if args.dobby:
        command.append(f"--sample={sample_name}")
    else:
        command.append(f"--id={sample_id}")

    failed = log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    if failed:
        raise RuntimeError("cellranger count failed")

    # Move outs folder to S3
    command = [
        "aws",
        "s3",
        "sync",
        "--no-progress",
        os.path.join(result_path, sample_id, "outs"),
        args.s3_output_dir,
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
