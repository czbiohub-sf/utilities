#!/usr/bin/env python

import argparse
import csv
import os
import pathlib
import subprocess
import posixpath

from utilities.log_util import get_logger, log_command
from utilities.references import (
    download_cellranger_reference,
    reference_genomes
)


# other helpful constants
CELLRANGER = "cellranger-arc"
S3_RETRY = 5


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=256000, storage=2000, ecr_image="multiomics"
    )


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = argparse.ArgumentParser(
        prog="run_cellranger_arc_count.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Run counts using cellranger arc",
    )

    parser.add_argument(
        "--taxon",
        required=True,
        choices=list(reference_genomes.keys()),
        help="Reference genome for the alignment run",
    )

    parser.add_argument(
        "--run_id",
        required=True,
        help="Name of the folder to write results to"
    )

    parser.add_argument(
        "--s3_libraries_csv_path",
        required=True,
        help="The csv with the s3 paths and metadata needed for cellranger arc count"
    )

    parser.add_argument(
        "--s3_output_path",
        required=True,
        help="The folder to store the alignment results",
    )

    parser.add_argument("--root_dir", default="/mnt")

    return parser


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()
    args = parser.parse_args()

    args.root_dir = pathlib.Path(args.root_dir)

    run_id = args.run_id

    if os.environ.get("AWS_BATCH_JOB_ID"):
        args.root_dir = args.root_dir / os.environ["AWS_BATCH_JOB_ID"]

    data_dir = args.root_dir / "data"
    data_dir.mkdir(parents=True)

    result_path = args.root_dir / "data" / run_id
    result_path.mkdir(parents=True)

    orignial_libraries_path = data_dir / "original_libraries.csv"
    libraries_path = data_dir / "libraries.csv"

    genome_dir = args.root_dir / "genome" / "reference"
    genome_dir.mkdir(parents=True)

    s3_cp(logger, args.s3_libraries_csv_path, str(orignial_libraries_path))

    ref_path = download_cellranger_reference(args.taxon, genome_dir, logger)

    with open(orignial_libraries_path, newline='') as csvfile, \
         open(libraries_path, 'w') as new_csv:
        headers = next(csvfile)
        new_csv.write(f"{headers}")

        for row in csv.reader(csvfile):
            s3_path_of_fastqs = row[0]
            sample_id = row[1]
            method = row[-1].replace(" ", "_")
            local_path = data_dir / posixpath.basename(s3_path_of_fastqs) / sample_id / method
            s3_sync(logger, s3_path_of_fastqs, str(local_path))
            row[0] = str(local_path)
            row_values = ",".join(row)
            new_csv.write(f"{row_values}\n")

    # Run cellranger
    os.chdir(result_path)
    command = [
        CELLRANGER,
        "count",
        f"--id={run_id}",
        f"--reference={ref_path}",
        f"--libraries={libraries_path}",
        "--localmem=64",
        "--localcores=16",
    ]

    failed = log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    if failed:
        raise RuntimeError("cellranger-arc count failed")

    local_output_path = posixpath.join(result_path, run_id, "outs")
    s3_output_path = posixpath.join(args.s3_output_path, run_id)

    s3_sync(logger, local_output_path, s3_output_path)


def s3_sync(logger, input, output):
    command = [
        "aws",
        "s3",
        "sync",
        "--no-progress",
        input,
        output,
    ]
    for _ in range(S3_RETRY):
        if not log_command(logger, command, shell=True):
            break
        logger.info(f"retrying sync")
    else:
        raise RuntimeError(f"couldn't sync output")

    
def s3_cp(logger, input, output):
    command = [
        "aws",
        "s3",
        "cp",
        input,
        output,
    ]
    for _ in range(S3_RETRY):
        if not log_command(logger, command, shell=True):
            break
        logger.info(f"retrying cp")
    else:
        raise RuntimeError(f"couldn't cp output")
    

if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)
    main(mainlogger)
