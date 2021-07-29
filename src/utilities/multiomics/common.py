import argparse
import csv
import os
import pathlib
import posixpath
import subprocess

from utilities.log_util import log_command
from utilities.references import (
    download_cellranger_reference,
    reference_genomes
)
from utilities.s3_util import s3_cp, s3_sync


def get_base_parser(prog, description):
    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=description,
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
        help="The s3 path to store the alignment results",
    )

    parser.add_argument("--root_dir", default="/mnt")

    return parser


# TODO(neevor): Clean up the number of args this function takes.
def process_libraries_file(original_libraries_path, libraries_path, data_dir, logger):
    with open(original_libraries_path, newline='') as csvfile, \
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


def get_default_requirements():
    return argparse.Namespace(
        vcpus=64, memory=256000, storage=2000, ecr_image="multiomics"
    )


def prepare_and_return_base_data_paths(run_id, args, logger):
    root_dir = pathlib.Path(args.root_dir)

    if os.environ.get("AWS_BATCH_JOB_ID"):
        root_dir = root_dir / os.environ["AWS_BATCH_JOB_ID"]

    data_dir = root_dir / "data"
    data_dir.mkdir(parents=True)

    result_path = root_dir / "data" / args.run_id
    result_path.mkdir(parents=True)

    original_libraries_path = data_dir / "original_libraries.csv"
    libraries_path = data_dir / "libraries.csv"

    genome_dir = root_dir / "genome" / "reference"
    genome_dir.mkdir(parents=True)

    s3_cp(logger, args.s3_libraries_csv_path, str(original_libraries_path))

    ref_path = download_cellranger_reference(args.taxon, genome_dir, logger)

    return {
        "root_dir": root_dir,
        "data_dir": data_dir,
        "result_path": result_path,
        "original_libraries_path": original_libraries_path,
        "libraries_path": libraries_path,
        "genome_dir": genome_dir,
        "ref_path": ref_path,
        "local_output_path": posixpath.join(result_path, run_id, "outs"),
        "output_dir": posixpath.join(args.s3_output_path, run_id)
    }


def process_results(logger,
                    command,
                    paths,
                    error_message):
    failed = log_command(
        logger,
        command,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,
    )

    if failed:
        raise RuntimeError(error_message)

    s3_sync(logger, paths["local_output_path"], paths["output_dir"])
