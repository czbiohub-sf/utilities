#!/usr/bin/env python

import argparse
import os
import pathlib
import subprocess
import posixpath

from utilities.log_util import get_logger, log_command
from utilities.multiomics.common import (
    get_base_parser,
    get_default_requirements,
    process_libraries_file
)
from utilities.references import (
    download_cellranger_reference,
    reference_genomes
)
from utilities.s3_util import s3_cp, s3_sync


CELLRANGER = "/bin/cellranger"  # NOTE(neevor): I'm not very happy to have to have the /bin so I would like to be able to get rid of it.
get_default_requirements = get_default_requirements


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = get_base_parser(prog="run_cellranger_cite_seq.py",
                             description="Run cite-seq using cellranger count")

    parser.add_argument(
        "--s3_feature_ref_path",
        required=True,
        help="The s3 path of the feature reference csv",
    )

    parser.add_argument(
        "--expect_cells",
        default=3000,
        help="Expected number of cells",
    )

    return parser


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()
    args = parser.parse_args()

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

    feature_ref = result_path / "feature_ref.csv"

    s3_cp(logger, args.s3_libraries_csv_path, str(original_libraries_path))
    s3_cp(logger, args.s3_feature_ref_path, str(feature_ref))

    ref_path = download_cellranger_reference(args.taxon, genome_dir, logger)

    process_libraries_file(original_libraries_path, libraries_path, data_dir, logger)

    # Run cellranger
    os.chdir(result_path)
    command = [
        CELLRANGER,
        "count",
        f"--id={args.run_id}",
        f"--transcriptome={ref_path}",
        f"--libraries={libraries_path}",
        f"--feature-ref={feature_ref}",
        f"--expect-cells={args.expect_cells}",
        "--localmem=256",
        "--localcores=64",
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

    local_output_path = posixpath.join(result_path, args.run_id, "outs")
    s3_output_path = posixpath.join(args.s3_output_path, args.run_id)

    s3_sync(logger, local_output_path, s3_output_path)


if __name__ == "__main__":
    print("Starting job processing")
    mainlogger, log_file, file_handler = get_logger(__name__)
    main(mainlogger)
