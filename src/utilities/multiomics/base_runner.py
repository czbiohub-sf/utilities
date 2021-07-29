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


def run(cellranger):
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
