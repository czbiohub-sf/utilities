#!/usr/bin/env python

import os
import posixpath

import pandas as pd

from utilities.log_util import get_logger
from utilities.multiomics.common import (
    get_base_parser,
    get_default_requirements,
    prepare_and_return_base_data_paths,
    process_libraries_file,
    process_results
)
from utilities.s3_util import s3_cp


CELLRANGER = "/bin/cellranger-arc"  # NOTE(neevor): I'm not very happy to have to have the /bin so I would like to be able to get rid of it.
get_default_requirements = get_default_requirements


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    return get_base_parser(prog="run_cellranger_arc_count.py",
                           description="Run aggr using cellranger arc")


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()
    args = parser.parse_args()

    run_id = args.run_id

    paths = prepare_and_return_base_data_paths(run_id, args, logger)

    os.chdir(str(paths["local_output_path"]))
    command = [
        CELLRANGER,
        "aggr",
        f"--id={run_id}",
        f"--csv={args.s3_libraries_csv_path}",
        f"--reference={paths['ref_path']}",
        "--normalize=depth",
        "--localmem=256",
        "--localcores=64",
    ]

    process_results(logger,
                    command,
                    paths,
                    "cellranger-arc aggr failed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)
    main(mainlogger)
