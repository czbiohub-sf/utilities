#!/usr/bin/env python

import os

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
                           description="Run counts using cellranger arc")


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()
    args = parser.parse_args()

    run_id = args.run_id
    library = args.s3_libraries_csv_path

    paths = prepare_and_return_base_data_paths(run_id, args, logger)

    original_libraries_path = paths["data_dir"] / "original_libraries.csv"
    libraries_path = paths["data_dir"] / "libraries.csv"

    s3_cp(
        logger,
        library,
        str(original_libraries_path)
    )

    process_libraries_file(
        original_libraries_path,
        libraries_path,
        paths["data_dir"],
        logger
    )

    os.chdir(str(paths["result_path"]))
    command = [
        CELLRANGER,
        "count",
        f"--id={run_id}",
        f"--reference={paths['ref_path']}",
        f"--libraries={libraries_path}",
        "--localmem=256",
        "--localcores=64",
    ]

    process_results(logger,
                    command,
                    paths,
                    "cellranger-arc count failed")


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)
    main(mainlogger)
