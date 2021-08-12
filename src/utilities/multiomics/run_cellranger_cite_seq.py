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

    run_id = args.run_id

    paths = prepare_and_return_base_data_paths(run_id,
                                               args,
                                               logger)
    paths["feature_ref"] = paths["result_path"] / "feature_ref.csv"
    paths["original_libraries_path"] = paths["data_dir"] / "original_libraries.csv"
    paths["libraries_path"] = paths["data_dir"] / "libraries.csv"
    paths["sync_to_s3"] = paths["result_path"] / args.run_id / "outs"

    s3_cp(
        logger,
        args.s3_libraries_csv_path,
        str(paths["original_libraries_path"])
    )

    s3_cp(logger, args.s3_feature_ref_path, str(paths["feature_ref"]))

    process_libraries_file(paths["original_libraries_path"],
                           paths["libraries_path"],
                           paths["data_dir"],
                           logger)

    # Run cellranger
    os.chdir(str(paths["local_output_path"]))
    command = [
        CELLRANGER,
        "count",
        f"--id={args.run_id}",
        f"--transcriptome={paths['ref_path']}",
        f"--libraries={paths['libraries_path']}",
        f"--feature-ref={paths['feature_ref']}",
        f"--expect-cells={args.expect_cells}",
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
