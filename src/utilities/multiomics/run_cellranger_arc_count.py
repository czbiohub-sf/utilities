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
                           description="Run counts using cellranger arc")


def main(logger):
    """ Download reference genome, run alignment jobs, and upload results to S3.

        logger - Logger object that exposes the interface the code directly uses
    """

    parser = get_parser()
    args = parser.parse_args()

    run_id = args.run_id

    paths = prepare_and_return_base_data_paths(run_id, args, logger)
    paths["aggr_libraries_path"] = paths['data_dir'] / "aggr_libraries.csv"

    aggr_df = pd.DataFrame({
        "library_id": [],
        "atac_fragments": [],
        "per_barcode_metrics": [],
        "gex_molecule_info": []
    })

    for library in args.s3_libraries_csv_path:
        library_base = os.path.splitext(posixpath.basename(library))[0]
        library_results_path = paths["result_path"] / library_base
        library_results_path.mkdir(parents=True)

        original_libraries_path = paths["data_dir"] / library_base / "original_libraries.csv"
        libraries_path = paths["data_dir"] / library_base / "libraries.csv"

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

        os.chdir(library_results_path)
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

        aggr_df.loc[len(aggr_df.index)] = [
            library_base,
            f"{library_results_path}/{run_id}/outs/atac_fragments.tsv.gz",
            f"{library_results_path}/{run_id}/outs/per_barcode_metrics.csv",
            f"{library_results_path}/{run_id}/outs/gex_molecule_info.h5"
        ]

    aggr_df.to_csv(str(paths["aggr_libraries_path"]), index=False)

    os.chdir(paths["result_path"])
    command = [
        CELLRANGER,
        "aggr",
        f"--id={run_id}",
        f"--csv={str(paths['aggr_libraries_path'])}",
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
