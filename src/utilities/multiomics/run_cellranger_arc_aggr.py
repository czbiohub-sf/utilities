#!/usr/bin/env python

import csv
import os
import posixpath

from utilities.log_util import get_logger
from utilities.multiomics.common import (
    get_base_parser,
    get_default_requirements,
    prepare_and_return_base_data_paths,
    process_results
)
from utilities.s3_util import s3_cp, s3_sync


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

    data_dir = paths["data_dir"]

    original_libraries_path = paths["data_dir"] / run_id / "original_libraries.csv"
    libraries_path = paths["data_dir"] / run_id / "libraries.csv"

    s3_cp(
        logger,
        args.s3_libraries_csv_path,
        str(original_libraries_path)
    )

    with open(original_libraries_path) as csvfile, \
         open(libraries_path, 'w') as new_csv:
        header = next(csvfile)
        new_csv.write(f"{header}")

        for row in csv.reader(csvfile):
            atac_fragments = row[1]
            per_barcode_metrics = row[2]
            gex_molecule_info = row[3]

            # TODO(neevor): clean up duplicated logic.
            atac_fragments_local_path = data_dir / posixpath.basename(atac_fragments)
            per_barcode_metrics_local_path = data_dir / posixpath.basename(per_barcode_metrics)
            gex_molecule_info_local_path = data_dir / posixpath.basename(gex_molecule_info)

            s3_sync(logger, atac_fragments, str(atac_fragments_local_path))
            s3_sync(logger, per_barcode_metrics, str(per_barcode_metrics_local_path))
            s3_sync(logger, gex_molecule_info, str(gex_molecule_info_local_path))

            row[1] = str(atac_fragments_local_path)
            row[2] = str(per_barcode_metrics_local_path)
            row[3] = str(gex_molecule_info_local_path)

            row_values = ",".join(row)
            new_csv.write(f"{row_values}\n")


    os.chdir(str(paths["local_output_path"]))
    command = [
        CELLRANGER,
        "aggr",
        f"--id={run_id}",
        f"--csv={str(libraries_path)}",
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
