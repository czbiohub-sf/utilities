#!/usr/bin/env python

import csv
import os
import posixpath

from utilities.log_util import get_logger
from utilities.multiomics.common import (
    get_base_parser,
    get_default_requirements,
    prepare_and_return_base_data_paths,
    process_results,
    run_command,
    sync_results
)
from utilities.s3_util import s3_cp


CELLRANGER = "/bin/cellranger-arc"  # NOTE(neevor): I'm not very happy to have to have the /bin so I would like to be able to get rid of it.
get_default_requirements = get_default_requirements


def get_parser():
    """ Construct and return the ArgumentParser object that parses input command.
    """

    parser = get_base_parser(prog="run_cellranger_arc_count.py",
                             description="Run aggr using cellranger arc")

    parser.add_argument("--neurips", action='store_true')

    return parser


def create_annotated_results(
    feature_matrix_path,
    libraries_csv,
    output_path):

    import scanpy as sc
    import pandas as pd
    import numpy as np

    adata = sc.read_10x_h5(feature_matrix_path, gex_only=False)

    # Load aggr.csv and make 1-indexed
    aggr = pd.read_csv(libraries_csv)
    aggr.index = pd.Index(np.arange(aggr.shape[0])+1)

    # Get site and donor for each GEM well
    sites = {}
    donors = {}
    for ix, lib_id in aggr['library_id'].iteritems():
        *_, site, donor = lib_id.split('_')
        sites[ix] = site
        donors[ix] = donor

    # Get GEM well for each cell barcode
    gem_well = [ix.split('-')[1] for ix in adata.obs.index]

    # Append site and donor information
    adata.obs['site'] = [sites[int(g_well)] for g_well in gem_well]
    adata.obs['donor'] = [donors[int(g_well)] for g_well in gem_well]

    # Write output
    adata.write_h5ad(output_path, compression=9)


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
            library = row[0]
            atac_fragments = row[1]
            atac_fragments_tbi = f"{atac_fragments}.tbi"
            per_barcode_metrics = row[2]
            gex_molecule_info = row[3]

            # TODO(neevor): clean up duplicated logic.
            atac_fragments_local_path = data_dir / library / posixpath.basename(atac_fragments)
            atac_fragments_tbi_local_path = data_dir / library / f"{posixpath.basename(atac_fragments)}.tbi"
            per_barcode_metrics_local_path = data_dir / library / posixpath.basename(per_barcode_metrics)
            gex_molecule_info_local_path = data_dir / library / posixpath.basename(gex_molecule_info)

            s3_cp(logger, atac_fragments, str(atac_fragments_local_path))
            s3_cp(logger, atac_fragments_tbi, str(atac_fragments_tbi_local_path))
            s3_cp(logger, per_barcode_metrics, str(per_barcode_metrics_local_path))
            s3_cp(logger, gex_molecule_info, str(gex_molecule_info_local_path))

            row[1] = str(atac_fragments_local_path)
            row[2] = str(per_barcode_metrics_local_path)
            row[3] = str(gex_molecule_info_local_path)

            row_values = ",".join(row)
            new_csv.write(f"{row_values}\n")

    os.chdir(str(paths["result_path"]))
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

    run_command(logger,
                command,
                "cellranger-arc aggr failed")

    if args.neurips:
        create_annotated_results(
            f"{paths['local_output_path']}/filtered_feature_bc_matrix.h5",
            libraries_path,
            f"{paths['local_output_path']}/filtered_feature_bc_matrix.annotated.h5ad"
        )

    sync_results(logger, paths)


if __name__ == "__main__":
    mainlogger, log_file, file_handler = get_logger(__name__)
    main(mainlogger)
