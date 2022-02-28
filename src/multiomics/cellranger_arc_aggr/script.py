#!/usr/bin/env python

import csv
import os
import posixpath
import subprocess

import scanpy as sc
import pandas as pd
import numpy as np

### VIASH START
par = {
    "data": "path/to/data",
    "peaks": "path/to/bedfile",
    "neurips": True,
    "run_id": "1",
    "reference_genome": "path/to/reference_genome",
    "output": "path/to/results",
}
### VIASH END

def create_annotated_results(
        feature_matrix_path,
        libraries_csv,
        output_path):

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

run_id = par["run_id"]

data_dir = par["data_dir"]

original_libraries_path = par["data_dir"] / run_id / "original_libraries.csv"
libraries_path = par["data_dir"] / run_id / "libraries.csv"


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

        row[1] = str(atac_fragments_local_path)
        row[2] = str(per_barcode_metrics_local_path)
        row[3] = str(gex_molecule_info_local_path)

        row_values = ",".join(row)
        new_csv.write(f"{row_values}\n")

os.chdir(str(par["output"]))
command = [
    "cellranger-arc", "aggr",
    f"--id={run_id}",
    f"--csv={str(libraries_path)}",
    f"--reference={paths['ref_path']}",
    "--normalize=depth",
    "--localmem=256",
    "--localcores=64",
]

if par["peaks"]:
    command.append(f"--peaks={par['peaks']}")

# run cellranger-arc aggr
with subprocess.Popen(
    command,
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
) as p:
    for line in p.stdout:
        print(line.decode(), end="")

if p.returncode > 0:
    raise RuntimeError(f"cellranger-arc aggr failed with exit code {p.returncode}")

if par["neurips"]:
    create_annotated_results(
        f"{par['output']}/filtered_feature_bc_matrix.h5",
        libraries_path,
        f"{par['output']}/filtered_feature_bc_matrix.annotated.h5ad"
    )
