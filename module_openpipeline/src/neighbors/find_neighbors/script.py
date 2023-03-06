import mudata as mu
import scanpy as sc
import logging
from sys import stdout

## VIASH START
par = {
    "input": "work/09/6b10f377b0c86a9da1024f8b9140c0/pbmc_1k_protein_v3_mms.harmonypy.output",
    "output": "output.h5mu",
    "metric": 'cosine',
    "num_neighbors": 15,
    "modality": "rna",
    "obsm_input": "X_pca",
    "uns_output": "neighbors",
    "obsp_distances": "distances",
    "obsp_connectivities": "connectivities"
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Computing a neighborhood graph on modality %s", mod)
adata = mdata.mod[mod]
neighbors = sc.Neighbors(adata)
neighbors.compute_neighbors(
    n_neighbors=par["num_neighbors"],
    use_rep=par["obsm_input"],
    metric=par["metric"],
    random_state=par["seed"],
    method="umap"
)

adata.uns[par["uns_output"]] = {
    'connectivities_key': par["obsp_connectivities"],
    'distances_key': par["obsp_distances"],
    'params': {
        'n_neighbors': neighbors.n_neighbors,
        'method': "umap",
        'random_state': par["seed"],
        'metric': par["metric"],
        'use_rep': par["obsm_input"]
    }
}

adata.obsp[par["obsp_distances"]] = neighbors.distances
adata.obsp[par["obsp_connectivities"]] = neighbors.connectivities

logger.info("Writing to %s", par["output"])
mdata.write_h5mu(filename=par["output"], compression="gzip")
