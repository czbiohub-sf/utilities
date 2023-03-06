import scanpy as sc
import mudata as mu
import logging
from anndata import AnnData
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_key": "pca",
    "num_components": 25,
    "layer": None,
    "obsm_output": "X_pca",
    "var_input": "filter_with_hvg",
    "varm_output": "varm_output",
    "uns_output": "pca_variance",
    "overwrite": True,
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
mdata = mu.read_h5mu(par["input"])

logger.info("Computing PCA components for modality '%s'", par['modality'])
data = mdata.mod[par['modality']]
if par['layer'] and par['layer'] not in data.layers:
    raise ValueError(f"{par['layer']} was not found in modality {par['modality']}.")
layer = data.X if not par['layer'] else data.layers[par['layer']]
adata_input_layer = AnnData(layer)
adata_input_layer.var.index = data.var.index

use_highly_variable = False
if par["var_input"]:
    if not par["var_input"] in data.var.columns:
        raise ValueError(f"Requested to use .var column {par['var_input']} "
                         "as a selection of genes to run the PCA on, "
                         f"but the column is not available for modality {par['modality']}")
    use_highly_variable = True
    adata_input_layer.var['highly_variable'] = data.var[par["var_input"]]

# run pca
output_adata = sc.tl.pca(
    adata_input_layer,
    n_comps=par["num_components"],
    copy=True,
    use_highly_variable=use_highly_variable
)

# store output in specific objects

check_exist_dict = {
    "obsm_output": ("obs"),
    "varm_output": ("varm"),
    "uns_output": ("uns")
}
for parameter_name, field in check_exist_dict.items():
    if par[parameter_name] in getattr(data, field):
        if not par["overwrite"]:
            raise ValueError(f"Requested to create field {par[parameter_name]} in .{field} "
                            f"for modality {par['modality']}, but field already exists.")
        del getattr(data, field)[par[parameter_name]]

data.obsm[par["obsm_output"]] = output_adata.obsm['X_pca']
data.varm[par["varm_output"]] = output_adata.varm['PCs']
data.uns[par["uns_output"]] = { "variance": output_adata.uns['pca']['variance_ratio'],
                                "variance_ratio": output_adata.uns['pca']['variance_ratio'] }


logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"], compression="gzip")

logger.info("Finished")