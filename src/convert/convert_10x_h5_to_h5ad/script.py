import scanpy as sc

## VIASH START
par = {
  'input': 'input.h5',
  'output': 'output.h5ad'
}
## VIASH END

adata = sc.read_10x_h5(par["input"])
adata.var_names_make_unique()
adata.write_h5ad(par["output"])
