import muon as mu

## VIASH START
par = {
  'input': 'input.h5',
  'output': 'output.h5mu'
}
## VIASH END

mdata = mu.read_10x_h5(par["input"])
mdata.var_names_make_unique()
mdata.write_h5mu(par["output"])
