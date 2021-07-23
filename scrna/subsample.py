import numpy as np

def q(x):
    return npaarray(list(x))

def subsample(adata,clname1,clname2=None,nc=10000,MIN=15):
    """
    This function goes through each cell type label in `clname1` and then makes sure to select a proportional number
    of each type of cell from `clname2` within that cluster.
    
    So, for example, let's say clname1 = tissue annotations, and clname2 = subtype annotations. 
    This function would go tissue by tissue, and sample a proportional number of cells from each present subtype
    annotation present within that tissue. This ensures all cell types are included when downsampling.
    
    If `clname1=clname2` then it just samples a proportional number of cells from each cluster.
        
    Parameters
    -----------
    adata : AnnData object
    
    clname1 : str
        Annotation key 1

    clname2 : str, optional, default None
        Annotation key 2. If None, `clname2` is set equal to `clname1`
        
    nc : int, optional, default 10000
        Number of cells to target for downsampling. The final number of cells may be slightly lower.
        Tune `nc` accordingly if you'd like to hit an exact target.
        
    MIN : int, optional, default 15
        Minimum cluster size after downsampling
    
    Returns
    -------
    AnnData - A new subsetted AnnData object.
    """
    if clname2 is None:
        clname2 = clname1
    
    frac = nc/adata.shape[0]
    cu = q(adata.obs[clname1])
    cuu = np.unique(cu)
    CELLS=[]
    obsn = q(adata.obs_names)
    for i in range(cuu.size):
        print(i)
        ix = np.where(cu==cuu[i])[0]
        a1 = adata[ix,:]
        lc = q(a1.obs[clname2])
        lcu,lcuc = np.unique(lc,return_counts=True)
        z=0
        cells=[]
        for j in range(lcu.size):
            CELLS.append(q(obsn[ix[np.where(lc==lcu[j])[0][np.random.choice(lcuc[j],replace=False,size=int(lcuc[j]*frac))]]]))
            cells.append(CELLS[-1])
        z = np.concatenate(cells).size
        obsnc = obsn[ix]

        if min(MIN,obsnc.size) - z > 0:
            xx = obsnc[np.in1d(obsnc,np.concatenate(cells),invert=True)]
            CELLS.append(np.random.choice(xx,replace=False,size = min(xx.size,min(30,obsnc.size) - z)))
            print(len(CELLS[-1]))
    CELLS=np.concatenate(CELLS)
    assert np.unique(CELLS,return_counts=True)[1].max()==1
    CELLS = q(adata.obs_names)[np.in1d(q(adata.obs_names),CELLS)]
    return adata[CELLS].copy()