# from https://github.com/theislab/scvelo/issues/255#issuecomment-739995301
import scanpy as sc
import pandas as pd

def fix_cellxgene(file, file_out=None):
    file_out = file if file_out is None else file_out
    scrna_adata = sc.read_h5ad(file)
    # find all columns containing the word "cluster"
    cluster_cols = scrna_adata.obs.filter(regex='cluster').columns
    # convert all cluster metadata columns to categorical variables
    scrna_adata.obs[cluster_cols] = scrna_adata.obs[cluster_cols].astype('category')
    # for reading into cellxgene
    scrna_adata.var_names = scrna_adata.var.features
    scrna_adata.__dict__['_raw'].__dict__['_var'] = scrna_adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    # del(scrna_adata.var['_index'])
    scrna_adata.write_h5ad(file_out)
