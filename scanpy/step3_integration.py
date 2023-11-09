import scanpy as sc 
import anndata as an
import pandas as pd
import numpy as np
import scanorama
import os
'''
Scanorama to integrate/batch correct
This takes ~100G memory
More efficient ways?
'''
adatas = sc.read_h5ad('merge_pca.h5ad')
print("data loaded")
adatas_int = scanorama.integrate_scanpy(adatas, batch_size=100000, dimred=50)
print("intergation done")
adatas_int.write('integrated.h5ad', compression="gzip")
print("file saved")