import scanpy as sc 
import anndata as an
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scanorama
import os

# Get file paths
raw_counts_file=[x[0] for x in os.walk("/data2/rawCount")][1:]
raw_counts_file

# Read in raw matrix
raw_count = []
for i, dir in enumerate(raw_counts_file):
    print(i)
    raw_count.append(sc.read_10x_mtx(dir, cache=True))

# Preprocess annData
for i, adata in enumerate(raw_count):
    print(i)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, inplace=True)

# Read in Sample name
sample = []
for i in raw_counts_file:
    sample.append(i.split("/")[3])

# Add sample ID to datasets
for i, x in enumerate(raw_count):
    x.obs["sample"] = sample[i]

# Preprocess with PCA
sc.pp.pca(adatas)

# Save data
adatas.write('merge_pca.h5ad', compression="gzip")
