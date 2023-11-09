'''
R script with Seurat
to parse count matrix
and integrate/batch correction 
'''
library(Seurat)
library(parallel)

# Find count matrix 
h5_files <- list.dirs(path = "/data2/rawCount")
#h5_files <- h5_files[grepl("*raw_feature*", h5_files)]
h5_files<-h5_files[2:17]

# read in the matrix file
sc_mat <- lapply(h5_files, function(x) {
  Read10X(x)
})

# Create a list of seurat object
#options(Seurat.object.assay.version = "v5")
sc_obj <- mclapply(sc_mat, mc.cores = 3, function(x) {
  CreateSeuratObject(counts = x) 
})

# add labeling
for (i in 1:length(sc_obj)) {
  sc_obj[[i]]@meta.data$sample <- h5_files[[i]]
}

