'''
Seurat based Workflow
on single cell LTR read through project
tianyu.lu@pennmedicine.upenn.edu
'''
library(Seurat)
library(parallel)
library(dplyr)

# Find count matrix 
h5_files <- list.files(path = "count", recursive = T, pattern = "filtered_feature_bc_matrix.h5")

# Create a copy
for (i in 1:length(h5_files)) {
  splitName <- strsplit(h5_files[[i]], "/")[[1]]
  targetFile<-paste(c("filterCount/data/", splitName[1],"_", splitName[3]), collapse = "")
  sourceFile<-paste(c("count", h5_files[[i]]), collapse = "/")
  file.copy(from = sourceFile, to = targetFile)
}

# read in the matrix file
matrixFile <- list.files(path = "data", full.names = T)
sc_mat <- lapply(matrixFile, function(x) {
  Read10X_h5(x)
})

# Create a list of seurat object
#options(Seurat.object.assay.version = "v5")
sc_obj <- mclapply(sc_mat, mc.cores = 3, function(x) {
  CreateSeuratObject(counts = x) 
})

# label samples in seurat object
sampleName <- lapply(matrixFile, function(x) {
  paste(strsplit(strsplit(x, "/")[[1]][2],"_")[[1]][1:2], collapse = "_")
})

for (i in 1:length(sc_obj)) {
  sc_obj[[i]]@meta.data$sample <- sampleName[[i]]
}

rm(sc_mat)

# Preprocess
for (i in 1:length(sc_obj)) {
  print(i)
  sc_obj[[i]] <- sc_obj[[i]] %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA() 
}

# integrate
anchor <- FindIntegrationAnchors(sc_obj, reduction = "rpca")
integrated <- IntegrateData(anchorset = anchor)

# Dimensional Reduction
integrated <- integrated %>%
  ScaleData() %>%
  RunPCA(reduction.name = "int_pca", reduction.key = "intPC_",) %>%
  RunUMAP(dims = 1:20, reduction = "int_pca", reduction.name = "int_umap")
  
integrated <- RunUMAP(integrated, dims = 1:20, reduction = "int_pca", reduction.name = "int_umap")

# Plot
DimPlot(integrated, reduction = "int_umap")

# Save
saveRDS(integrated, file="processedData/fitered_integrated_dr.rds")

