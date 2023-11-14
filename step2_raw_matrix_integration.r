'''
R script to parse count matrix
'''
library(Seurat) 
library(parallel)
library(dplyr)
library(reticulate)

# Find count matrix 
h5_files <- list.dirs(path = "scRNAdata")
#h5_files <- h5_files[grepl("*raw_feature*", h5_files)]
h5_files<-h5_files[2:17]

# read in the matrix file
sc_mat <- lapply(h5_files, function(x) {
  Read10X(x)
})

# Create a list of seurat object
sc_obj <- mclapply(sc_mat, mc.cores = 2,function(x) {
  CreateSeuratObject(counts = x) 
})

sampleName <- lapply(h5_files, function(x){
  strsplit(x, "/")[[1]][2]
})

# add labeling
for (i in 1:length(sc_obj)) {
  sc_obj[[i]]@meta.data$sample <- sampleName[[i]]
}

features <- SelectIntegrationFeatures(object.list = sc_obj)
save(features, file="int_feature.rds")

'''
for (i in 1:length(sc_obj)) {
  print(i)
  sc_obj[[i]] <- sc_obj[[i]] %>% 
    NormalizeData() %>%
    FindVariableFeatures(., nfeatures = 1000) %>%
    ScaleData() %>%
    RunPCA()
  sc_obj[[i]]$RNA$scale.data <- NULL
  print("remove cache")
  gc()
  print("clear")
}
'''

# Remove scaled data for memory efficiency
for (i in 1:length(obj)) {
  obj[[i]]$RNA$scale.data <- NULL
}
gc()

# Anchoring for intergation 
obj<-readRDS("objList_pca.rds")
features <- SelectIntegrationFeatures(obj, nfeatures = 1200)
saveRDS(features, file = "feat.rds")
anchor<-FindIntegrationAnchors(object.list = obj, anchor.features = features,
                               reduction = "rpca")
saveRDS(anchor, file = "anchor.rds")

# Integrate dataset 
anchor<-readRDS("anchor_t.rds")
integrated <- IntegrateData(anchor)

# Merge the objects
#merge_obj <- merge(sc_obj[[1]], y = sc_obj[2:16], project = "scRNA")

# Alternative way to integrate datasets
# Run with SCTransform
sc_list <- lapply(h5_files, function(x) {
  x <- SCTransform(x)
  x <- RunPCA(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, normalization.method = "SCT", reduction = "rpca")
ifnb.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Convert to annData
integrated <- readRDS("integrated.rds")
SaveH5Seurat(integrated, filename = "integrated.h5Seurat")
Convert("integrated.h5Seurat", dest = "h5ad")

# Filter extremely low quality cell
# 10427353 > 5710359
integrated <- subset(integrated, subset = nFeature_RNA > 1 & nCount_RNA > 1)

# Scale and PCA
setwd("/mnt/isilon/vellalabserver/tianyu/single_cell")
integrated <- ScaleData(integrated, assay = "integrated")
integrated <- RunPCA(integrated, assay = "integrated", reduction.name = "intPCA", reduction.key = "intPC")
integrated$integrated$scale.data <- NULL
saveRDS(integrated, file="integrated.rds")

'''
Notes for reticulate
use_condaenv("basic")
reticulate::py_install(packages ="umap-learn")
reticulate::use("r-reticulate")
'''
#library(future)
#future::plan("multisession", workers = 8)
#integrated <- RunUMAP(integrated, dims = 1:20, reduction = "intPCA", umap.method="umap-learn")
integrated <- RunUMAP(object=integrated, dims = 1:20, ,reduction = "intPCA", re)
saveRDS(integrated, file="integrated_umap.rds")

