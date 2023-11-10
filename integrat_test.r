library(Seurat) 
library(parallel)
library(dplyr)

setwd("/mnt/isilon/vellalabserver/tianyu/single_cell")

obj<-readRDS("objList_pca.rds")
features <- SelectIntegrationFeatures(obj, nfeatures = 1200)
anchor<-FindIntegrationAnchors(object.list = obj, anchor.features = features,
                               reference = c(1, 2), reduction = "rpca")
saveRDS(anchor, file = "anchor_t.rds")