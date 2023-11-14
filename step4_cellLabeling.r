'''
Parse cell barcodes files to group by integration sites
Identify LTR cells on PC spaces
Compared distance in PC spaces between different group of LTR-cells
'''

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

integrated <- readRDS("filterCount/processedData/fitered_integrated.rds")

'''
bc <- unique(u$cellBarcode)
bc <- lapply(bc, function(x){
  paste(x, "-1_1", sep = "")
})
'''

# Get unique cell barcodes with duplicated integration site
result<-read.csv("result.csv")
duplicatePos <- result[duplicated(result$posid) | duplicated(result$posid, fromLast = T), ] # Finds duplicated integration sites
dup_bc <- unique(duplicatePos$cellBarcode)

# Rename cell ID in Seurat
# To match ambiguously regardless of samples
cellName <- lapply(colnames(x = integrated$integrated), function(x){
  strsplit(x, "-")[[1]][1]
})
#renamed.int <- RenameCells(integrated,new.names = cellName)
integrated@meta.data$ltr <- cellName %in% dup_bc
#integrated@meta.data$size <- ifelse(integrated@meta.data$ltr, 10, 0.1)
integrated@meta.data$ltr <- ifelse(integrated@meta.data$ltr, "1-LTR", "0-Normal")
#DimPlot(integrated, reduction = "int_umap", group.by = "ltr", pt.size =integrated@meta.data$size, raster = F)
Idents(integrated)<-integrated@meta.data$ltr
DimPlot(integrated, reduction = "int_umap", cells.highlight = WhichCells(integrated, idents = "1-LTR"), 
        raster = F)+NoLegend()+ggtitle("755 Cells with duplicated LTR readthrough UMAP")
DimPlot(integrated, reduction = "int_pca", cells.highlight = WhichCells(integrated, idents = "1-LTR"), 
        raster = F)+NoLegend()+ggtitle("755 Cells with duplicated LTR readthrough PCA")

table(integrated@meta.data$ltr)

# Group LTR-cells by identical integration sites
# Only select ones that has more than 1 cell
grouped <- result %>% 
  distinct(cellBarcode, .keep_all = TRUE) %>%
  group_by(posid) %>%
  summarize(bc_list = list(unique(cellBarcode))) %>%
  filter(lengths(bc_list) > 1)

# Label cells with groups
for (i in 1:dim(grouped)[1]) {
  print(i)
  integrated@meta.data$group <- ifelse(cellName %in% unlist(grouped[[2]][i]), i, integrated@meta.data$group ) 
}
integrated@meta.data <- integrated@meta.data %>%
  mutate(group = ifelse(is.na(group), 0, group))

highlight<-list()
for (i in 1:11) {
  highlight[[i]] <- WhichCells(integrated, idents = i)
}

colorSet <- brewer.pal(n =11, name = "Paired")
DimPlot(integrated, reduction = "int_umap",group.by = "group", 
        cells.highlight = highlight, cols.highlight = colorSet,
        raster = F)+ggtitle("Grouped by integration sites UMAP")

DimPlot(integrated, reduction = "int_pca",group.by = "group", 
        cells.highlight = highlight, cols.highlight = colorSet,
        raster = F)+ggtitle("Grouped by integration sites PCA")
