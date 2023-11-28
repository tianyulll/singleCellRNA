"
Parse cell barcodes files to group by integration sites
Identify LTR cells on PC spaces
Compared distance in PC spaces between different group of LTR-cells
"

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

integrated <- readRDS("filterCount/processedData/fitered_integrated_dr.rds")
result<-read.csv("result_marked.tsv", sep = '\t')

"
# Get unique cell barcodes with duplicated integration site
# duplicatePos <- result[duplicated(result$posid) | duplicated(result$posid, fromLast = T), ] # Finds duplicated integration sites
# dup_bc <- unique(duplicatePos$cellBarcode)
# Rename cell ID in Seurat
# To match ambiguously regardless of samples
# cellName <- lapply(colnames(x = integrated$integrated), function(x){
#   strsplit(x, "-")[[1]][1]
# })
#renamed.int <- RenameCells(integrated,new.names = cellName)
"

integrated@meta.data$ltr <- colnames(integrated) %in% result$labeled
integrated@meta.data$ltr <- ifelse(integrated@meta.data$ltr, "1-LTR", "0-Normal")
table(integrated@meta.data$ltr)

# Check all cells with LTR read through
Idents(integrated)<-integrated@meta.data$ltr
DimPlot(integrated, reduction = "int_umap", cells.highlight = WhichCells(integrated, idents = "1-LTR"), 
        raster = F)+NoLegend()+ggtitle("198 Cells with duplicated LTR readthrough UMAP")
DimPlot(integrated, reduction = "int_pca", cells.highlight = WhichCells(integrated, idents = "1-LTR"), 
        raster = F)+NoLegend()+ggtitle("198 Cells with duplicated LTR readthrough PCA")


# Group LTR-cells by identical integration sites
# Only select ones that has more than 1 cell
grouped <- result %>% 
  distinct(labeled, .keep_all = TRUE) %>%
  group_by(posid) %>%
  summarize(bc_list = list(unique(labeled))) %>%
  filter(lengths(bc_list) > 1)

# Label cells with groups
integrated@meta.data$group <- 0
for (i in 1:dim(grouped)[1]) {
  print(i)
  integrated@meta.data$group <- ifelse(colnames(integrated) %in% unlist(grouped[[2]][i]), i, integrated@meta.data$group) 
}

# Plot labeled cells in reduced dimensions
Idents(integrated) <- integrated@meta.data$group
highlight<-list()
n <- length(table(Idents(integrated)))-1
for (i in 1: n) {
  highlight[[i]] <- WhichCells(integrated, idents = i)
}
colorSet <- brewer.pal(n = n, name = "Paired")

DimPlot(integrated, reduction = "int_umap",group.by = "group", 
        cells.highlight = highlight, cols.highlight = colorSet,
        raster = F)+ggtitle("Grouped by integration sites UMAP")
DimPlot(integrated, reduction = "int_pca",group.by = "group", 
        cells.highlight = highlight, cols.highlight = colorSet,
        raster = F)+ggtitle("Grouped by integration sites PCA")


# Calculate PC distance 
pcTable <- integrated$int_pca@cell.embeddings[, 1:2]

# Calculate Averaged Euclidean distance in PC1-PC2 per cell
calculate_chunk_distances <- function(df, chunk_size) {
  total_distances <- numeric(0)
  num_chunks <- ceiling(nrow(df) / chunk_size)
  
  for (i in 1:num_chunks) {
    print(paste("in iteration", i))
    start_index <- (i - 1) * chunk_size + 1
    end_index <- min(i * chunk_size, nrow(df))
    
    chunk_distances <- dist(df[start_index:end_index, ])
    total_distances <- c(total_distances, rowMeans(as.matrix(chunk_distances)))
  }
  return(total_distances)
}

distList <- calculate_chunk_distances(pcTable, 10000)
plot(density(distList), main = "Density Plot of Average Distances per Cell",
     xlab = "Averaged Euclidean Distance in PC", ylab = "Density")
mean(distL)

euclideanDist <- function(df) {
  return(mean(as.numeric(dist(df))))
}

# Same integration site cells Euclidean distance in PC1-PC2
sameSite <- unlist(grouped[[2]][1:2])
euclideanDist(pcTable[sameSite[1:2],])
euclideanDist(pcTable[sameSite[3:5],])

# Randomly select points for PC distance
randomDist <- function(df, num_rows, time, seed){
  set.seed(seed)
  distList <- list()
  counter <- time/10
  for (i in 1:time){
    if (i %% counter == 0){
      print(paste("in iteration", i))
    }
    ran_df <- df[sample(nrow(df), num_rows, replace = F), ]
    distList<-c(distList, mean(as.matrix(dist(ran_df))))
  }
  return(mean(unlist(distList)))
}

randomDist(pcTable, 10, 100000, seed = 0)




