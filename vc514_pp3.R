library(umap)
library(Seurat)
library(glmGamPoi)
library(tidyverse)
library(ggplot2)
library(dplyr)

root_dir <- "~/Desktop/Rsession/VisiumDir/"

#root_dir <- "/home/vc514/VisiumDir"
num_objects <- 11
seurat_objects <- list()
for (i in 1:num_objects) {
  inside_folder <- paste0("brain", i, "_spaceranger_out")
  data_dir <- file.path(root_dir, inside_folder)
  h5_mat_name <- "./filtered_feature_bc_matrix.h5"
  slice <- "slice_human"
  seurat_obj <- Seurat::Load10X_Spatial(
    data.dir = data_dir,
    filename = h5_mat_name,
    assay = "Spatial", 
    slice = slice,
    filter.matrix = TRUE, 
    to.upper = FALSE
  )
  seurat_objects[[i]] <- seurat_obj
}

length(seurat_objects)
seurat_objects

# copy spatial assay to RNA
for (i in 1:length(seurat_objects)){
  seurat_objects[[i]]@assays$RNA <- seurat_objects[[i]]@assays$Spatial
}

# Then - typically I would run SCT on each individual

#Created an empty list that stores all the integrated objects from the list seurat_objects
ifnb.list <- list()

#for loop helps in iterating all the objects
for (i in seq_along(seurat_objects)) {
  stim <- SCTransform(seurat_objects[[i]], vst.flavor = "v2",assay ="Spatial",verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)
  object_name <- paste0("object", i)
  ifnb.list[[object_name]] <- stim
}
saveRDS(ifnb.list, file = "brain_1.rds")
ifnb.list <- readRDS("brain_1.rds")


#Select integration features
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
# Prepare for integration
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

# add in anchor based integration
ifnb.anchors <- FindIntegrationAnchors(ifnb.list, normalization.method = "SCT", anchor.features = features)
combined <- IntegrateData(anchorset = ifnb.anchors, normalization.method = "SCT")
#we can integrate data with genes list ,so that it appears in the plots,for that use below line to create the seurat object 
#combined <- IntegrateData(anchorset = ifnb.anchors, normalization.method = "SCT",features.to.integrate = rownames(ifnb.list[[1]]))

#combined
saveRDS(combined, file = "brain.rds")
#loading the integrated seurat file to do analysis as a single cell pipeline
combined.sct <- readRDS("/home/vc514/brain_combined.rds")
brain.combined.sct <- RunPCA(combined.sct, verbose = FALSE)
brain.combined.sct <- RunUMAP(brain.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
brain.combined.sct <- FindNeighbors(brain.combined.sct, reduction = "pca", dims = 1:30)
brain.combined.sct <- FindClusters(brain.combined.sct, resolution = 0.3)
plot1 <- DimPlot(brain.combined.sct, reduction = "umap")

# Violin plots for cluster-specific genes
#plot2 <- VlnPlot(brain.combined.sct, features = c("CSNK2A1", "CSNK2A2", "CSNK2B"), group.by = )
#TopMarkers(brain.combined.sct, group.by = "seurat_clusters", n = 10) # Identify top markers
#plot3 <- DoHeatmap(brain.combined.sct, features = brain.combined.sct@assays$RNA@var.top10, group.colors = c("blue", "red"))
#png("plot2.png")

brain.combined.sct@reductions
table(brain.combined.sct@active.ident)
brain_obj_markers <- FindAllMarkers(brain.combined.sct, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)

brain.combined.sct
Barcodes <- as.data.frame(colnames(brain.combined.sct))
colnames(Barcodes) <- c('Barcodes')
colnames(Barcodes)
for (i in 1:length(rownames(Barcodes))){
    Barcodes$slide_id[i] <- unlist(strsplit(Barcodes$Barcodes[i], split = '_', fixed = TRUE))[length(unlist(strsplit(Barcodes$Barcodes[i], split = '_', fixed = TRUE)))]
}
rownames(Barcodes) <- Barcodes$Barcodes
Barcodes$Barcodes <- NULL
Barcodes
brain.combined.sct <- AddMetaData(brain.combined.sct,metadata = Barcodes) 
table(brain.combined.sct$slide_id,brain.combined.sct$seurat_clusters) 
#To see how the object is divided into clusters and which images are used  
#colnames(ifnb.list[[1]])
comb_meta <- brain.combined.sct@meta.data
slide1_meta <- subset(comb_meta,comb_meta$slide_id == 1)
slide1_meta  
for (i in 1:length(rownames(slide1_meta))){
  slide1_meta$newrownames[i] <- unlist(strsplit(rownames(slide1_meta)[i], split = '_', fixed = TRUE))[1]
}
rownames(slide1_meta) <- slide1_meta$newrownames
slide1_meta$newrownames <- NULL
ifnb.list[[1]] <- AddMetaData(ifnb.list[[1]],metadata = slide1_meta)
ifnb.list[[1]]


brain_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top20 <- brain_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

fs <- c("CSNK2A1","CSNK2A2","CSNK2B")  

plot4 <- DotPlot(brain.combined.sct, features = fs) + theme_minimal()

#plot3 <- DoHeatmap(brain.combined.sct, features = top20$gene) + NoLegend()
png("plot4.png")

#p2 <- SpatialFeaturePlot(brain.combined.sct, fs, alpha = c(0.8, 1))
#png("p2.png")



