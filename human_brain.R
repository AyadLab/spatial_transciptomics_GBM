#packages needed
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(limma)
#paths
root_dir <- "~/Desktop/Rsession/GBM/Human" 

spatial_folder_name <- "./spatial.zip"

h5_mat_name <- "./filtered_feature_bc_matrix.h5"
unzip(spatial_folder_name)
unzip(spatial_folder_name, list = TRUE) 

#Created a seurat object and called this variable as brain_data
bh_data <- Seurat::Load10X_Spatial(
  data.dir = root_dir,
  filename = h5_mat_name,
  assay = "Spatial", 
  slice = "slice3",
  filter.matrix = TRUE, 
  to.upper = FALSE
)
#checking the data
class(bh_data[[]])
bh_data
#finding dimensions which shows us how many genes and spots are present in row col format ?
nrow(x = bh_data)
ncol(x = bh_data)
#no of unique genes in the sample = ? #1063#
#sum(bh_data$nFeature_Spatial ==  colSums(bh_data@assays$Spatial@counts > 0))
bh_data[["percent.mt"]] <- PercentageFeatureSet(bh_data, 
                                                   pattern = "^mt-")
#shows us plots to deal with any quality control issues
VlnPlot(
  bh_data, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot1 <- FeatureScatter(
  bh_data, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + NoLegend()
plot1

bh_subset <- subset(
  bh_data, 
  subset = nFeature_Spatial < 8000 & nFeature_Spatial > 1000 & 
    nCount_Spatial < 50000 )

print(paste("Filter out", ncol(bh_data) - ncol(bh_subset), 
            "samples because of the outlier QC metrics, with", ncol(bh_subset),
            "samples left."))
SpatialFeaturePlot(
  bh_subset, features = c("nFeature_Spatial", "nCount_Spatial")) &
  theme(legend.position = "bottom")  

#scaling,normalization,dimension reduction
brain_norm <- SCTransform(bh_subset, assay = "Spatial", verbose = FALSE)
brain_obj <- RunPCA(brain_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
brain_obj <- FindNeighbors(brain_obj, reduction = "pca", dims = 1:30)
#for community detection
brain_obj <- FindClusters(brain_obj, verbose = FALSE)
# UMAP input uses PCA 
brain_obj <- RunUMAP(brain_obj, reduction = "pca", dims = 1:30)
#plot shows how many clusters are formed
p3_dim <- DimPlot(brain_obj, reduction = "umap", label = TRUE) + NoLegend()
p4_cluster <- SpatialDimPlot(brain_obj, label = TRUE, label.size = 3) + NoLegend()
p3_dim + p4_cluster

brain_obj@reductions
table(brain_obj@active.ident)
cluster1_markers <- FindMarkers(brain_obj, ident.1 =1, min.pct = 0.25)
head(cluster1_markers, n = 20)
SpatialFeaturePlot(object = brain_obj, 
                   features = rownames(cluster1_markers)[19:21], 
                   alpha = c(0.1, 1), ncol = 3)

brain_obj_markers <- FindAllMarkers(brain_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold = 0.25)
brain_obj_markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

top20 <- brain_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)
DoHeatmap(brain_obj, features = top20$gene) + NoLegend()

fs <- c("CSNK2A1","CSNK2A2","CSNK2B")  

features <- c("CSNK2A1","CSNK2A2","COPS2","TP53","CSNK2B","PTEN","SSRP1")  # Specify the feature names
#single cell distribution in each cluster
VlnPlot(object = brain_obj,features = c("TP53")) + geom_boxplot()
DotPlot(brain_obj, features = features)
SpatialFeaturePlot(brain_moransi, features = features, ncol = 3, alpha = c(0.1, 1))
RidgePlot(brain_obj, features = fs, ncol = 2)

#gene expression found in human brain
SpatialFeaturePlot(brain_obj,fs,alpha = c(0.1, 1))

#lower expression were down weighted here
p2 <- SpatialFeaturePlot(brain_obj, fs, alpha = c(0.8, 1))
p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain_obj, idents = c(1, 2, 3, 4,
                                                                                     5)), facet.highlight = TRUE, ncol = 3)
#interactive plots
SpatialDimPlot(brain_obj, interactive = TRUE)

SpatialFeaturePlot(brain_obj, features = "CSNK2B", interactive = TRUE)

SpatialFeaturePlot(object = brain_obj, 
                   features = fs, 
                   alpha = c(0.1, 1), ncol = 3)






