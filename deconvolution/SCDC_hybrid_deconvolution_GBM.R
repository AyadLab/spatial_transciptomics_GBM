#SCDC deconvolutiona and exploring the seurat v5 package
#Deconvolution helps us to visualize how many cell types are present in each spot 
library(Seurat)
library(dplyr)
library(cowplot)
library(viridis)
#install scdc package
# if (!("xbioc" %in% rownames(inst))) {
#   remotes::install_github("renozao/xbioc", dependencies = FALSE)
# }
# if (!("SCDC" %in% rownames(inst))) {
#   remotes::install_github("meichendong/SCDC", dependencies = FALSE)
# }

suppressPackageStartupMessages(require(SCDC))
suppressPackageStartupMessages(require(Biobase))
#bainn is our spatial object and reference object gives us cell type information as it is supervised algorithm
brain <- readRDS("vc514_sct.rds")
selected_slide_ids <- c(6, 8,9, 11)
selected_slide_ids <- c(7)
subset_brain_7 <- brain[, brain$slide_id %in% selected_slide_ids]
reference <- readRDS("isosceles_suterJohnsonMerge_fourStates.RDS")

# After subsetting, we renormalize 
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
reference <- SCTransform(reference, ncells = 3000, verbose = FALSE, method = "poisson") %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

#3
# the annotation is stored in the 'subclass' column of object metadata
s1 <- DimPlot(reference, group.by = "CellType", label = TRUE)
pdf("dimplot_ref.pdf")
s1
dev.off()
anchors <- FindTransferAnchors(reference = reference, query = subset_brain_8, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = reference@meta.data$CellType, prediction.assay = TRUE,
                                  weight.reduction = subset_brain_8[["pca"]], dims = 1:30)
subset_brain_8[["predictions"]] <- predictions.assay
DefaultAssay(subset_brain_8) <- "Spatial"
#Transfer predictions from reference to spatial object brain and visualize the plotting through spatialfeatureplot
s2_11 <- SpatialFeaturePlot(subset_brain_8, features = c("T-Cells","Myeloid"), pt.size.factor = 1.5, ncol = 2 , crop = TRUE)

pdf("spatial_feature_plot_primary.pdf",width = 25,height =16)
s2_11
dev.off()
#See all the features spatial imaging at once
features <- unique(reference@meta.data$CellType)
spatial_plot_ls <- list()
for (feature in features) {
  plot <- SpatialFeaturePlot(
    subset_brain_7, 
    features = feature, 
    pt.size.factor = 2, 
    ncol = 3, 
    crop = TRUE
  )
  spatial_plot_ls[[feature]] <- plot
}
combined_plot <- plot_grid(plotlist = spatial_plot_ls, nrow =3)
pdf("spatial_feature_plot_slice_human_n.pdf",width =30,height =10)
combined_plot
dev.off()

reference@active.assay = "RNA"

markers_sc <- FindAllMarkers(reference, only.pos = TRUE, logfc.threshold = 0.1,
                             test.use = "wilcox", min.pct = 0.05, min.diff.pct = 0.1, max.cells.per.ident = 200,
                             return.thresh = 0.05, assay = "RNA")

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(subset_brain_7), ]


# Select top 20 genes per cluster, select top by first p-value, then absolute
# diff in pct, then quota of pct.
markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                               1))
markers_sc %>%
  group_by(cluster) %>%
  top_n(-100, p_val) %>%
  top_n(50, pct.diff) %>%
  top_n(20, log.pct.diff) -> top20
m_feats <- unique(as.character(top20$gene))

eset_SC <- ExpressionSet(assayData = as.matrix(allen_reference@assays$RNA@counts[m_feats,
]), phenoData = AnnotatedDataFrame(allen_reference@meta.data))
eset_ST <- ExpressionSet(assayData = as.matrix(subset_brain_7@assays$Spatial@counts[m_feats,
]), phenoData = AnnotatedDataFrame(subset_brain_7@meta.data))

deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "CellType",
                                     ct.sub = as.character(unique(eset_SC$CellType)))
#prop gives us proportions of cell types in each spot i.e RNA sequence
prop <- deconvolution_crc$prop.est.mvw

#visualization done through SPOTLIGHT method 
neighbourhood_interaction <- plotInteractions(prop, which = "heatmap", metric = "prop")

ct <- colnames(prop)
prop[prop < 0.1] <- 0

# Define color palette
# (here we use 'paletteMartin' from the 'colorBlindness' package)
paletteMartin <- c(
  "#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", 
  "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", 
  "#920000", "#924900", "#db6d00", "#24ff24")

pal <- colorRampPalette(paletteMartin)(length(ct))
names(pal) <- ct

scatterpie_image <- plotSpatialScatterpie(
  x = subset_brain_7,
  y = prop,
  cell_types = colnames(prop),
  img = FALSE,
  scatterpie_alpha = 1,
  pie_scale = 0.3) +
  scale_fill_manual(
    values = pal,
    breaks = names(pal))

pdf("gbm_recurrent_neftelstate_.pdf")
scatterpie_image
dev.off()

pdf("plotinteractions_gbm_recurrent_GBM.pdf")
neighbourhood_interaction
dev.off()
