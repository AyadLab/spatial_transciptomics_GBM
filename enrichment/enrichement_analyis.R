#Enrichment analysis on our integrated seurat object called data
#packages to be loaded
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(png)
library(DESeq2)
library(RColorBrewer)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

#to find heatmap between cluster and gene and understand the relationship between them
data <- readRDS("~/Desktop/vc514_sct.rds")
markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers
top_markers <- markers %>% group_by(cluster) %>% top_n(5, wt = avg_log2FC)
top_markers$gene
heatmap_gene <- DoHeatmap(data,
                          features = top_markers$gene,
                          group.by = "seurat_clusters")
pdf("~/Desktop/heatmap_(cluster-gene).pdf",height =18 ,width = 25)
print(heatmap_gene)
dev.off()

##to find heatmap between cluster and BP 
cluster_results <- list()

for (i in 1:10) {
  all_de_markers <- FindMarkers(data,ident.1 = i,only.pos = TRUE)
  cluster_results[[i]] <- all_de_markers
}

# Get all marker genes across all clusters
all_markers <- do.call(rbind, cluster_results)

# Extract gene names from marker genes
gene_names <- rownames(all_markers)
gene_names
# Perform GO enrichment analysis
enrichment_data_list <- list()

for (i in 1:10) {
  cluster_id <- paste0(i)
  genes <- rownames(cluster_results[[i]])
  enrich_result <- enrichGO(gene = genes, keyType = "SYMBOL", ont = "BP", OrgDb = org.Hs.eg.db)
  cluster_num <- as.integer(strsplit(cluster_id, "")[[1]][1])
  df <- as.data.frame(enrich_result@result)
  df$Cluster <- cluster_num
  enrichment_data_list[[cluster_id]] <- df
  cat("cluster number ::", cluster_id, "\n")
  print(enrich_result)
  cat("\n")
}

enrichment_data <- do.call(rbind, enrichment_data_list)
enrichment_data <- enrichment_data[order(enrichment_data$Cluster, enrichment_data$pvalue), ]
top_5_enriched <- by(enrichment_data, enrichment_data$Cluster, function(x) head(x, 5))
top_5_enriched <- do.call(rbind, top_5_enriched)
top_5_enriched

cluster_data <- top_5_enriched %>%
  select(Cluster = Cluster,ID,Biological_Process = Description, p.adjust) %>%
  group_by(Cluster) %>%
  top_n(-5, p.adjust) 
cluster_data

heatmap_data <- as.matrix(t(table(cluster_data$Biological_Process,cluster_data$Cluster)))
heatmap_data

pdf("~/Desktop/heatmapBP1.pdf")

ps <- heatmap(
  heatmap_data,
  Rowv = TRUE,   
  Colv = TRUE, 
  scale = "none",  
  col = colorRampPalette(c("blue", "white", "red"))(100), 
  main = "Cluster-Biological Process Heatmap",
  xlab = NULL,   
  ylab = "Cluster",
  cexRow = 0.5,  
  cexCol = 0.5,  
  offset = c(1, 0.5), 
  margins = c(8, 5),  
  las = 2,  
  key = TRUE  
)
print(ps)
dev.off()
