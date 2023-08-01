#Stddeconvolve analysis
#Loading all the results for visualization of single image
#Packages required
library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(STdeconvolve)
library(topicmodels)
library(ldatuning)
library(LDAvis)
library(stringr)
library(cowplot)

#loading rds file which is the normalized integrated seurat object
brain <- readRDS("~/Desktop/vc514_sct.rds")
image_index <- 0

# Extract counts for the specific image
counts1 <- brain@assays$RNA@counts[, image_index]
counts1
cd <- brain@assays$RNA@counts

counts <- cleanCounts(counts = cd,
                        min.lib.size = 100,
                        min.reads = 1,
                        min.detected = 1,
                        verbose = TRUE)

dim(counts)
dim(brain)
coords <- brain@images$slice_human@coordinates
pos <- coords[c('imagerow', 'imagecol')]
colnames(pos) <- c('x','y')
annot  <- as.data.frame(brain$seurat_clusters)
annot <- subset(annot,rownames(annot) %in% rownames(pos))
annot2 <- as.factor(annot$`brain$seurat_clusters`)
names(annot2) <- rownames(annot) 
dim(annot2)
dim(pos)

#saveRDS(combined, file = "brain.rds")
#brain_seurat <- readRDS("brain.rds")
deconProp <- readRDS("~/Desktop/vc514_deconP.rds")
deconGexp <- readRDS("~/Desktop/vc514_deconG.rds")
#Visualization

plot <- vizAllTopics(deconProp, pos,r=0.4)
print(plot)
plot
png("celltype.png")
dev.off()

#intersection <- intersect(rownames(pos), rownames(deconProp))
#intersection
plot2 <- vizAllTopics(deconProp, pos, 
             groups = annot2, 
             group_cols = rainbow(length(levels(annot2))),
             r=0.4)
print(plot2)
plot2
pdf("~/Desktop/celltype1.pdf",height = 3.5,width = 3.5)
plot2
dev.off()

dim(deconProp)
dim(annot2)
dim(pos)

pos_list <- list()

# Define the image names
image_names <- c('slice_human', 'slice_human.1', 'slice_human.2', 'slice_human.3', 'slice_human.4',
                 'slice_human.5', 'slice_human.6', 'slice_human.7', 'slice_human.8', 'slice_human.9',
                 'slice_human.10')

# Loop over the image names and extract the positions
for (image_name in image_names) {
  coords <- brain@images[[image_name]]@coordinates
  pos1 <- coords[c('imagerow', 'imagecol')]
  colnames(pos1) <- c('x', 'y')
  pos_list[[image_name]] <- pos
}
pos_list[['slice_human.2']]
pos_list
viridis(length(levels(annot2)))
----
max(pos)/nrow(pos) * 4
plt1 <- vizAllTopics(theta = deconProp,
                      pos = pos,
                      groups = annot2, 
                      group_cols = viridis(length(levels(annot2rm)), option = 'turbo'),
                      r = max(0.4, max(pos)/nrow(pos) * 15),
                      lwd = 0.5,
                      showLegend = TRUE,
                      plotTitle = NA) +
  ggplot2::guides(fill=ggplot2::guide_legend(ncol=2)) +
  
  ggplot2::geom_rect(data = data.frame(pos),
                     ggplot2::aes(xmin = min(x)-90, xmax = max(x)+90,
                                  ymin = min(y)-90, ymax = max(y)+90),
                     fill = NA, color = "black", linetype = "solid", size = 0.5) +
  
  ggplot2::theme(
    plot.background = ggplot2::element_blank()
  ) +
    ggplot2::guides(colour = "none")
pdf("~/Desktop/cell_10.pdf",height =14 ,width = 14)  
plt
dev.off()
#-----#
# Combine x and y values from pos_list into a single dataframe
combined_pos <- do.call(rbind, pos_list)
# Extract image names from combined_pos
rownames(combined_pos) <- sub(".*\\.", "", rownames(combined_pos))
combined_pos
#---------#
##trying to creta a function that takes all the positions
DecovPlot <- function(pos_image_n, annot2) {
  plt <- vizAllTopics(theta = deconProp,
                      pos = pos_image_n,
                      groups = annot2,
                      group_cols = viridis::viridis(length(levels(annot2)), option = 'turbo'),
                      r = max(0.4, max(pos_image_n$y) / nrow(pos_image_n) * 15),
                      lwd = 0.5,
                      showLegend = TRUE,
                      plotTitle = NA) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::geom_rect(data = data.frame(pos_image_n),
                       ggplot2::aes(xmin = min(x) - 90, xmax = max(x) + 90,
                                    ymin = min(y) - 90, ymax = max(y) + 90),
                       fill = NA, color = "black", linetype = "solid", size = 0.5) +
    ggplot2::theme(plot.background = ggplot2::element_blank()) +
    ggplot2::guides(colour = "none")
  
  return(plt)
}
# plots_decov : list to store plots(creating list of lists)#
plots_decov <- list()
# Looping the image names and create plots for each image#
for (image_name in image_names) {
  pos_image_n <- pos_list[[image_name]]
  plot <- DecovPlot(pos_image_n, annot2)
  plots_decov[[image_name]] <- plot
}
#cow plot




#####

pcs.info <- stats::prcomp(t(log10(as.matrix(counts)+1)), center=TRUE)
nPcs <- 7 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]  
-----------------------

DecovPlot(pos_list[["slice_human.4"]],annot2)
image_nos <- 1 : length(pos_list)
labels <- Paste("image",image_nos)
plotgrid(plotlist = plots_decov,labels = labels,ncol = 4)
#-----------
dim(annot2)
dim(pos_list[["slice_human"]])
pos_list[["slice_human"]]
