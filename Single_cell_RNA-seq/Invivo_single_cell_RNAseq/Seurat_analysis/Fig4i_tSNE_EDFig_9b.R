# Objective: Perform tSNE analysis and identify markers of specific subgroups using Seurat
# Use normalized data from Seurat_step1 

# Generate a Seurat object that contains tSNE and marker expression data
# Generate tSNE plots for age and identified clusters
# Generate heatmap with toop10 markers for identified clusters
####################################################################################
# Required packages 
library(Seurat)         # Seurat_2.3.4
library(ggplot2)        # ggplot2_3.0.0
library(dplyr)          # dplyr_0.7.6
sessionInfo()
####################################################################################
####################################################################################
setwd("~/Dropbox/Code_checking_SM/For_Matt/Seurat/")
setwd("~/Desktop/Dropbox/For_Matt/")

# Load data from step1
fibroblast.combined2 <- readRDS(file = "MATT_fibroblast.combined.rds")

# Perform an integrated analysis
fibroblast.combined2 <- RunTSNE(fibroblast.combined2, 
                               reduction.use = "cca.aligned", 
                               dims.use = 1:30, 
                               do.fast = T)

# Identify clusters - change resolution (mainly) and dims.used to alter the number of clusters identified
fibroblast.combined2 <- FindClusters(fibroblast.combined2, 
                                    reduction.type = "cca.aligned", 
                                    resolution = 0.30, 
                                    dims.use = 1:30)


# tSNE plot looking at distribution of young and old cells
p1 <- TSNEPlot(fibroblast.combined2, 
               do.return = T, 
               pt.size = 0.3,
               group.by = "old",
               colors.use = c("darkblue","orange"))

# Save tSNE plots with age 
pdf("tSNE_plot_age.pdf", width =4.3, height =5 )
plot_grid(p1)
dev.off()

# tSNE plot looking at the identified clusters 
p2 <- TSNEPlot(fibroblast.combined2, 
               do.label = T, 
               do.return = T, 
               pt.size = 0.3)


# Save tSNE plots with identifiedclusters 
pdf("tSNE_plot_celltypes.pdf", width =4.3, height =5 )
plot_grid(p2)
dev.off()

# save the clusters and age info to be used with PAGODA and 
# for calculating proportion of young and old in identified clusters
Seurat_aging <- p1$data$ident
save(Seurat_aging, file="Seurat_aging.RData")

Seurat_cluster2_dims30 <- p2$data$ident
save(Seurat_cluster2_dims30, file="Seurat_cluster2_dims30.RData")

####################################################################################
## Identify conserved cell type markers for each subgroup
fib.markers <- FindAllMarkers(object = fibroblast.combined2, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              thresh.use = 0.25)

# save table with markers of each cluster
write.csv(fib.markers, file="Conserved_markers_Seurat_2subgroups.csv")

# Extract the top 10 markers for each class and visualize by heatmap
top10 <- fib.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
pdf("Heatmap_Seurat_SP_markers_2_clusters_top10genes.pdf", width =6, height =4.2 )
#png("tSNE_Seurat_SP_markers_2_clusters_top10genes.png", width =400, height =400)
DoHeatmap(object = fibroblast.combined2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

# save the workspace
saveRDS(fibroblast.combined2, file = "fibroblast_combined2.rds")


