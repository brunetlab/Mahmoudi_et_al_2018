# Objective: Generate heatmap plots for pagoda analysis of young and old fibroblast in vivo using KEGG Activation geneset
# Fig. 4j and Extended Data Fig. 9c  
# Use data generated using the PAGODA_step1_core_analysis.R and PAGODA_step2_core_analysis.R scripts:
## young_old_norm_counts_gene5500_mito10perc.Rdata (normalized gene counts)
## knn.error.models.RData
## varinfo.RData 
## geneset used (KEGG_2016_disease_aging_activation_lower_case_for_PAGODA.RData)
## pwpca_top.aspects_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_varnorm_knn.error.models.RData
## clpca_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_novel_gene_clusters.RData

################################################################################################################
# packages required
library(Biobase)       # Biobase_2.40.0
library(scde)          # scde_2.8.0
library(RColorBrewer)  # RColorBrewer_1.1-2
sessionInfo()
################################################################################################################
setwd("/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/")

# Load the files required 
load("Single_cell_RNA-seq/Invivo_single_cell_RNAseq/Essentials/Expression_data_PAGODA/young_old_raw_counts_gene5500_mito10perc.Rdata")
load("/PAGODA/knn.error.models/knn.error.models.RData") 
load("/PAGODA/Varinfo/varinfo.RData") 
load("Single_cell_RNA-seq/Invivo_single_cell_RNAseq/Essentials/Pathways/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA.RData")
load("/PAGODA/KEGG_disease_aging/pwpca_top.aspects_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_varnorm_knn.error.models.RData")
load("/PAGODA/KEGG_disease_aging/clpca_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_novel_gene_clusters.RData")

# ## Determine overall cell clustering (hclust) based on weighted correlation of genes underlying the top aspects transcriptional heterogeneity.
# ## For some reason, return.table and return.genes have to be FALSE in order to get the list with all results for clustering:
tam <- pagoda.top.aspects(pwpca,
                          return.table = FALSE,
                          return.genes = FALSE,
                          plot = FALSE,
                          z.score = 4)
save(tam, file = "tam_top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.RData")

hc <- pagoda.cluster.cells(tam,
                           varinfo = varinfo_gene5500_mito10perc,
                           verbose = 1,
                           include.aspects = FALSE,
                           return.details = TRUE) ## Whether to return also a the distance matrix and gene values
## clustering cells based on 174 genes
#### 195 genes!
save(hc, file = "cluster.cells_return.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.RData")

## I need to run it again with return.details = FALSE in order to pass it to the plot functions:
rm(hc)
hc <- pagoda.cluster.cells(tam,
                           include.aspects = FALSE,
                           varinfo = varinfo_gene5500_mito10perc)
save(hc, file = "cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.RData")

# Input color inforamtion to be used by PAGODA
# load Seurat clusters : Seurat_cluster2_dim30 and input color of choice
load("Single_cell_RNA-seq/Invivo_single_cell_RNAseq/Essentials/Expression_data_PAGODA/Seurat_cluster2_dims30.RData")

Seurat <- as.character(Seurat_cluster2_dims30)
Seurat[Seurat==0] <- "#F8766D"
Seurat[Seurat==1] <- "#00BFC4"

# Input aging clusters
age <- c(rep("orange",length(grep("Young", colnames(young_old_raw_counts_gene5500_mito10perc)))),
         rep("darkblue",length(grep("Old", colnames(young_old_raw_counts_gene5500_mito10perc)))))

# decide on number of PAGODA clusters and colors
n.clusters <- 6
col.pal <- brewer.pal(n.clusters, "Paired")
pagoda <- col.pal[cutree(hc, n.clusters)]
# merge all info into 1 data frame
col.cols <- rbind(pagoda, Seurat, age)

## Visualize significant aspects of heterogeneity
# Show top aspects 
pdf("PAGODAplot_top.aspects_KEGG_DISEASE_AGING_ALL.pdf", width = 10, height = 10)
png("PAGODAplot_top.aspects_KEGG_DISEASE_AGING_ALL.png", width = 600, height = 600)
tamr <- pagoda.reduce.loading.redundancy(tam,
                                         pwpca,
                                         cell.clustering = hc,
                                         distance.threshold = 0.01, ## similarity threshold for grouping interdependent aspects, default: 0.01
                                         abs = TRUE, ## Whether to use absolute correlation.
                                         plot = TRUE,
                                         col.cols = col.cols,
                                         box = TRUE,
                                         labCol = NA,
                                         top=10,
                                         cols = colorRampPalette(c("cornflowerblue", "white", "palevioletred4"), ## low neutral high
                                                                 space = "Lab")(1024),
                                         margins = c(0.5, 15))
dev.off()

save(tamr, file = "tamr_top.aspects_KEGG_DISEASE_AGING_ALL.RData")
