# Objective: Perform Pathway analysis on DEGs (FDR < 0.05, and >1.5 Fold Change) using Fischer's exact test
# Generate a heatmaps for the different classes of pathways
# Figure XX
# Use DEG table (with fold change and FDR-values) created by DEseq2  

# Will use KEGG pathways downloaded from GSEA Broad Institute

###########################################################################################
# Need the following package
library('pheatmap')      # v1.0.7
library('DESeq2')        # v1.20.0
###########################################################################################
# Set working directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Heatmaps/Fig.2h_heatmap_human/')

#############################################################################################################
# Load vst values from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Human_healthy/DESeq2_vst_human_Fibs_healthy_all_genes_linear_model.Rda")

# Make a Young/Old signature expression based on all the genes that are significantly lower/hihger with age
# Load DEG data from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Human_healthy/DESeq2_DEG_human_Fibs_healthy_all_genes_linear_model.Rda")
# Make it into a data frame
genes<- as.data.frame(resOrdered)

# Extract all DEGs that change with age 
deg <- genes[which(genes$padj < 0.05 ),] #323
vstMat.deg <- vstMat[which(row.names(vstMat) %in% row.names(deg)),]

# reorder samples based on age 
colnames(vstMat.deg)
vstMat.deg.order <- vstMat.deg[,c(grep("_0Y", colnames(vstMat.deg)),
                                  grep("_1Y", colnames(vstMat.deg)),
                                  grep("_29Y", colnames(vstMat.deg)),
                                  grep("_3OY", colnames(vstMat.deg)),
                                  grep("_43Y", colnames(vstMat.deg)),
                                  grep("_66Y", colnames(vstMat.deg)),
                                  grep("_69Y", colnames(vstMat.deg)),
                                  grep("_71Y", colnames(vstMat.deg)),
                                  grep("_87Y", colnames(vstMat.deg)),
                                  grep("_89Y", colnames(vstMat.deg)))]

# generate a heatmap
pdf("Heatmap_human fibroblasts_1801004.pdf", width =4.3, height =5, onefile=F  )
pheatmap(vstMat.deg.order, 
         cluster_col=F, 
         cluster_row=T, 
         scale="row", 
         annotation_legend=TRUE,
         clustering_distance_cols="correlation", 
         clustering_method="complete",
         fontsize_col=6,
         fontsize_row=8,
         show_rownames=F,
         show_colnames=T,
         fontface="bold",
         treeheight_col=20,
         treeheight_row=0,
         cellwidth = 10, 
         cellheight = NA,
         border_color = NA,
         color = colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100))
dev.off()
