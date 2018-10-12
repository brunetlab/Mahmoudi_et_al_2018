# Objective: Generate heatmap based on all DEGs identified between young and old fibroblasts in vitro
# Figure 2c
# Use DESeq2 core DEG analysis DESeq2_analysis_DEGs_VST_values.R

###########################################################################################
# Need the following package
library('pheatmap')    # v1.0.10
###########################################################################################
# Set working directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Heatmaps/Fig2c_heatmap_fibroblast_in_vitro/')

#### Extract all DEGs between young and old that goes up with age  
# Upload DEG data from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro/DESeq2_DEG_all_genes.Rda")
# Make it into a data frame
genes<- as.data.frame(resOrdered)
# Make gene names in CAPTIAL letter
row.names(genes) <- toupper(row.names(genes))
# Extract all DEGs that are FDR <0.05 and >1.5 in Fold Change, higher in young
deg <- genes[which(genes$padj < 0.05 & abs(genes$log2FoldChange) > 0.58),] # 1022 DEGs

#### Extract expression values (vst) for all DEGs between young and old that goes up with age  
# Upload vst values from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany//Core_Analysis/Fibroblast_in_vitro/DESeq2_vst_corrected_values.Rda")

# Extract vst values for DEGs
vst.deg <-as.data.frame(v.data[which(toupper(row.names(v.data)) %in% row.names(deg)),])

# Generate heatmap 
pdf("Heatmap_DESeq2_BATCH_Corrected_vst_DEGs_fibs_in_vitro_080618.pdf", width =4.3, height =5, onefile=F  )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
pheatmap(vst.deg,
         border_color = NA,
         cluster_col=F, 
         cluster_row=T, 
         scale="row", 
         annotation_legend=TRUE,
         clustering_distance_cols="manhattan", 
         clustering_method="complete",
         fontsize_col=8,
         fontsize_row=3.5,
         show_rownames=F,
         show_colnames=T,
         fontface="bold",
         treeheight_col=10,
         treeheight_row=0,
         cellwidth = NA, 
         cellheight = NA,
         color = colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100))
dev.off()
