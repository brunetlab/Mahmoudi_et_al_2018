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
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Heatmaps/Fig3i_heatmap_CKs_in_vivo/')

#############################################################################################################
# Upload vst values from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vivo/DESeq2_DEG_young_old_in_vivo_all_genes.Rda")

# Make a Young/Old signature expression based on all the genes that are significantly lower/hihger with age
# Upload DEG data from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vivo/DESeq2_young_old_in_vivo_vst_values.Rda")
# Make it into a data frame
genes<- as.data.frame(resOrdered)
# Extract all DEGs that gow up with afe 
deg.up <- genes[which(genes$padj < 0.05 & genes$log2FoldChange < 0),] #189
deg.down <- genes[which(genes$padj < 0.05 & genes$log2FoldChange > 0),] #144

young.up.id <- v.data[which(row.names(v.data) %in% row.names(deg.down)),]
young.up.sig <- colSums(young.up.id)
old.up.id <- v.data[which(row.names(v.data) %in% row.names(deg.up)),]
old.up.sig <- colSums(old.up.id)

v.data <- rbind(v.data,
                young.up.sig,
                old.up.sig)
tail(v.data)

#############################################################################################################
# Minimum subset to changes with age across both populations

# Alternative name
# Ccl8 = Mcp2
# Ccl11 = Eotaxin
# Ccl2 = Mcp1
# Ccl7 = Mcp3
# Figf = VEGFd
# Cxcl12 = Sdf1

subset.small <- c("Wnt2", "Tgfb3", "Lif", 
                  "Tgfb2","Il33","Il6", "Fgf7","Figf","Ccl2","Ccl7","Il1b","Ccl8", "Cxcl12", "Il15", "Tnf",
                  "Ccl11", "Il7", "Vegfc","young.up.sig","old.up.sig")



# Extract vst for ck genes 
vst.subset.small <- v.data[which(toupper(row.names(v.data)) %in% toupper(subset.small)),]

pdf("Heatmap_CK_in_vivo_fibroblasts_180805.pdf", width =4.3, height =5, onefile=F  )
pheatmap(vst.subset.small, 
         cluster_col=T, 
         cluster_row=T, 
         scale="row", 
         annotation_legend=TRUE,
         clustering_distance_cols="correlation", 
         clustering_method="complete",
         fontsize_col=6,
         fontsize_row=8,
         show_rownames=T,
         show_colnames=T,
         fontface="bold",
         treeheight_col=20,
         treeheight_row=0,
         cellwidth = 10, 
         cellheight = NA,
         border_color = NA,
         color = colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100))
dev.off()
