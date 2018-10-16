# Objective: Make heatmap for Seurat marker genes, Fibroblast actiavtion and TNF siganling siganture 
# Extended Data Fig. 9d-f
# Use data generated using the PAGODA_step1_core_analysis.R, PAGODA_step2_core_analysis.R and Fig_4j_EDFig_9c.R script scripts:
## 
## varinfo.RData 
## geneset used (KEGG_2016_disease_aging_activation_lower_case_for_PAGODA.RData)

# Marker genes comes from: Seurat_Single_Cell_Analysis.R

################################################################################################################
# packages required
library(Biobase)
library(scde)
library(RColorBrewer)
library(Rtsne)
################################################################################################################
# Load files required 
setwd("/PAGODA/Heatmap_expression_transcriptional_signatures/")

# Load required data
load("/PAGODA/knn.error.models/knn.error.models.RData") 
load("/PAGODA/Varinfo/varinfo.RData") 
load("/PAGODA/KEGG_disease_aging/pwpca_top.aspects_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_varnorm_knn.error.models.RData")
load("/PAGODA/KEGG_disease_aging/clpca_KEGG_DISEASE_AGING_ACTIVATION_SENESCENCE_novel_gene_clusters.RData")
load("/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/Expected_Data/cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.repeat_set_seed.RData")

# load kegg pathways 
load("/Single_cell_RNA-seq/Invivo_single_cell_RNAseq/Essentials/Pathways/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA.RData")
gene.set.data <- as.data.frame(t(gene.sets))

# Input color inforamtion to be used by PAGODA
# decide on number of PAGODA clusters and colors
n.clusters <- 6
col.pal <- brewer.pal(n.clusters, "Paired")
pagoda <- col.pal[cutree(hc, n.clusters)]

# Visualize gene expression within PAGODA plot
# Fibroblast activation genes
##############################################################################################
genes.FA <- gene.set.data$`Fibroblast activation`
genes.FA <- na.omit(genes.FA)

# Extract expression data of genes of interest
mat.FA <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(genes.FA)),]
row.names(mat.FA)
# determine range of expression 
range(mat.FA)
#[1]  -21.920585   5.956637

# Balance expression for the sake of color gradient
mat.FA[mat.FA < -5.956637] <- -5.956637
mat.FA[mat.FA > 5.956637] <- 5.956637

# Generate heatmap 
pdf("Heatmap_Fibroblast_activation_genes.pdf", width =4.3, height =5, onefile=F  )
heatmap(mat.FA[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()

##############################################################################################
# TNF signalling genes
##############################################################################################
# identify the top 30 genes driving TNF signlaing
tnf <- pagoda.show.pathways("TNF signaling pathway", 
                     varinfo_gene5500_mito10perc, 
                     kegg_aging.env, 
                     cell.clustering = hc,
                     n.genes =30,
                     margins = c(1,5), 
                     show.cell.dendrogram = TRUE, 
                     showRowLabels = TRUE, 
                     showPC = TRUE,
                     return.details = TRUE)
row.names(tnf$rotation)
# extract the genes
genes.tnf <- row.names(tnf$rotation)

# heatmap for subset of cytokine genes
mat.tnf <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(genes.tnf)),]

# determine range of expression
range(mat.tnf)
#[1]  -5.234634  8.902922

# Balance expression for the sake of coclor gradient
mat.tnf[mat.tnf < -5.234634] <- -5.234634
mat.tnf[mat.tnf > 5.234634] <- 5.234634

# Generate heatmap
pdf("Heatmap_TNF_signaling_genes_top30.pdf", width =4.3, height =5, onefile=F  )
heatmap(mat.tnf[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()

##############################################################################################
# cytokine-cytokine interaction genes
##############################################################################################
ck <- pagoda.show.pathways("Cytokine-cytokine receptor interaction", 
                           varinfo_gene5500_mito10perc, 
                           kegg_aging.env, 
                           cell.clustering = hc,
                           n.genes =30,
                           margins = c(1,5), 
                           show.cell.dendrogram = TRUE, 
                           showRowLabels = TRUE, 
                           showPC = TRUE,
                           return.details = TRUE)
row.names(ck$rotation)
# extract the genes
genes.ck <- row.names(ck$rotation)

# heatmap for subset of cytokine genes
mat.ck <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(genes.ck)),]

# determine range of expression
range(mat.ck)
#[1]  -5.234634  9.280114

# Balance expression for the sake of coclor gradient
mat.ck[mat.ck < -5.234634] <- -5.234634
mat.ck[mat.ck > 5.234634] <- 5.234634

# Generate heatmap
pdf("Heatmap_Cytokine-cytokine_receptor_interaction_genes_top30.pdf", width =4.3, height =5, onefile=F  )
heatmap(mat.ck[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()


