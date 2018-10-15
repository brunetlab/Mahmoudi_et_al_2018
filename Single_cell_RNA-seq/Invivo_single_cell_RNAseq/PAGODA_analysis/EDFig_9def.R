# Objective: Make heatmap for Seurat marker genes, Fibroblast actiavtion and TNF siganling siganture 
# Use data clustering based on PAGODA
# Use normalized gene expression from PAGODa (var.info)
# Use geneset used (KEGG_2016_disease_aging_activation_lower_case_for_PAGODA_090518.RData)

# Marker genes comes from: Seurat_Single_Cell_Analysis.R

################################################################################################################
# packages required
library(Biobase)
library(scde)
library(RColorBrewer)
library(Rtsne)
################################################################################################################
# Load files required 
setwd("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/Heatmap_expression_transcriptional_signatures/Expected_results/")

# load required data 
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/code_checking_repeat/PAGODA_set_seed/Varinfo/varinfo.RData")
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/code_checking/R_analysis_gene5500_mito10perc/Pathways/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA_090518.RData")
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/Expected_Data/cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.repeat_set_seed.RData")
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/code_checking/R_analysis_gene5500_mito10perc/PAGODA/Expected_results/Varinfo/varinfo.RData")

# load kegg pathways 
gene.sets <- read.table("~/Dropbox/Code_checking_SM/For_Matt/Essentials/Pathways/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA_090518.txt", sep="\t", header = F, quote="",row.names = 1)
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
mat.FA[mat.FA < -6] <- -6
mat.FA[mat.FA > 6] <- 6

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
#png("Heatmap_Fibroblast_activation_genes.png", width=4, height=5, units="in", res=300)
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
#png("Heatmap_Fibroblast_activation_genes.png", width=4, height=5, units="in", res=300)
heatmap(mat.ck[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()


##############################################################################################
# Seurat marker genes
##############################################################################################
# Based on Seurat analysis
# Seurat marker genes
markers <- c("GAPDH", "MIF", "PKM", "LDHA", "ENO1", 
             "IL1RL1","TIMP1","ACTA2","FABP5","SPP1",
             "CLEC3B","SEPP1","SERPING1","APOD","IGFBP4",
             "GAS6","IGF1","IFI27L2A","PI16","SFRP2")

# Extract marker genes
mat.marker <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(markers)),]

# determine range of expression
range(mat.marker)
#[1]  -8.546065  9.517986

# Balance expression for the sake of color gradient
mat.marker[mat.marker < -8.546065] <- -8.546065
mat.marker[mat.marker > 8.546065] <- 8.546065

pdf("Heatmap_Seurat_markers.pdf", width =4.3, height =5, onefile=F  )
heatmap(mat.marker[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()


##############################################################################################
# cytokine-cytokine interaction and TNF signaling pathway genes
##############################################################################################
ck_tnf <- pagoda.show.pathways(c("Cytokine-cytokine receptor interaction","TNF signaling pathway"), 
                           varinfo_gene5500_mito10perc, 
                           kegg_aging.env, 
                           cell.clustering = hc,
                           n.genes =8,
                           margins = c(1,5), 
                           show.cell.dendrogram = TRUE, 
                           showRowLabels = TRUE, 
                           showPC = TRUE,
                           return.details = TRUE)
row.names(ck_tnf$rotation)
# extract the genes
genes.ck_tnf <- row.names(ck_tnf$rotation)

# heatmap for subset of cytokine genes
mat.ck_tnf <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(genes.ck_tnf)),]

# determine range of expression
range(mat.ck_tnf)
#[1]  -5.234634  8.902922

# Balance expression for the sake of coclor gradient
mat.ck_tnf[mat.ck_tnf < -5.234634] <- -5.234634
mat.ck_tnf[mat.ck_tnf > 5.234634] <- 5.234634

# Generate heatmap
pdf("Heatmap_Cytokine-cytokine_receptor_interaction__TNF_signaling_pathway_genes_top8.pdf", width =4.3, height =5, onefile=F  )
#png("Heatmap_Fibroblast_activation_genes.png", width=4, height=5, units="in", res=300)
heatmap(mat.ck_tnf[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()


##############################################################################################
# Complement and coagulation cascades and NF−kappa B signaling pathway genes 
##############################################################################################
CC_NFKB <- pagoda.show.pathways(c("Complement and coagulation cascades","NF-kappa B signaling pathway"), 
                               varinfo_gene5500_mito10perc, 
                               kegg_aging.env, 
                               cell.clustering = hc,
                               n.genes =8,
                               margins = c(1,5), 
                               show.cell.dendrogram = TRUE, 
                               showRowLabels = TRUE, 
                               showPC = TRUE,
                               return.details = TRUE)
row.names(CC_NFKB$rotation)
# extract the genes
genes.CC_NFKB <- row.names(CC_NFKB$rotation)

# heatmap for subset of cytokine genes
mat.CC_NFKB <- varinfo_gene5500_mito10perc$mat[which(toupper(row.names(varinfo_gene5500_mito10perc$mat)) %in% toupper(genes.CC_NFKB)),]

# determine range of expression
range(mat.CC_NFKB)
#[1]  -4.437894  5.496788

# Balance expression for the sake of coclor gradient
mat.CC_NFKB[mat.CC_NFKB < -4.437894] <- -4.437894
mat.CC_NFKB[mat.CC_NFKB > 4.437894] <- 4.437894

# Generate heatmap
pdf("Heatmap_complement_and_coaulation_NFKb_genes_top8.pdf", width =4.3, height =5, onefile=F  )
#png("Heatmap_Fibroblast_activation_genes.png", width=4, height=5, units="in", res=300)
heatmap(mat.CC_NFKB[,hc$labels], 
        Colv=as.dendrogram(hc), 
        Rowv=T, 
        scale="none",
        na.rm =T,
        labCol=F,
        col=colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100),
        ColSideColors=pagoda)
dev.off()
