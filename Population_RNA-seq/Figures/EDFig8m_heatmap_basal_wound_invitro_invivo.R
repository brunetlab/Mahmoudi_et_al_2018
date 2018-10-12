# Objective: Generate heatmap based on all pathways that are significatn between comparison of young and old
# in vivo and in vitro, as well as young an dold basal vs wound
# Extended Data Fig. 8m
# Use DESeq2 core DEG analysis DESeq2_analysis_DEGs_VST_values.R

###########################################################################################
# Need the following package
library('pheatmap')    # v1.0.10
###########################################################################################
# Set working directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Heatmaps/ED_Fig8m_heatmap_pathways/')

#### Extract all DEGs between young and old that goes up with age  
# load pathway enrichment analyses results
load("~/Dropbox/Code_checking_SM/For_Brittany/Tables/Katja_Extended_Data_Table_invivo_old_basal_vs_wounds/kegg.results.pathways.kegg.results.old_basal_wounds.RData")
old_basal_wounds <- as.data.frame(kegg.results.old_basal_wounds$pathway.results)
old_basal_wounds.o <- old_basal_wounds[order(rownames(old_basal_wounds)),]
old_basal_wounds.o$name <- rownames(old_basal_wounds.o)

load("~/Dropbox/Code_checking_SM/For_Brittany/Tables/Katja_Extended_Data_Table_invivo_young_basal_vs_wounds/kegg.results.pathways.kegg.results.young_basal_wounds.RData")
young_basal_wounds <- as.data.frame(kegg.results.young_basal_wounds$pathway.results)
young_basal_wounds.o <- young_basal_wounds[order(rownames(young_basal_wounds)),]
young_basal_wounds.o$name <- rownames(young_basal_wounds.o)

load("~/Dropbox/Code_checking_SM/For_Brittany/Tables/Katja_Extended_Data_Table_young_old_invitro/kegg.results.pathways.kegg.results.young_old_invitro.RData")
young_old_invitro <- as.data.frame(kegg.results.young_old_invitro$pathway.results)
young_old_invitro.o <- young_old_invitro[order(rownames(young_old_invitro)),]
young_old_invitro.o$name <- rownames(young_old_invitro.o)

load("~/Dropbox/Code_checking_SM/For_Brittany/Tables/KATJA_Pathway_analysis_Fib_young_old_in_vivo/kegg.results.pathways.kegg.results.young_old_invivo.RData")
young_old_invivo <- as.data.frame(kegg.results.young_old_invivo$pathway.results)
young_old_invivo.o <- young_old_invivo[order(rownames(young_old_invivo)),]
young_old_invivo.o$name <- rownames(young_old_invivo.o)


# merge all dataframes together
all.wound.data <- merge(young_basal_wounds.o, 
                        old_basal_wounds.o, by= "name")
all.basal.data <- merge(young_old_invitro.o, 
                        young_old_invivo.o,  by= "name")
all.data <- merge(all.wound.data, 
                  all.basal.data, by= "name")
row.names(all.data) <- row.names(young_old_invivo.o)

# limit analysis to pathways that are at least significant in one condition, subset data
all.data.subset <- all.data[,grep("adj.p.val", colnames(all.data))]
row.names(all.data.subset) <- row.names(young_old_invivo.o)
all.data.keep.id <- all.data.subset[rowSums(all.data.subset < 0.05) >= 2,]

# subset data.frame based on this list
all.data.sig <- all.data[which(row.names(all.data) %in% row.names(all.data.keep.id)),]
colnames(all.data.sig)

# subset data to only include t-stats
heatmap.data <- as.data.frame(c(all.data.sig[,c(5,10)],-all.data.sig[,c(15,20)]))
row.names(heatmap.data) <- row.names(all.data.sig)
colnames(heatmap.data) <- c("Young basal vs wound", "Old basal vs wound",
                            "Young vs old in vitro", "Young vs old in vivo")

# Generate heatmap based on the t-stat values
pdf("Heatmap_common_pathways_fibs_aging_wounding_invivo_invitro_092818.pdf", width =4.3, height =5, onefile=F  )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
pheatmap(heatmap.data,
         border_color = NA,
         cluster_col=T, 
         cluster_row=T, 
         scale="row", 
         annotation_legend=TRUE,
         clustering_distance_cols="manhattan", 
         clustering_method="complete",
         fontsize_col=8,
         fontsize_row=3.5,
         show_rownames=T,
         show_colnames=T,
         fontface="bold",
         treeheight_col=10,
         treeheight_row=0,
         cellwidth = NA, 
         cellheight = NA,
         color = colorRampPalette(c("darkblue","dodgerblue4", "white","brown1","darkred"))(100))
dev.off()
