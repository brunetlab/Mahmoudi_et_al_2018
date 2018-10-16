# Objective: Determine the proportion of young and old cells in the identified PAGODA clusters
# Use clustering data from PAGODA Results_KEGG_disease_aging_correctingCELLCYCLE

# Generate a barplot to visualizes the differences in subpopulations with age
################################################################################################################
# packages required

###############################################################################################################
# Load files required 
setwd("/PAGODA/Porportions_young_old/")

# PAGODA clusters
load("/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/Expected_Data/cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.repeat_set_seed.RData")

# Expression data with age information
load("/PAGODA/code_checking/expression_data/young_old_raw_counts_gene5500_mito10perc.Rdata")

################################################################################################################
# Define the different subpopulations (6 main ones based on PAGODA)
sg <- as.factor(cutree(hc, k=6))
names(sg) <- hc$labels

# Calculate proportion of young and old cells in the different clusters 
# make a data fram with necessary information
test <- data.frame(matrix(ncol = ncol(young_old_raw_counts_gene5500_mito10perc),nrow = 2))
test[1,] <-  c(rep("Young",length(grep("Young", colnames(young_old_raw_counts_gene5500_mito10perc)))),
               rep("Old",length(grep("Old", colnames(young_old_raw_counts_gene5500_mito10perc)))))
test[2,] <- sg

# subset old cells in the 4 identified clusters
old <- as.data.frame(test[,grep("Old",test[1,])]) # 1444
Cluster_o1 <- as.data.frame(old[,which(old[2,] == "1")]) # 250
Cluster_o2 <- as.data.frame(old[,which(old[2,] == "2")]) # 191
Cluster_o3 <- as.data.frame(old[,which(old[2,] == "3")]) # 264
Cluster_o4 <- as.data.frame(old[,which(old[2,] == "4")]) # 354
Cluster_o5 <- as.data.frame(old[,which(old[2,] == "5")]) # 266
Cluster_o6 <- as.data.frame(old[,which(old[2,] == "6")]) # 119

# Proportion (%) of old in clusters
o1 <- length(Cluster_o1[,grep("Old", Cluster_o1[1,])]) / length(test[,grep("Old", test[1,])]) # 0.1731302
o2 <- length(Cluster_o2[,grep("Old", Cluster_o2[1,])]) / length(test[,grep("Old", test[1,])]) # 0.1322715
o3 <- length(Cluster_o3[,grep("Old", Cluster_o3[1,])]) / length(test[,grep("Old", test[1,])]) # 0.1828255
o4 <- length(Cluster_o4[,grep("Old", Cluster_o4[1,])]) / length(test[,grep("Old", test[1,])]) # 0.2451524
o5 <- length(Cluster_o5[,grep("Old", Cluster_o5[1,])]) / length(test[,grep("Old", test[1,])]) # 0.1842105
o6 <- length(Cluster_o6[,grep("Old", Cluster_o6[1,])]) / length(test[,grep("Old", test[1,])]) # 0.08240997

# subset young cells in the 4 identified clusters
young <- as.data.frame(test[,grep("Young",test[1,])]) # 1592
Cluster_y1 <- as.data.frame(young[,which(young[2,] == "1")]) # 312
Cluster_y2 <- as.data.frame(young[,which(young[2,] == "2")]) # 202
Cluster_y3 <- as.data.frame(young[,which(young[2,] == "3")]) # 193
Cluster_y4 <- as.data.frame(young[,which(young[2,] == "4")]) # 650
Cluster_y5 <- as.data.frame(young[,which(young[2,] == "5")]) # 172
Cluster_y6 <- as.data.frame(young[,which(young[2,] == "6")]) # 63

# Proportion of young in clusters
y1 <- length(Cluster_y1[,grep("Young", Cluster_y1[1,])]) / length(test[,grep("Young", test[1,])]) # 0.1959799
y2 <- length(Cluster_y2[,grep("Young", Cluster_y2[1,])]) / length(test[,grep("Young", test[1,])]) # 0.1268844
y3 <- length(Cluster_y3[,grep("Young", Cluster_y3[1,])]) / length(test[,grep("Young", test[1,])]) # 0.1212312
y4 <- length(Cluster_y4[,grep("Young", Cluster_y4[1,])]) / length(test[,grep("Young", test[1,])]) # 0.4082915
y5 <- length(Cluster_y5[,grep("Young", Cluster_y5[1,])]) / length(test[,grep("Young", test[1,])]) # 0.1080402
y6 <- length(Cluster_y6[,grep("Young", Cluster_y6[1,])]) / length(test[,grep("Young", test[1,])]) # 0.03957286

# make a table with fold changes in populations with age
young.pop <- c(y1,y2,y3,y4,y5,y6)
old.pop <- c(o1,o2,o3,o4,o5,o6)
fc.change <- old.pop/young.pop
log2fc.change <- log2(old.pop/young.pop)

# make a dataframe for the data
table.pop <- data.frame(old.pop,young.pop,fc.change,log2fc.change)                                
# reorder so that it comes in the order of the PAGODA plot
table.pop.order <- table.pop[c(3,5,2,4,6,1),]

# make a barplot to visualize the data
pdf("Barplot_changes_in_cell_number_with_age_PAGODA_clustering_6main_clusters.pdf", width =4.3, height =5, onefile=F  )
barplot(table.pop.order$log2fc.change,
        ylim = c(-1,1),
        ylab = "Log2 Enrichment",
        main = "Change in cell number with age")
dev.off()
