# Objective: Determine the proportion of young and old cells in the identified PAGODA clusters
# Use clustering data from PAGODA Results_KEGG_disease_aging_correctingCELLCYCLE

# Generate a barplot to visualizes the differences in subpopulations with age
################################################################################################################
# packages required

###############################################################################################################
# Load files required 
setwd("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/Porportions_young_old/Expected_data/")

# PAGODA clusters
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/Expected_Data/cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.repeat_set_seed.RData")

# Expression data with age information
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/code_checking/expression_data/young_old_raw_counts_gene5500_mito10perc.Rdata")

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


###################################################################################################################
###################################################################################################################
# Determine statistically if there is a general difference in populations with age using chi-square test

# First normalize the number of young cells to match the old ones (as there are less old than young cells) 
y1.norm <- as.integer((length(Cluster_y1[,grep("Young", Cluster_y1[1,])])/ncol(young)) * ncol(old))
y2.norm <- as.integer((length(Cluster_y2[,grep("Young", Cluster_y2[1,])])/ncol(young)) * ncol(old))
y3.norm <- as.integer((length(Cluster_y3[,grep("Young", Cluster_y3[1,])])/ncol(young)) * ncol(old))
y4.norm <- as.integer((length(Cluster_y4[,grep("Young", Cluster_y4[1,])])/ncol(young)) * ncol(old))
y5.norm <- as.integer((length(Cluster_y5[,grep("Young", Cluster_y5[1,])])/ncol(young)) * ncol(old))
y6.norm <- as.integer((length(Cluster_y6[,grep("Young", Cluster_y6[1,])])/ncol(young)) * ncol(old))

# Make a data frame with all data for the test
chisq.data <- as.data.frame(matrix(nrow = 2, ncol = 6))

chisq.data[,1] <- c(y1.norm,length(Cluster_o1[,grep("Old", Cluster_o1[1,])]))
chisq.data[,2] <- c(y2.norm,length(Cluster_o2[,grep("Old", Cluster_o2[1,])]))
chisq.data[,3] <- c(y3.norm,length(Cluster_o3[,grep("Old", Cluster_o3[1,])]))
chisq.data[,4] <- c(y4.norm,length(Cluster_o4[,grep("Old", Cluster_o4[1,])]))
chisq.data[,5] <- c(y5.norm,length(Cluster_o5[,grep("Old", Cluster_o5[1,])]))
chisq.data[,6] <- c(y6.norm,length(Cluster_o6[,grep("Old", Cluster_o6[1,])]))
chisq.data

# Add row and col names
row.names(chisq.data) <- c("young", "old")
colnames(chisq.data) <- c("Sg1","Sg2","Sg3","Sg4","Sg5","Sg6")

# Run Chi-square test
chisq <- chisq.test(chisq.data)
chisq
# data:  chisq.data
# X-squared = 129.42, df = 5, p-value < 2.2e-16

# Calculate for each population, if they are more different in old than expected by chance
# total number of cells subgroup 1
sb1 <- y1.norm + length(Cluster_o1[,grep("Old", Cluster_o1[1,])])
res.sb1 <- prop.test(x = length(Cluster_o1[,grep("Old", Cluster_o1[1,])]), 
                     n = sb1, 
                     p = 0.5, correct = FALSE,
                 alternative = "two.sided") 
res.sb1$p.value
# p-value = 0.1653272

sb2 <- y2.norm + length(Cluster_o2[,grep("Old", Cluster_o2[1,])])
res.sb2 <- prop.test(x = length(Cluster_o2[,grep("Old", Cluster_o2[1,])]), 
                     n = sb2, 
                     p = 0.5, correct = FALSE,
                     alternative = "two.sided") 
res.sb2$p.value # 0.6791157

sb3 <- y3.norm + length(Cluster_o3[,grep("Old", Cluster_o3[1,])])
res.sb3 <- prop.test(x = length(Cluster_o3[,grep("Old", Cluster_o3[1,])]), 
                     n = sb3, 
                     p = 0.5, correct = FALSE,
                     alternative = "two.sided") 
res.sb3$p.value # 2.159395e-05

sb4 <- y4.norm + length(Cluster_o4[,grep("Old", Cluster_o4[1,])])
res.sb4 <- prop.test(x = length(Cluster_o4[,grep("Old", Cluster_o4[1,])]), 
                     n = sb4, 
                     p = 0.5, correct = FALSE,
                     alternative = "two.sided") 
res.sb4$p.value # 1.968748e-14

sb5 <- y5.norm + length(Cluster_o5[,grep("Old", Cluster_o5[1,])])
res.sb5 <- prop.test(x = length(Cluster_o5[,grep("Old", Cluster_o5[1,])]), 
                     n = sb5, 
                     p = 0.5, correct = FALSE,
                     alternative = "two.sided") 
res.sb5$p.value # 8.56908e-08

sb6 <- y6.norm + length(Cluster_o6[,grep("Old", Cluster_o6[1,])])
res.sb6 <- prop.test(x = length(Cluster_o6[,grep("Old", Cluster_o6[1,])]), 
                     n = sb6, 
                     p = 0.5, correct = FALSE,
                     alternative = "two.sided") 
res.sb6$p.value # 2.962167e-06

table.pop.order <- table.pop[c(3,5,2,4,6,1),]

# Correct for FDR
p.adjust(c(res.sb3$p.value,res.sb5$p.value,res.sb2$p.value,res.sb4$p.value,res.sb6$p.value,res.sb1$p.value), method = "BH")
# [1] 3.239092e-05 2.570724e-07 6.791157e-01 1.181249e-13 5.924334e-06 1.983926e-01