# Objective: Determine the proportion of young and old cells in the identified Seurat clusters
# Use clustering data from Seurat step2

# Generate a barplot to visualizes the differences in subpopulations with age
####################################################################################
# Required packages 
library("Seurat")

###################################################################################
####################################################################################
setwd("~/Dropbox/Code_checking_SM/For_Matt/Seurat/")

# load Seurat object with clustering info
fibroblast_combined2 <- readRDS("~/Dropbox/Code_checking_SM/For_Matt/Seurat/fibroblast_combined2.rds")
fibroblast_combined2 <- readRDS(file = "fibroblast_combined2.rds")

# tSNE aging information
p1 <- TSNEPlot(fibroblast_combined2, 
               do.return = T, 
               pt.size = 0.3,
               group.by = "old",
               colors.use = c("darkblue","orange"))


# tSNE cluster information
p2 <- TSNEPlot(fibroblast_combined2, 
               do.label = T, 
               do.return = T, 
               pt.size = 0.3)


# Calculate proportion of cells in cluster 1 versus cluster 2 
# make a data fram with the required information
test <- data.frame(matrix(ncol = length(fibroblast_combined2@ident),nrow = 2))
test[1,] <- p1$data$ident
test[2,] <- p2$data$ident
head(test)

# Subset old cells in the 2 identified clusters
old <- as.data.frame(test[,grep("Old",test[1,])]) # 1444
Cluster1 <- as.data.frame(old[,which(old[2,] == "0")]) # 784
Cluster2 <- as.data.frame(old[,which(old[2,] == "1")]) # 660

# Proportion of old cells in cluster 1 and 2
o1 <- length(Cluster1[,grep("Old", Cluster1[1,])]) / length(test[,grep("Old", test[1,])]) # 0.5429363 
o2 <- length(Cluster2[,grep("Old", Cluster2[1,])]) / length(test[,grep("Old", test[1,])]) # 0.4570637

# subset young cells in the 2 identified clusters
young <- as.data.frame(test[,grep("Young",test[1,])]) # 1592
Cluster1 <- as.data.frame(young[,which(young[2,] == "0")]) # 1090
Cluster2 <- as.data.frame(young[,which(young[2,] == "1")]) # 502

# Proportion of young in cluster 1 and 2
y1 <- length(Cluster1[,grep("Young", Cluster1[1,])]) / length(test[,grep("Young", test[1,])]) # 0.6846734 
y2 <- length(Cluster2[,grep("Young", Cluster2[1,])]) / length(test[,grep("Young", test[1,])]) # 0.3153266

# make a table with fold changes in populations with age
young.pop <- c(y1,y2)
old.pop <- c(o1,o2)
fc.change <- old.pop/young.pop
log2fc.change <- log2(old.pop/young.pop)

# make a dataframe for the data
table.pop <- data.frame(old.pop,young.pop,fc.change,log2fc.change)                                

# make a barplot to visualize the results
pdf("Barplot_changes_in_cell_number_with_age_Seurat_clustering_2clusters.pdf", width =4.3, height =5, onefile=F  )
barplot(table.pop$log2fc.change,
        ylim = c(-0.6,0.6),
        ylab = "Log2 Enrichment",
        main = "Change in cell number with age",
        col = c("#F8766D","#00BFC4"))
dev.off()


