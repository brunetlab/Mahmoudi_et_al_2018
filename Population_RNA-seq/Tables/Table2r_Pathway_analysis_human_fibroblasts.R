# Objective: Perform Pathway analysis on DEGs (FDR < 0.05) using Fischer's exact test
# Figure 2h
# Use DEG table (with fold change and FDR-values) created by DEseq2  

# Will use KEGG pathways downloaded from GSEA Broad Institute

###########################################################################################
# Need the following package
sessionInfo()
###########################################################################################
# Set working directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Tables/Pathway_analysis_human_healthy_fib/')

# Load DEG data from DESeq2
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Human_healthy/DESeq2_DEG_human_Fibs_healthy_all_genes_linear_model.Rda")

# Make it into a data frame
genes<- as.data.frame(resOrdered)
# Make gene names in CAPTIAL letter
row.names(genes) <- toupper(row.names(genes))

# Upload kegg pathways (the .csv version) 
gene.sets <- read.table("~/Dropbox/Code_checking_SM/For_Brittany/Essentials/Pathway_genesets/KEGG_2016_disease_activation_UPPER_case.txt", sep="\t", header = F, quote="",row.names = 1)
gene.set.data <- as.data.frame(t(gene.sets))

########################################################################
# Extract genes that are DOWN with age
deg.down <- genes[which(genes$padj < 0.05 & genes$log2FoldChange < -0),] #240 DOWN with age

# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.down <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.down) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.down) <- colnames(gene.set.data)

# In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in the KEGG pathway are found in my gene list (background) (kegg_genes)
# 2) how many genes in the KEGG pathway are found DEG in my data set (deg_kegg_genes)
# 3) how many genes in the KEGG pathway are found in my gene list but are not DEG in my data set (deg_non_kegg_genes)
# 4) how many DEG genes are not found in the KEGG pathway (non_deg_kegg_genes)
# 5) how many genes are not DEG nor found in the KEGG pathway (non_deg_non_kegg_genes)
# Use the numbers to perform fischer's t-test, and save the p-value in the data frame 
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(toupper(rownames(genes)) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(toupper(rownames(deg.down)) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(deg.down[which(toupper(rownames(deg.down)) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(deg.down) - deg_kegg_genes 
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.down)) - kegg_genes 
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data.down$p.values[i] <- ftest$p.value
  my.data.down$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.down$FDR.BH <- p.adjust(my.data.down$p.values, "BH")
# Order by FDR adjusted p-values
my.data.downOrdered <- my.data.down[order(my.data.down$FDR.BH),]
head(my.data.downOrdered, n=10)

# Write a table with the final results 
write.csv(my.data.downOrdered, file="Kegg2016_pathways_DOWN_with_age_human.csv")

#####################################################################################################
# Extract genes that are UP with age
deg.up <- genes[which(genes$padj < 0.05 & genes$log2FoldChange > 0),] #83 UP with age
 
# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.up <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.up) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.up) <- colnames(gene.set.data)
 
# In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in the KEGG pathway are found in my gene list (background) (kegg_genes)
# 2) how many genes in the KEGG pathway are found DEG in my data set (deg_kegg_genes)
# 3) how many genes in the KEGG pathway are found in my gene list but are not DEG in my data set (deg_non_kegg_genes)
# 4) how many DEG genes are not found in the KEGG pathway (non_deg_kegg_genes)
# 5) how many genes are not DEG nor found in the KEGG pathway (non_deg_non_kegg_genes)
# Use the numbers to perform fischer's t-test, and save the p-value in the data frame 
# Finally correct for FDR using benjimini_hochberg
 
for(i in 1:ncol(gene.set.data)){
 kegg_genes <- nrow(genes[which(toupper(rownames(genes)) %in% gene.set.data[,i]),])
 deg_kegg_genes <- nrow(genes[which(toupper(rownames(deg.up)) %in% gene.set.data[,i]),])
 deg_kegg_genes.names <- as.data.frame(row.names(deg.up[which(toupper(rownames(deg.up)) %in% gene.set.data[,i]),]))
 deg_non_kegg_genes <- nrow(deg.up) - deg_kegg_genes 
 non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
 non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg.up)) - kegg_genes 
 ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
 my.data.up$p.values[i] <- ftest$p.value
 my.data.up$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.up$FDR.BH <- p.adjust(my.data.up$p.values, "BH")
# Order by FDR adjusted p-values
my.data.upOrdered <- my.data.up[order(my.data.up$FDR.BH),]
head(my.data.upOrdered, n=10)
 
# Write a table with the final results 
write.csv(my.data.upOrdered, file="Kegg2016_pathways_UP_with_age_human.csv")


#####################################################################################################
# Extract genes that are CHANGE with age
deg <- genes[which(genes$padj < 0.05),] #323 UP with age

# Make a data frame with 3 columns for p-value, FDR-corrected p-values , and the string of genes within the pathway
my.data.change <- as.data.frame(matrix(0, ncol = 3, nrow = ncol(gene.set.data)))
colnames(my.data.change) <- c("p.values", "FDR.BH", "genes")
rownames(my.data.change) <- colnames(gene.set.data)

# In a loop go though each and one of the KEGG pathways, and determine
# 1) how many genes in the KEGG pathway are found in my gene list (background) (kegg_genes)
# 2) how many genes in the KEGG pathway are found DEG in my data set (deg_kegg_genes)
# 3) how many genes in the KEGG pathway are found in my gene list but are not DEG in my data set (deg_non_kegg_genes)
# 4) how many DEG genes are not found in the KEGG pathway (non_deg_kegg_genes)
# 5) how many genes are not DEG nor found in the KEGG pathway (non_deg_non_kegg_genes)
# Use the numbers to perform fischer's t-test, and save the p-value in the data frame 
# Finally correct for FDR using benjimini_hochberg

for(i in 1:ncol(gene.set.data)){
  kegg_genes <- nrow(genes[which(toupper(rownames(genes)) %in% gene.set.data[,i]),])
  deg_kegg_genes <- nrow(genes[which(toupper(rownames(deg)) %in% gene.set.data[,i]),])
  deg_kegg_genes.names <- as.data.frame(row.names(deg[which(toupper(rownames(deg)) %in% gene.set.data[,i]),]))
  deg_non_kegg_genes <- nrow(deg) - deg_kegg_genes 
  non_deg_kegg_genes <- kegg_genes - deg_kegg_genes 
  non_deg_non_kegg_genes <- (nrow(genes) -nrow(deg)) - kegg_genes 
  ftest <- fisher.test(matrix(c(deg_kegg_genes,deg_non_kegg_genes,non_deg_kegg_genes,non_deg_non_kegg_genes),nrow=2,ncol=2),alternative="greater")
  my.data.change$p.values[i] <- ftest$p.value
  my.data.change$genes[i] <- paste(deg_kegg_genes.names[1:nrow(deg_kegg_genes.names),], collapse=", ")
}
my.data.change$FDR.BH <- p.adjust(my.data.change$p.values, "BH")
# Order by FDR adjusted p-values
my.data.changeOrdered <- my.data.change[order(my.data.change$FDR.BH),]
head(my.data.changeOrdered, n=10)

# Write a table with the final results 
write.csv(my.data.changeOrdered, file="Kegg2016_pathways_CHANGE_with_age_human.csv")
