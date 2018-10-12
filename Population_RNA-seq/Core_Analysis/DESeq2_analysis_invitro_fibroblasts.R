# Objective: To analyze transcriptome of young and old fibroblasts in culture using DESeq2
# Determine DEGs and derive batch-corrected VST values for clustering, PCA, heatmaps

# Upload HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have across all samples less than 0.3 fpkm
# Normalize and aqcuire vst value using the DESeq2 package, remove batch effect using limma package
# 
# Output:
# 1) Table with all DEGs (padjust < 0.05, abs(LogFC) > 0.58) -> Extended Data Table
# 2) .Rdata table with all genes for DEG analysis -> used for different pathway enrichment analyses
# 3) .RData table with batch-corrected VST values -> used for PCA, heatmaps
###########################################################################################
# Required packages
library('edgeR')    # v3.12.1 
library("DESeq2")   # v1.20.0
library('limma')    # v3.26.9
sessionInfo()
######################################################################################################################################
###############################    removing lowly expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro/')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Fibroblast_in_vitro/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1a <- HTseq[,grep("Fib", colnames(HTseq))]
colnames(smCountTable.1a)

# reorder so that young (3m) comes first
smCountTable.1 <- cbind(smCountTable.1a[,grep("Fib_3", colnames(smCountTable.1a))], smCountTable.1a[,grep("Fib_29", colnames(smCountTable.1a))])
colnames(smCountTable.1)

## Remove all genes where at least 1 samples have more than 0.3 FPKM
# Upload gene length info
load("~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Gene_sizes/exonic_gene_sizes_mm9.RData")

# extract the genes that exist in your list that you want the gene length for
ind <- match(rownames(smCountTable.1), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
length <- data.frame(gene.symbol=rownames(smCountTable.1),
                     exonic.size=exonic.gene.sizes.ord)
length1<- as.numeric(length[,2])

# Calculate fpkm using edgeR function, extract all genes that have at least >0.3 fpkm in 1 sample
fpkm <- rpkm(smCountTable.1, gene.length=length1,normalized.lib.sizes=F, log=FALSE)
keep <- rowSums(fpkm > 0.3) >= 1
smCountTable.filt <- smCountTable.1[keep,] # maintain 14310 genes for further analysis

# specify batch and condtion and RE
condition = c(rep('Fib_young',8), rep('Fib_Old',10))
batch =factor(c("B1","B1","B3",
                "B2","B1","B3",
                "B2","B2","B1",
                "B2","B3","B3",
                "B2","B3","B3",
                "B1","B1","B2"))

# Make table with counts, condtion and batch, and store data in DESeqDataSet (dds). Run DESeq2 
exp.data <- data.frame(
  row.names=colnames(smCountTable.filt),
  condition = condition,
  batch = batch)

dds <- DESeqDataSetFromMatrix(countData = smCountTable.filt,
                              colData = exp.data,
                              design = ~ batch + condition )
dds
dds <- DESeq(dds)

###############################   Differential Gene Expression Analysis  ----------------------------------------------------------
# First look at the different varibles included in the model
resultsNames(dds)
# [1] "Intercept"     "batch_B2_vs_B1"   "batch_B3_vs_B1"   "condition_Fib_young_vs_Fib_Old"

# Choose the variable you are interested in
res <- results(dds, name = "condition_Fib_young_vs_Fib_Old", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!

res <- results(dds, cooksCutoff=T, independentFiltering = T)
# order by FDR
resOrdered <- res[order(res$ padj),]
summary(resOrdered)
# out of 14310 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1433, 10%
# LFC < 0 (down)     : 1833, 13%
# outliers [1]       : 4, 0.028%
# low counts [2]     : 278, 1.9%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05) # 2563 DEGs at FDR < 0.05
summary(resSig)
# out of 2565 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 1103, 43%
# LFC < 0 (down)     : 1462, 57%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with a log2 fold change of >0.58 (corresponds to 1.5 fold change)
resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.58) # 885 DEGs at FDR < 0.05 and 1.5 FC
summary(resSIG_FC)
# out of 1022 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 389, 38%
# LFC < 0 (down)     : 633, 62%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write tables
write.table(resSIG_FC, file="DESeq2_DEGs_1_5FC.csv", sep=',')
write.table(resOrdered, file="DESeq2_all_genes.csv", sep=',')

# save DEG table
save(resOrdered, file = "DESeq2_DEG_all_genes.Rda")
###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)

# batch-correct using limma function
v.data <- removeBatchEffect(vstMat,batch=batch)

# save VST table
save(v.data, file = "DESeq2_vst_corrected_values.Rda")


