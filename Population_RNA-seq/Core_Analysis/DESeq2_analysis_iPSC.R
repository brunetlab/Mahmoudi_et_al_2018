# Objective: To analyze transcriptome of young and old derived iPSCs using the DESeq2 package
# Determine DEGs and derive VST values for clustering, PCA, heatmaps

# Upload HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have across all samples less than 0.3 fpkm
# Normalize and aqcuire vst value using the DESeq2 package, remove batch effect using limma package
# 
# Output:
# 1) Table with all DEGs (padjust < 0.05, abs(LogFC) > 0.58) -> Extended Data Table
# 2) .Rdata table with all genes for DEG analysis -> used for different pathway enrichment analyses
# 3) .RData table with batch-corrected VST values -> used for PCA, heatmaps
###########################################################################################
# Need the following packages
library('edgeR')    # v3.12.1 
library("DESeq2")   # v1.20.0
library('limma')    # v3.26.9
sessionInfo()
######################################################################################################################################
###############################    removing low expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/iPSC_in_vitro/')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Fibroblast_iPSC/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1a <- HTseq[,grep("iPSC", colnames(HTseq))]
colnames(smCountTable.1a)

# reorder so that young (3m) comes first
smCountTable.1 <- smCountTable.1a[,c(grep("_3m", colnames(smCountTable.1a)),grep("_29m", colnames(smCountTable.1a)) )]
colnames(smCountTable.1)

# Rename the files so that it is not so long
colnames(smCountTable.1) <- c("Young iPSC 1", "Young iPSC 2","Young iPSC 3","Young iPSC 4","Young iPSC 5",
                              "Old iPSC 1",  "Old iPSC 2",  "Old iPSC 3",  "Old iPSC 4",  "Old iPSC 5",  "Old iPSC 6")
colnames(smCountTable.1)

# Remove the last 5 lines as these are not genes
smCountTable.1 <- smCountTable.1[1:23365,]

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
smCountTable.filt <- smCountTable.1[keep,] # maintain 14043 genes for further analysis

# specifiy batch and condtion
condition = c('iPSC_A3m','iPSC_A3m','iPSC_A3m','iPSC_A3m','iPSC_A3m',
              'iPSC_B_29m','iPSC_B_29m','iPSC_B_29m','iPSC_B_29m','iPSC_B_29m','iPSC_B_29m')

# Make table with counts, condtion and batch and store data in DESeqDataSet (dds). Run DESeq2 
exp.data <- data.frame(
  row.names=colnames(smCountTable.filt),
  condition = condition)
dds <- DESeqDataSetFromMatrix(countData = smCountTable.filt,
                              colData = exp.data,
                              design = ~ condition)
dds
dds <- DESeq(dds)
###############################   Differential Gene Expression Analysis  ----------------------------------------------------------
# First look at the different varibles included in the model
resultsNames(dds)
#[1] "Intercept"                        "condition_iPSC_B_29m_vs_iPSC_A3m"

# Choose the variable you are interested in
res <- results(dds, name = "condition_iPSC_B_29m_vs_iPSC_A3m", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!

# order by FDR
resOrdered <- res[order(res$ padj),]
# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05) # 5 DEGs at FDR < 0.05
summary(resSig)
# out of 5 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 3, 60% 
# LFC < 0 (down)   : 2, 40% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with a log2 fold change of >0.58 (corresponds to 1.5 fold change)
resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.58) # 5 DEGs at FDR < 0.05 and 1.5 FC
summary(resSIG_FC)
head(resSIG_FC)
# out of 5 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3, 60%
# LFC < 0 (down)     : 2, 40%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


# write tables
write.table(resSIG_FC, file="DESeq2_iPSC_DEGs_1_5FC.csv", sep=',')
write.table(resOrdered, file="DESeq2_iPSC_all_genes.csv", sep=',')

# save DEG table 
save(resOrdered, file = "DESeq2_DEG_iPSC_all_genes.Rda")
###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
v.data <- assay(vst)

# save VST table 
save(v.data, file = "DESeq2_iPSC_vst_values.Rda")


