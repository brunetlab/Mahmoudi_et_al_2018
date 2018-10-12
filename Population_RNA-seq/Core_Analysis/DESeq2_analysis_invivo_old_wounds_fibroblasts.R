# Objective: To compare transcriptomes of in vivo FACS sorted old THY1+ and THY1- fibroblasts before and wounds 7days after recovery using the DESeq2 package
# Determine DEGs between basal and wounded THY1+ and THY1- cells and derive batch-corrected VST values for clustering, PCA, heatmaps

# Upload HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have across all samples less than 0.3 fpkm
# Take into account batch and celltype (if they are THY1+ or THY1-, as we are only interested in the aging aspects here)
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
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vivo_old_wounds//')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Fibroblast_in_vivo/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

# Restrict analysis to the same bacth 
smCountTable.1a <- HTseq[,grep("Fib_O1|Fib_O2|Fib_WO", colnames(HTseq))]
colnames(smCountTable.1a)

# reorder so that young (3m) comes first
# smCountTable.1a <- smCountTable.1a[,c(grep("Fib_Y", colnames(smCountTable.1a)),
#                                       grep("Fib_WY", colnames(smCountTable.1a)))]
colnames(smCountTable.1a)

# Remove the last 5 lines as these are not genes
smCountTable.1 <- smCountTable.1a[1:23365,]

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
smCountTable.filt <- smCountTable.1[keep,] # maintain 15961 genes for further analysis

# specify condtion, batch and cell type (if they are THY1pos or THY1neg)
condition = c(rep("basal",4),
              rep("wound",6))
# batch =factor(c("B1","B1","B1","B1","B1","B1",
#                 "B2","B2","B2","B2","B2",
#                 "B1","B1","B1","B1",
#                 "B2","B2","B2","B2","B2","B2"))
celltype <- factor(c("Neg","Pos","Neg","Pos",
                     "Neg","Pos","Neg","Pos","Neg","Pos"))

# Make table with counts, condtion, batch and celltype and store data in DESeqDataSet (dds). Run DESeq2 
exp.data <- data.frame(
  row.names=colnames(smCountTable.filt),
  condition = condition,
  celltype = celltype)
exp.data
dds <- DESeqDataSetFromMatrix(countData = smCountTable.filt,
                              colData = exp.data,
                              design = ~ celltype + condition)
dds
dds <- DESeq(dds)

###############################   Differential Gene Expression Analysis  ----------------------------------------------------------
# First look at the different varibles included in the model
resultsNames(dds)
# "Intercept"              "celltype_Pos_vs_Neg"    "condition_wound_vs_basal"

# Choose the variable you are interested in
res <- results(dds, name = "condition_wound_vs_basal", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!

# order by FDR
resOrdered <- res[order(res$padj),]
# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05) 
summary(resSig)
# out of 3981 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2235, 56%
# LFC < 0 (down)     : 1746, 44%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with a log2 fold change of >0.58 (corresponds to 1.5 fold change)
resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.58) 
summary(resSIG_FC)
# out of 3957 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 2219, 56%
# LFC < 0 (down)     : 1738, 44%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write tables
write.table(resSIG_FC, file="DESeq2_old_wound_vs_basal_DEGs_1_5FC.csv", sep=',')
write.table(resOrdered, file="DESeq2_old_wound_vs_basal_genes.csv", sep=',')

# save DEG table 
save(resOrdered, file = "DESeq2_DEG_old_wound_vs_basal_all_genes.Rda")

###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)

# # batch-correct using limma function
# v.data <- removeBatchEffect(vstMat,batch=batch)

# save VST table 
save(vstMat, file = "DESeq2_young_old_in_vivo_vst_values.Rda")

