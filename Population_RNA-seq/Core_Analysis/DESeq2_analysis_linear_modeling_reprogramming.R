# Objective: Regression modeling of young and old fibroblast transcriptomes to identify genes that correlate with reprogramming using DEseq2

# Load HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have less than 0.3 fpkm across all samples 
# Include in the model reprogramming efficiency 
# Determine the genes that show correlation with RE

# Output:
# 1) Table with all DEGs (padjust < 0.05, abs(LogFC) > 0.58) -> Extended Data Table
# 2) .Rdata table with all genes for DEG analysis -> used for different pathway enrichment analyses
# 3) .RData table with batch-corrected VST values -> used for PCA, heatmaps

###########################################################################################
# Need the following packages
library('edgeR')    # v3.10.2 
library("DESeq2")   # v1.20.0
library('limma')    # v3.24.15

######################################################################################################################################
###############################    removing low expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro_linear_modeling_RE/')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany//Essentials/HTSeq_reads/Fibroblast_in_vitro/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1a <- HTseq[,grep("Fib", colnames(HTseq))]
colnames(smCountTable.1a)

# reorder so that young (3m) comes first, and remove Fib_29M6C11L_BAD_htseq_count.txt and Fib_29M1C10R_GOOD_htseq_count.txt
smCountTable.1 <- smCountTable.1a[,c(grep("Fib_3", colnames(smCountTable.1a)),
                                     grep("Fib_29", colnames(smCountTable.1a)))]
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

RE = c(1.02, 0.83,1.26,
       1.00, 1.01,0.74,
       0.62, 1.10,0.66,
       1.56, 0.25,1.33,
       0.55, 1.61,0.19,
       2.79, 1.50,0.73)

# Make table with counts, batch and RE, and store data in DESeqDataSet (dds). Run DESeq2 
exp.data <- data.frame(
  row.names=colnames(smCountTable.filt),
  condition = condition,
  batch = batch,
  RE = RE)

dds <- DESeqDataSetFromMatrix(countData = smCountTable.filt,
                              colData = exp.data,
                              design = ~ batch + RE )
dds
dds <- DESeq(dds)

###############################   Differential Gene Expression Analysis  ----------------------------------------------------------
# First look at the different varibles included in the model
resultsNames(dds)
# "Intercept"      "batch_B2_vs_B1" "batch_B3_vs_B1" "RE"

# Choose the variable you are interested in
res <- results(dds, name = "RE", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!


# order by FDR
resOrdered <- res[order(res$ padj),]
summary(resOrdered)
# out of 14310 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 0, 0% 
# LFC < 0 (down)   : 0, 0% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 0, 0% 
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05)
summary(resSig)
# out of 0 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 0, NaN% 
# LFC < 0 (down)   : 0, NaN% 
# outliers [1]     : 0, NaN% 
# low counts [2]   : 0, NaN% 
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write tables
write.table(resOrdered, file="DESeq2_Linear_modeling_RE_young_old.csv", sep=',')

# save DEG table 
save(resOrdered, file = "DESeq2_Linear_modeling_RE_young_old.Rda")

###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)

# batch-correct using limma function
v.data <- removeBatchEffect(vstMat,batch=batch)

# save VST table 
save(v.data,file="DESeq2_Linear_modeling_RE_young_old_batch_corrected_vst.Rda")

