# Objective: To compare the transcriptomes of the top 5 reprogrammers versus bottom 5 reporgammers in my datasets using DESeq2 

# Upload HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have across all samples less than 0.3 fpkm
# Normalize and aqcuire vst value using the DESeq2 package, remove batch effect using limma package

# Output:
# 1) Table with all DEGs (padjust < 0.05, abs(LogFC) > 0.58) -> Extended Data Table
# 2) .Rdata table with all genes for DEG analysis -> used for different pathway enrichment analyses
# 3) .RData table with batch-corrected VST values -> used for PCA, heatmaps

###########################################################################################
# Need the following packages
library('edgeR')    # v3.10.2 
library("DESeq2")   # v1.20.0
library('limma')    # v3.24.15
sessionInfo()
######################################################################################################################################
###############################    removing low expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro_Top5_vs_bottom5/')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Fibroblast_in_vitro_Top5_vs_bottom5/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1 <- HTseq[,grep("Fib", colnames(HTseq))]
colnames(smCountTable.1)

## Remove all genes where at least 1 samples have more than 0.3 FPKM
# Upload gene length info
load("~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Gene_sizes//exonic_gene_sizes_mm9.RData")

# extract the genes that exist in your list that you want the gene length for
ind <- match(rownames(smCountTable.1), names(exonic.gene.sizes))
exonic.gene.sizes.ord <- exonic.gene.sizes[ind]
length <- data.frame(gene.symbol=rownames(smCountTable.1),
                     exonic.size=exonic.gene.sizes.ord)
length1<- as.numeric(length[,2])

# Calculate fpkm using edgeR function, extract all genes that have at least >0.3 fpkm in 1 sample
fpkm <- rpkm(smCountTable.1, gene.length=length1,normalized.lib.sizes=F, log=FALSE)
keep <- rowSums(fpkm > 0.3) >= 1
smCountTable.filt <- smCountTable.1[keep,] # maintain 13919 genes for further analysis
colnames(smCountTable.filt)

# specify batch and condtion and RE
condition = c("bad", "bad", "bad", "bad","bad",
              "good", "good","good","good","good")

batch =factor(c("Exp3","Exp3","Exp2",
                "Exp2","Exp1","Exp3",
                "Exp2","Exp1","Exp3",
                "Exp1"))

# Make table with counts, condition and batch and store data in DESeqDataSet (dds). Run DESeq2 
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
# [1] "Intercept"             "batch_Exp2_vs_Exp1"    "batch_Exp3_vs_Exp1"    "condition_good_vs_bad"

# Choose the variable you are interested in
res <- results(dds, name = "condition_good_vs_bad", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!

# order by FDR
resOrdered <- res[order(res$ padj),]
summary(resOrdered)
# out of 13919 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3, 0.022%
# LFC < 0 (down)     : 4, 0.029%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05) # 5 DEGs at FDR < 0.05
summary(resSig)
# out of 6 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3, 50%
# LFC < 0 (down)     : 3, 50%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Keep all genes with a log2 fold change of >0.58 (corresponds to 1.5 fold change)
resSIG_FC <- subset(resSig, abs(log2FoldChange) > 0.58) # 0 DEGs at FDR < 0.05 and 1.5 FC
summary(resSIG_FC)
# out of 6 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 3, 50%
# LFC < 0 (down)     : 3, 50%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# write tables
write.table(resOrdered, file="DESeq2_TopBottom5_RE_young_old.csv", sep=',')

# save DEG table
save(resOrdered, file = "DESeq2_TopBottom5_RE_young_old_all_genes.Rda")

###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)

# batch-correct using limma function
v.data <- removeBatchEffect(vstMat,batch=batch)

# save VST table
save(v.data,file="DESeq2_TopBottom5_RE_young_old_batch_corrected_vst..Rda")

