# Objective: To analyze RNAseq (EMTAB3037) of human primary fibroblast of donors from different ages using the DESeq2 package
# Determine DEGs and derive VST values for clustering, PCA, heatmaps

# Upload HTseq counts data of tophat2 mapped files (made on the cluster).
# Remove genes that have across all samples less than 0.3 fpkm
# Normalize and aqcuire vst value using the DESeq2 package, remove batch effect using limma package

###########################################################################################
# Need the following packages
library('edgeR')    # v3.22.3 
library("DESeq2")   # v1.20.0
sessionInfo()
######################################################################################################################################
###############################    removing low expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Human_healthy/')

# Upload all HTseq count files and merge fibroblast data into a data frame
HTseqpath="~/Dropbox/Code_checking_SM/For_Brittany/Essentials/HTSeq_reads/Human_aging_data/"
list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=TRUE)
HTseqfiles <-list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=T)
HTseq <- do.call("cbind", lapply(HTseqfiles, read.csv, header = TRUE, row.names=1, sep='\t'))
HTseqfiles_names <- list.files(".txt",path=HTseqpath,recursive=TRUE, full.names=F)
colnames(HTseq) <- HTseqfiles_names

smCountTable.1a <- HTseq[,grep("ERR", colnames(HTseq))]
colnames(smCountTable.1a)
smCountTable.1a <- smCountTable.1a[1:26363,]
tail(smCountTable.1a) 

# Rename the files so that it is not so long
colnames(smCountTable.1a) <- c("HF_69Y","HF_1Y", "HF_0Y","HF_30Y",
                               "HF_0Y","HF_74Y","HF_43Y","HF_29Y",
                               "HF_45Y","HF_87Y", "HF_71Y", "HF_89Y",
                               "HF_1Y","HF_36Y", "HF_65Y", "HF_66Y", 
                               "HF_45Y","HF_50Y", "HF_55Y", "HF_29Y")
colnames(smCountTable.1a)

# Select for only the heathy samples based on supplementary table
smCountTable.1.healthy <- smCountTable.1a[, -grep("HF_74Y|HF_65Y|HF_55Y|HF_50Y|HF_45Y|HF_36Y", colnames(smCountTable.1a))]
colnames(smCountTable.1.healthy)

## Remove all genes where at least 1 samples have more than 0.3 FPKM
# Upload gene length info
load("~/Dropbox/RNAseq/Pathway_terms/gene_length_hg19.RData")
exonic.gene.sizes.table <- exonic.gene.sizes.table[2:26364,]

## Remove all genes where at least 1 samples have more than 0.3 FPKM
# # extract the genes that exist in your list that you want the gene length for
ind <- match(rownames(smCountTable.1.healthy), names(exonic.gene.sizes.table))
exonic.gene.sizes.ord <- exonic.gene.sizes.table[1]
length <- data.frame(gene.symbol=rownames(smCountTable.1.healthy),
                     exonic.size=exonic.gene.sizes.ord)
length1<- as.numeric(length[,2])

# Calculate fpkm using edgeR function, extract all genes that have at least >0.3 fpkm in 1 sample
rpkm <- rpkm(smCountTable.1.healthy, gene.length=exonic.gene.sizes.table, normalized.lib.sizes=F, log=FALSE)
keep <- rowSums(rpkm > 0.3) >= 1
smCountTable.1.healthy.filt <- smCountTable.1.healthy[keep,] # maintain 14933 genes for further analysis

# specifiy batch and condtion, (for some reason had to add 0.001 to make them numeric)
Age = c(69.001, 1.001, 0.001, 30.001,	
        0.001, 43.001, 29.001,	
        87.001, 71.001, 89.001,
        1.001, 66.001, 29.001)

# Make table with counts, condtion and batch and store data in DESeqDataSet (dds). Run DESeq2 
# Make sure batch effect comes before condition in formula
exp.data.healthy <- data.frame(
                    row.names=colnames(smCountTable.1.healthy.filt),
                    Age = Age)
exp.data.healthy
dds <- DESeqDataSetFromMatrix(countData = smCountTable.1.healthy.filt,
                              colData = exp.data.healthy,
                              design = ~Age )
dds
dds <- DESeq(dds)
##########################################
# First look at the different varibles included in the model
resultsNames(dds)
# [1] "Intercept" "Age" 

# Choose the variable you are interested in
res <- results(dds, name = "Age", cooksCutoff=T, independentFiltering = T)

# QC check for distribution of p-values
hist(res$pvalue) # Looks good!

# order by FDR
resOrdered <- res[order(res$pvalue),]
# Keep all genes with FDR < 0.05
resSig <- subset(resOrdered, padj < 0.05) 
summary(resSig)
# out of 323 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 83, 26%
# LFC < 0 (down)     : 240, 74%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 2)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# Wont filter based on fold change as this is linear modeling

# write tables
write.table(resOrdered, file="DESeq2_human_Fibs_healthy_all_genes_linear_model.csv", sep=',')
save(resOrdered, file = "DESeq2_DEG_human_Fibs_healthy_all_genes_linear_model.Rda")

###############################    Extracting batch-corrected vst-values ----------------------------------------------------------
## Extract the log2 normalized variance stabilized values (vst) and remove batch-effect by removeBatchEffect from limma
vst <- varianceStabilizingTransformation(dds, blind =F)
vstMat <- assay(vst)

# write vst-vlaues to table 
write.table(vstMat, file="DESeq2_vst_human_Fibs_healthy_all_genes_linear_model.csv", sep=',')
save(vstMat, file = "DESeq2_vst_human_Fibs_healthy_all_genes_linear_model.Rda")

