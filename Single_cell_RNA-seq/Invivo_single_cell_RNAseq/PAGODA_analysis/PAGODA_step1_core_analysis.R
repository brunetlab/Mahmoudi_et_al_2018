# Objective: perform pagaoda analysis of single cell data of young and old fibroblasts 
# Fig. 4 and Extended Data Fig. 9 
# Use expression data filtered by Seurat using Seurat_filter_normalize_merge_single_cell_RNAseq_data_EDFig_9a.R script: young_old_raw_counts_gene5500_mito10perc.RData
# I follow the tutorial of pagoda on: http://hms-dbmi.github.io/scde/pagoda.html

################################################################################################################
# packages required
library(Biobase)
library(scde)
library(RColorBrewer)
################################################################################################################
load("Single_cell_RNA-seq/Invivo_single_cell_RNAseq/Essentials/Expression_data_PAGODA/young_old_raw_counts_gene5500_mito10perc.Rdata")

## Fitting error models
# Use k-nearest neighbor model fitting procedure to construct error models for individual cells
# with age-groups as cell types:

# Use the following parameters:
# k = 1/4 of most similar cells (estimating ~4 subpopulations)
# n.cores = 10 (number of cores to be used for analysis, this only works for cluster) 
# min.size.entries = 2000 (min number of genes to use for model fitting)
# min.count.threshold = 1 (min number of reads for the gene to be initially classified as a non-failed measurement),
# min.nonfailed = 5 (requiring at least 5 non-failed measurements per gene)
# max.model.plots = 50 (maximum number of models to save plots for (saves time when there are too many cells))

# # Make new directories and setwd accordingly
dir.create("PAGODA")
setwd("./PAGODA/")
dir.create("knn.error.models")
setwd("./knn.error.models/")

knn <- knn.error.models(young_old_raw_counts_gene5500_mito10perc,
                        k = round(ncol(young_old_raw_counts_gene5500_mito10perc)/4),
                        n.cores = 10,
                        min.size.entries = 2000,
                        min.count.threshold = 1,
                        min.nonfailed = 5,
                        max.model.plots = 50)

save(knn, file = "knn.error.models.RData")

# ################################################################################################################
# ## Normalizing variance
# # Normalize out expected levels of technical and intrinsic biological noise
# 
# # Make new directory
setwd("/PAGODA/")
dir.create("Varinfo")
setwd("./Varinfo/")

# Use the following parameters:
# trim = 3/ncol(young_old_raw_counts_gene5500_mito10perc) (trim most extreme cells to reduce impact of outliers)
# max.adj.var = 10 (maximum value allowed for the estimated adjusted variance, )
# Save plot

pdf("variance_adjusted_variance_knn.error.models.pdf")
varinfo_gene5500_mito10perc <- pagoda.varnorm(knn,
                                counts = young_old_raw_counts_gene5500_mito10perc,
                                trim = 3,
                                max.adj.var = 10,
                                n.cores = 10,
                                plot = TRUE)
dev.off()

# list top overdispersed genes
sort(varinfo_gene5500_mito10perc$arv, decreasing = TRUE)[1:10]

## Controlling for sequencing depth
# Control for the gene coverage (estimated as a number of genes
# with non-zero magnitude per cell) and normalize out that aspect of cell heterogeneity

# Need the expression count with corresponding names for the varinfo
cd <- young_old_raw_counts_gene5500_mito10perc
varinfo_gene5500_mito10perc <- pagoda.subtract.aspect(varinfo_gene5500_mito10perc,
                                                      colSums(cd[, rownames(knn)]>0))

# Save varinfo
save(varinfo_gene5500_mito10perc, file = "varinfo.RData")
