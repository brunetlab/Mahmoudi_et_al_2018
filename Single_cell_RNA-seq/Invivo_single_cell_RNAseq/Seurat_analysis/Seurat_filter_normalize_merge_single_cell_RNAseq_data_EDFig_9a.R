# Objective: Process single cell data from young and old wounds at day 7 using Seurat
# Filter based on min and max number of genes expressed, and % mitochondrial gene expression 
# Normalize young and old data datsets, and merge the into one

# Generate a Seurat object that contain the merged normalized data. Will be used for downstream analysis such as tSNE
####################################################################################
# Required packages
library(Seurat)         # Seurat_2.3.4
library(ggplot2)        # ggplot2_3.0.0
library(dplyr)          # dplyr_0.7.6
sessionInfo()
####################################################################################
####################################################################################
setwd("~/Desktop/Dropbox/For_Matt/")

# load("~/Dropbox/RNAseq/DESeq2_analysis_Fib/180430_Single_cell_in_vivo/Seurat/Old_files/old_exprs_matrix.RData")
# load("~/Dropbox/RNAseq/DESeq2_analysis_Fib/180430_Single_cell_in_vivo/Seurat/Old_files/young_exprs_matrix.RData")
load("./Essentials/Expression_data_Seurat/young_exprs_matrix.RData")
load("./Essentials/Expression_data_Seurat/old_exprs_matrix.RData")


#Establishes Seurat matrix with cutoff of 300 genes per cell and >=5 cell expressing a gene
# Set up young fibroblast
young <- CreateSeuratObject(raw.data = young_exprs_matrix_agg_mod,
                           project = "10X_fibroblast_Young_Old_Fibroblast",
                           min.cells = 5)
young@meta.data$old <- "Young Fibroblast"

# Determine percent mitochondrial genes for subsetting
mito.genes.young <- grep(pattern = "mt-", x = rownames(x = young@data), value = TRUE)
percent.mito.young <- Matrix::colSums(young@raw.data[mito.genes.young, ])/Matrix::colSums(young@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
young <- AddMetaData(object = young,
                     metadata = percent.mito.young,
                     col.name = "percent.mito")

# Based on the QC data/plots determine filter data at min 500 and max 5500 genes expressed and mito content 10%
young <- FilterCells(object = young,
                     subset.names = c("nGene", "percent.mito"),
                     low.thresholds = c(500, -Inf),
                     high.thresholds = c(5500, 0.1))
young
# An object of class seurat in project 10X_fibroblast_Young_Old_Fibroblast
# 13998 genes across 1592 samples.

# Normalize the data
young <- NormalizeData(object = young,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)
young <- ScaleData(young,
                   vars.to.regress = c("nUMI", "percent.mito"),
                   display.progress = F)

# Set up old fibroblast
old <- CreateSeuratObject(raw.data = old_exprs_matrix_agg_mod,
                          project = "10X_fibroblast_Young_Old_Fibroblast",
                          min.cells = 5)
old@meta.data$old <- "Old Fibroblast"

# Determine percent mitochondrial genes for subsetting
mito.genes.old <- grep(pattern = "mt-", x = rownames(x = old@data), value = TRUE)
percent.mito.old <- Matrix::colSums(old@raw.data[mito.genes.old, ])/Matrix::colSums(old@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
old <- AddMetaData(object = old,
                   metadata = percent.mito.old,
                   col.name = "percent.mito")

# Based on the QC data/plots determine filter data at min 500 and max 5500 genes expressed and mito content 10%
old <- FilterCells(object = old,
                     subset.names = c("nGene", "percent.mito"),
                   low.thresholds = c(500, -Inf),
                   high.thresholds = c(5500, 0.1))
old
# An object of class seurat in project 10X_fibroblast_Young_Old_Fibroblast
# 14067 genes across 1444 samples.

old <- NormalizeData(old)
old <- ScaleData(old,
                 vars.to.regress = c("nUMI", "percent.mito"),
                 display.progress = F)

# Gene selection for input to CCA
young <- FindVariableGenes(young, do.plot = T)
old <- FindVariableGenes(old, do.plot = T)
g.1 <- head(rownames(young@hvg.info), 1000)
g.2 <- head(rownames(old@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(young@scale.data))
genes.use <- intersect(genes.use, rownames(old@scale.data))

# Perform a canonical correlation analysis (CCA)
fibroblast.combined <- RunCCA(young,
                              old,
                              genes.use = genes.use,
                              num.cc = 50)

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = fibroblast.combined,
              reduction.use = "cca",
              group.by = "old",
              pt.size = 0.5,
              fibroblast.combined = TRUE)
p2 <- VlnPlot(object = fibroblast.combined,
              features.plot = "CC2",
              group.by = "old",
              do.return = TRUE)
plot_grid(p2)

p3 <- MetageneBicorPlot(fibroblast.combined,
                        grouping.var = "old",
                        dims.eval = 1:50,
                        display.progress = FALSE)

PrintDim(object = fibroblast.combined,
         reduction.type = "cca",
         dims.print = 1:2,
         genes.print = 10)

DimHeatmap(object = fibroblast.combined,
           reduction.type = "cca",
           cells.use = 500,
           dim.use = 1:9,
           do.balanced = TRUE)

fibroblast.combined <- AlignSubspace(fibroblast.combined,
                                     reduction.type = "cca",
                                     grouping.var = "old",
                                     dims.align = 1:35)

saveRDS(fibroblast.combined, file = "MATT_fibroblast.combined.rds")


