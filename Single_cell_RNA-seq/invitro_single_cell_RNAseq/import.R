## Make expressionSet with count matrix.

setwd("~/Salah_Katja/SingleCells")
source("~/My_R_functions/make.transparent.R")
library(Biobase)
library(RColorBrewer)

#### Counts:

load("Data/DESeq2_Single_cell_reads.RData")

## Keep only genes that Salah used for his PCAs:
load("Data/DESeq2_DEG_all_genes_vst_values.RData")
genes <- rownames(vstMat)
counts <- as.matrix(smCountTable.1[genes,])

eset <- ExpressionSet(assayData = counts)

sum(exprs(eset) == 0) / ( nrow(eset) * ncol(eset) )
## 0.369

pData(eset)$age <- "young"
pData(eset)$age[grepl("29M", colnames(eset))] <- "old"
pData(eset)$col.age <- make.transparent("orange", alpha = 200)
pData(eset)$col.age[pData(eset)$age == "old"] <- make.transparent("darkblue", alpha = 200)
pData(eset)$mouse.ID <- NA
pData(eset)$mouse.ID[pData(eset)$age == "young"] <- as.integer(substr(colnames(eset), 7, 7)[pData(eset)$age == "young"])
pData(eset)$mouse.ID[pData(eset)$age == "old"] <- as.integer(substr(colnames(eset), 8, 8)[pData(eset)$age == "old"])
## Reprogramming efficiency:
pData(eset)$RE <- numeric(length = ncol(eset))
pData(eset)$RE[grepl("3M3c4", colnames(eset))] <- 1.03
pData(eset)$RE[grepl("3M6c4", colnames(eset))] <- 0.78
pData(eset)$RE[grepl("3M7c4", colnames(eset))] <- 1.11
pData(eset)$RE[grepl("29M2c4", colnames(eset))] <- 1.44
pData(eset)$RE[grepl("29M4c4", colnames(eset))] <- 0.43
pData(eset)$RE[grepl("29M8c4", colnames(eset))] <- 0.54

col.pal <- brewer.pal(8, "Accent")
pData(eset)$col.mouse.ID <- col.pal[pData(eset)$mouse.ID]

rownames(eset) <- toupper(rownames(eset))

save(eset,
     file = "Data/eset.RData")


#### VST normalized data from Salah:

load("Data/DESeq2_DEG_all_genes_vst_values.RData")

eset <- ExpressionSet(assayData = vstMat)

pData(eset)$age <- "young"
pData(eset)$age[grepl("29M", colnames(eset))] <- "old"
pData(eset)$col.age <- make.transparent("orange", alpha = 200)
pData(eset)$col.age[pData(eset)$age == "old"] <- make.transparent("darkblue", alpha = 200)
pData(eset)$mouse.ID <- NA
pData(eset)$mouse.ID[pData(eset)$age == "young"] <- as.integer(substr(colnames(eset), 7, 7)[pData(eset)$age == "young"])
pData(eset)$mouse.ID[pData(eset)$age == "old"] <- as.integer(substr(colnames(eset), 8, 8)[pData(eset)$age == "old"])
## Reprogramming efficiency:
pData(eset)$RE <- numeric(length = ncol(eset))
pData(eset)$RE[grepl("3M3c4", colnames(eset))] <- 1.03
pData(eset)$RE[grepl("3M6c4", colnames(eset))] <- 0.78
pData(eset)$RE[grepl("3M7c4", colnames(eset))] <- 1.11
pData(eset)$RE[grepl("29M2c4", colnames(eset))] <- 1.44
pData(eset)$RE[grepl("29M4c4", colnames(eset))] <- 0.43
pData(eset)$RE[grepl("29M8c4", colnames(eset))] <- 0.54

col.pal <- brewer.pal(8, "Accent")
pData(eset)$col.mouse.ID <- col.pal[pData(eset)$mouse.ID]

rownames(eset) <- toupper(rownames(eset))

save(eset,
     file = "Data/eset_vst.RData")
