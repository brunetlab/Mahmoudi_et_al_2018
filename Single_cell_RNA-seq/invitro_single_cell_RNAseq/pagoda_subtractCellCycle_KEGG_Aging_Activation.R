## As I learned from Peter Kharchenko and Jean Fan:
## The hierarchical clustering of the cells is based on the genes with are within top 2/3 of the loading magnitude for the PCs of the significantly overdispersed gene sets ($gw slot) --> This is done in pagoda.top.aspects()
## pagoda.cluster.cells() then adds another requirement that the genes should be overdispersed
## Peter usually uses include.aspects = TRUE for pagoda.cluster.cells()
## That all explains why the cell clustering is based on 600 genes for the GO terms but only on 200 genes for the KEGG pathways
## So I will now include all KEGG and all GO gene sets plus the aging signature from Salah from the bulk data.


library(Biobase)
library(scde)
setwd("~/Salah_Katja/SingleCells/Pagoda")
load("../Data/eset.RData")

## Created in the file knn.error.models.R:
## load knn:
load("Data/knn.error.models/knn.error.models.RData")
## Created in the file pagoda_KEGG_Aging.R:
## load varinfo.new:
load("Data/varnorm/subtract.aspect_varnorm_knn.error.models.RData")

library(GSA)
filename <- "../Data/KEGG_2016_Aging_Activation.gmt"
gs <- GSA.read.gmt(filename)
## number of gene sets
n <- length(gs$geneset.names)
## create environment
env <- new.env(parent=globalenv())
invisible(lapply(1:n,function(i) {
  genes <- as.character(unlist(gs$genesets[i]))
  name <- as.character(gs$geneset.names[i])
  assign(name, genes, envir = env)
}))
go.env <- env
class(go.env)

library(GSA)
filename <- "Data/de.novo.gene.sets_gene.clusters_subtract.aspect_varnorm_knn.error.models.gmt"
gs <- GSA.read.gmt(filename)
## number of gene sets
n <- length(gs$geneset.names)
## create environment
env <- new.env(parent=globalenv())
invisible(lapply(1:n,function(i) {
  genes <- as.character(unlist(gs$genesets[i]))
  name <- as.character(gs$geneset.names[i])
  assign(name, genes, envir = env)
}))
go.env.de.novo <- env
class(go.env.de.novo)


## After talking to Salah I will account for these:
## KEGG_cell_cycle, geneCluster.37, geneCluster.119

## First: Correct for KEGG_cell_cycle:

## Get batch signature:
batch.pattern <- pagoda.show.pathways("Cell cycle_Homo sapiens_hsa04110",
                                      varinfo.new,
                                      go.env,
                                      plot = FALSE)
## Subtract the pattern:
varinfo.corrected.pre <- pagoda.subtract.aspect(varinfo.new,
                                                batch.pattern)


## Second: Correct for geneCluster.37 and geneCluster.119:

## Get batch signature:
batch.pattern.new <- pagoda.show.pathways(c("geneCluster.37", "geneCluster.119"),
                                          varinfo.corrected.pre,
                                          go.env.de.novo,
                                          plot = FALSE)
## Subtract the pattern:
varinfo.corrected <- pagoda.subtract.aspect(varinfo.corrected.pre,
                                            batch.pattern.new)

save(varinfo.corrected,
     file = "Data/varnorm/subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")

sort(varinfo.corrected$arv, decreasing = TRUE)[1:10]
##   ACTB   ACTG2 ANGPTL7  ANKRD1    CCL2    CCL7  COL2A1    GPX3  IGFBP5   LECT1
##      10      10      10      10      10      10      10      10      10      10

## Takes long:
## Calculate weighted first PC magnitudes for each gene set in the provided environment:
pwpca <- pagoda.pathway.wPCA(varinfo.corrected,
                             go.env,
                             n.components = 1,
                             n.cores = 1)
## Warning message:
## In bwpca(mat[, lab, drop = FALSE], matw[, lab, drop = FALSE], npcs = n.components,  :
##   When called from R, the RNG seed has to be set at the R level via set.seed()
## dir.create("Data/pathway.wPCA")
save(pwpca,
     file = "Data/pathway.wPCA/pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_varnorm_knn.error.models.RData")


## Evaluate statistical significance of the observed overdispersion for each gene set (KEGG and Aging)
## dir.create("Data/top.aspects")
pdf("Data/top.aspects/top.aspects_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_varnorm_knn.error.models.pdf")
df <- pagoda.top.aspects(pwpca,
                         return.table = TRUE,
                         plot = TRUE,
                         z.score = 1.96)
dev.off()
## Each point on the plot shows the PC1 variance (lambda1) magnitude (normalized by set size) as a function of set size. The red lines show expected (solid) and 95% upper bound (dashed) magnitudes based on the Tracey-Widom model.
save(df,
     file = "Data/top.aspects/top.aspects_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_varnorm_knn.error.models.RData")
head(df)
