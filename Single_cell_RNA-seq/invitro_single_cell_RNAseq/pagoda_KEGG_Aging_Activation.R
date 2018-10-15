library(scde)
library(Biobase)
setwd("~/Salah_Katja/SingleCells/Pagoda")
load("../Data/eset.RData")


## Created in the file knn.error.models.R:
## load knn:
load("Data/knn.error.models/knn.error.models.RData")
## Created in the file pagoda_KEGG_Aging.R:
## load varinfo.new:
load("Data/varnorm/subtract.aspect_varnorm_knn.error.models.RData")

## List of overdispersed genes:
sort(varinfo.new$arv, decreasing = TRUE)[1:10]
##    ACTB   ACTG2 ANGPTL7  ANKRD1    CCL2    CCL7  COL2A1    GPX3  IGFBP5   LECT1
##      10      10      10      10      10      10      10      10      10      10


## KEGG and Aging gene sets:
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


## Takes long:
## Calculate weighted first PC magnitudes for each gene set in the provided environment:
pwpca <- pagoda.pathway.wPCA(varinfo.new,
                             go.env,
                             n.components = 1,
                             n.cores = 1)
## Warning message:
## In bwpca(mat[, lab, drop = FALSE], matw[, lab, drop = FALSE], npcs = n.components,  :
##   When called from R, the RNG seed has to be set at the R level via set.seed()
## dir.create("Data/pathway.wPCA")
save(pwpca,
     file = "Data/pathway.wPCA/pathway.wPCA_KEGG_Aging_Activation_varnorm_knn.error.models.RData")



## Evaluate statistical significance of the observed overdispersion for each gene set
## dir.create("Data/top.aspects")
pdf("Data/top.aspects/top.aspects_pathway.wPCA_KEGG_Aging_Activation_varnorm_knn.error.models.pdf")
df <- pagoda.top.aspects(pwpca,
                         return.table = TRUE,
                         plot = TRUE,
                         z.score = 1.96)
dev.off()
## Each point on the plot shows the PC1 variance (lambda1) magnitude (normalized by set size) as a function of set size. The red lines show expected (solid) and 95% upper bound (dashed) magnitudes based on the Tracey-Widom model.
save(df,
     file = "Data/top.aspects/top.aspects_pathway.wPCA_KEGG_Aging_Activation_varnorm_knn.error.models.RData")
head(df)


## Load de novo lists (clpca):
load("Data/gene.clusters/gene.clusters_subtract.aspect_varnorm_knn.error.models.RData")



## Write de novo gene sets into .gmt file to then create the environment together with the KEGG pathways:
create.gmt <- function(pathway.list, outfile){
  max.l <- max(sapply(pathway.list, length))
  l <- lapply(pathway.list, function(x){
    length(x) <- max.l
    return(x)
  })
  tab <- do.call("rbind", l)
  tab[is.na(tab)] <- ""
  tab <- cbind(rep("", nrow(tab)), tab)
  write.table(tab, outfile,
              sep = "\t", quote = FALSE, col.names = FALSE)
}
de.novo.pathways <- clpca$clusters

create.gmt(pathway.list = de.novo.pathways,
           outfile = "Data/de.novo.gene.sets_gene.clusters_subtract.aspect_varnorm_knn.error.models.gmt")
