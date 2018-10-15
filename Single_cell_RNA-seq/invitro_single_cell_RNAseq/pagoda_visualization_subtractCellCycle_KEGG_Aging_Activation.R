library(Biobase)
library(scde)
library(RColorBrewer)
setwd("~/Salah_Katja/SingleCells/Pagoda")
load("../Data/eset.RData")

## load pwpca:
load("Data/pathway.wPCA/pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_varnorm_knn.error.models.RData")
## load clpca:
load("Data/gene.clusters/gene.clusters_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")
## load varinfo.corrected:
load("Data/varnorm/subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")

## KEGG Pathways and Aging:
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

## De novo gene sets (clpca):
load("Data/gene.clusters/gene.clusters_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")
## read in Broad gmt format
library(GSA)
filename <- "Data/de.novo.gene.sets_gene.clusters_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.gmt"
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



### Visualize significant aspects of heterogeneity:

dir.create("Results/subtractedCellCycle_KEGG_Aging_Activation")

## Evaulate statistical significance of the observed overdispersion for each gene set (also de novo gene sets):
pdf("Results/subtractedCellCycle_KEGG_Aging/top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.pdf")
tam <- pagoda.top.aspects(pwpca,
                          clpca,
                          return.table = TRUE,
                          return.genes = TRUE, ## whether set of genes driving significant aspects should be returned
                          plot = TRUE,
                          z.score = 1.96)
dev.off()

head(tam)

## Determine overall cell clustering (hclust) based on weighted correlation of genes underlying the top aspects transcriptional heterogeneity.
## For some reason, return.table and return.genes have to be FALSE in order to get the list with all results for clustering:
tam <- pagoda.top.aspects(pwpca,
                          clpca,
                          return.table = FALSE,
                          return.genes = FALSE,
                          plot = FALSE,
                          z.score = 1.96)
save(tam,
     file = "Data/top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")
## Look at significant aspects:
sig.aspects <- rownames(tam$xv)
length(grep("geneCluster.", sig.aspects))
## 8
length(grep("_hsa", sig.aspects))
## 74
## Plus: Fibroblast_activation_genes, AGING_SIGNATURE_FIBROBLASTS

hc <- pagoda.cluster.cells(tam,
                           varinfo = varinfo.corrected,
                           verbose = 1,
                           include.aspects = TRUE,
                           return.details = TRUE) ## Whether to return also a the distance matrix and gene values
## clustering cells based on 248 genes and 84 aspect patterns
save(hc,
     file = "Data/cluster.cells_return.details_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")

## I need to run it again with return.details = FALSE in order to pass it to the plot functions:
rm(hc)
hc <- pagoda.cluster.cells(tam,
                           include.aspects = TRUE,
                           varinfo = varinfo.corrected)
save(hc,
     file = "Data/cluster.cells_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")

identical(colnames(tam$xv), colnames(eset))
## TRUE

n.clusters <- 2
col.pal <- brewer.pal(n.clusters, "Paired")
col.cols <- rbind(#groups = col.pal[cutree(hc, n.clusters)],
                  pData(eset)[, "col.mouse.ID"],
                  pData(eset)[, "col.age"])

dir.create("Results/subtractedCellCycle_KEGG_Aging_Activation")

## Reduce redundant aspects part 1:
## Combine pathways/aspects that are driven by the same sets of genes.
## Examines PC loading vectors underlying the identified aspects and clusters based on a product of loading and score correlation.
## Clusters of aspects driven by the same genes are determined based on the distance.threshold and collapsed.
pdf("Results/subtractedCellCycle_KEGG_Aging_Activation/reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.pdf",
    width = 10, height = 13)
tamr <- pagoda.reduce.loading.redundancy(tam,
                                         pwpca,
                                         clpca,
                                         cell.clustering = hc,
                                         distance.threshold = 0.01, ## similarity threshold for grouping interdependent aspects, default: 0.01
                                         abs = TRUE, ## Whether to use absolute correlation.
                                         plot = TRUE,
                                         col.cols = col.cols,
                                         cols = colorRampPalette(c("cornflowerblue", "white", "palevioletred4"), ## low neutral high
                                             space = "Lab")(1024),
                                         margins = c(7, 25))
dev.off()

save(tamr,
     file = "Data/reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")



## Reduce redundant aspects part 2:
## Combine aspects that show similar patterns (i.e. separate the same sets of cells).
## Examines PC loading vectors underlying the identified aspects and clusters based on score correlation.
## Clusters of aspects driven by the same patterns are determined based on the distance.treshold.
## green-to-orange color scheme shows low-to-high weighted PCA scores (aspect patterns), where generally orange indicates higher expression
pdf("Results/subtractedCellCycle_KEGG_Aging/reduce.redundancy_cluster.cells_reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.pdf",
    width = 10, height = 7)
tamr2 <- pagoda.reduce.redundancy(tamr,
                                  plot = TRUE,
                                  distance.threshold = 0.1, ## similarity threshold for grouping interdependent aspects, default: 0.2
                                  cell.clustering = hc,
                                  col.cols = col.cols,
                                  margins = c(7, 25))
dev.off()

save(tamr2,
     file = "Data/reduce.redundancy_cluster.cells_reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.RData")





#### View top aspects:

## (While each row here represents a cluster of pathways, the row names are assigned to be the top overdispersed aspect in each cluster)
pdf("Results/subtractedCellCycle_KEGG_Aging/view.aspects_top_reduce.redundancy_cluster.cells_reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_subtract.aspect_cellCycle_subtract.aspect_varnorm_knn.error.models.pdf",
    width = 10, height = 6)
pagoda.view.aspects(tamr2,
                    cell.clustering = hc,
                    top = 30, ## default: Inf
                    margins = c(7, 25),
                    col.cols = col.cols)
dev.off()



## Extract the names of the 20 aspects (from function 'pagoda.view.aspects'):
top <- 20
top.aspects <- tamr2
rcmvar <- apply(top.aspects$xv, 1, var)
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]
top.aspects$xv <- top.aspects$xv[vi, ]
top.aspects$xvw <- top.aspects$xvw[vi, ]
top.aspects$cnam <- top.aspects$cnam[vi]
top.aspects <- gsub("#PC1# ", "", as.character(rownames(top.aspects$xv)))


## Show genes in top pathways:

dir.create("Results/subtractedCellCycle_KEGG_Aging/top20")


for(i in 1:length(top.aspects)){
    aspect.i <- top.aspects[i]
    file.name <- paste0("Results/subtractedCellCycle_KEGG_Aging/top20/",
                        aspect.i,
                        ".pdf")
    pdf(file.name, width = 5, height = 10)
    if(grepl("geneCluster.", aspect.i)){
        pagoda.show.pathways(pathways = aspect.i,
                             varinfo = varinfo.corrected,
                             goenv = go.env.de.novo,
                             cell.clustering = hc,
                             n.genes = 30,
                             colcols = col.cols,
                             margins = c(1,5),
                             show.cell.dendrogram = TRUE,
                             showRowLabels = TRUE,
                             showPC = TRUE)
    } else {
        pagoda.show.pathways(pathways = aspect.i,
                             varinfo = varinfo.corrected,
                             goenv = go.env,
                             cell.clustering = hc,
                             n.genes = 30,
                             colcols = col.cols,
                             margins = c(1,5),
                             show.cell.dendrogram = TRUE,
                             showRowLabels = TRUE,
                             showPC = TRUE)
    }
    dev.off()
}

## TNFSF11, MET, TGFB2, CCL17, VEGFA, CXCL12, TGFB3, TSLP, CSF1, CCL2, CCL7, CXCL1, HGF, CCL8, IL6
