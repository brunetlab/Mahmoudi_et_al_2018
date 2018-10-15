library(Biobase)
library(scde)
library(RColorBrewer)
setwd("~/Salah_Katja/SingleCells/Pagoda")
load("../Data/eset.RData")

## load pwpca:
load("Data/pathway.wPCA/pathway.wPCA_KEGG_Aging_Activation_varnorm_knn.error.models.RData")
## load clpca:
load("Data/gene.clusters/gene.clusters_subtract.aspect_varnorm_knn.error.models.RData")
## load varinfo.new:
load("Data/varnorm/subtract.aspect_varnorm_knn.error.models.RData")


## KEGG Pathways:
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



## De novo gene sets:
#load("Data/de.novo.gene.sets_gene.clusters_subtract.aspect_varnorm_knn.error.models.gmt")
## read in Broad gmt format
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



### Visualize significant aspects of heterogeneity:


## Evaulate statistical significance of the observed overdispersion for each gene set (also de novo gene sets):
pdf("Results/top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.pdf")
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
     file = "Data/top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.RData")
hc <- pagoda.cluster.cells(tam,
                           varinfo = varinfo.new,
                           include.aspects = TRUE,
                           verbose = 1,
                           return.details = TRUE) ## Whether to return also a the distance matrix and gene values
## clustering cells based on 405 genes and 97 aspect patterns
save(hc,
     file = "Data/cluster.cells_return.details_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.RData")


## I need to run it again with return.details = FALSE in order to pass it to the plot functions:
rm(hc)
hc <- pagoda.cluster.cells(tam,
                           include.aspects = TRUE,
                           varinfo = varinfo.new)
save(hc,
     file = "Data/cluster.cells_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.RData")



n.clusters <- 5
col.pal <- brewer.pal(n.clusters, "Paired")
col.cols <- rbind(#groups = col.pal[cutree(hc, n.clusters)],
                  pData(eset)[colnames(tam$xv), "col.mouse.ID"],
                  pData(eset)[colnames(tam$xv), "col.age"])

dir.create("Results/KEGG_Aging_Activation")

## Reduce redundant aspects part 1:
## Combine pathways/aspects that are driven by the same sets of genes.
## Examines PC loading vectors underlying the identified aspects and clusters based on a product of loading and score correlation.
## Clusters of aspects driven by the same genes are determined based on the distance.threshold and collapsed.
pdf("Results/KEGG_Aging_Activation/reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.pdf",
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
                                         margins = c(7, 30))
dev.off()

save(tamr,
     file = "Data/reduce.loading.redundancy_top.aspects_gene.clusters_pathway.wPCA_KEGG_Aging_Activation_subtract.aspect_varnorm_knn.error.models.RData")



## Extract the names of all aspects (from function 'pagoda.view.aspects'):
top <- Inf
top.aspects <- tam
rcmvar <- apply(top.aspects$xv, 1, var)
vi <- order(rcmvar, decreasing = TRUE)[1:min(length(rcmvar), top)]
top.aspects$xv <- top.aspects$xv[vi, ]
top.aspects$xvw <- top.aspects$xvw[vi, ]
top.aspects$cnam <- top.aspects$cnam[vi]
top.aspects <- gsub("#PC1# ", "", as.character(rownames(top.aspects$xv)))


## Show genes in top pathways:

dir.create("Results/KEGG_Aging_Activation/Heatmaps")

for(i in 1:length(top.aspects)){
  aspect.i <- top.aspects[i]
  aspect.i.file <- gsub("/", ".", aspect.i)
  file.name <- paste0("Results/KEGG_Aging_Activation/Heatmaps/",
                      aspect.i.file,
                      ".pdf")
  pdf(file.name, width = 5, height = 4)
  if(grepl("geneCluster.", aspect.i)){
    pagoda.show.pathways(pathways = aspect.i,
                         varinfo = varinfo.new,
                         goenv = go.env.de.novo,
                         cell.clustering = hc,
                         n.genes = 20,
                         colcols = col.cols,
                         margins = c(1,5),
                         show.cell.dendrogram = TRUE,
                         showRowLabels = TRUE,
                         showPC = TRUE)
  } else {
    pagoda.show.pathways(pathways = aspect.i,
                         varinfo = varinfo.new,
                         goenv = go.env,
                         cell.clustering = hc,
                         n.genes = 20,
                         colcols = col.cols,
                         margins = c(1,5),
                         show.cell.dendrogram = TRUE,
                         showRowLabels = TRUE,
                         showPC = TRUE)
  }
  dev.off()
}

### Correct for: Cell cycle_Homo sapiens_hsa04110, geneCluster.37, gene.cluster.119

## write out the genes in geneCluster.37 and geneCluster.119:

cluster.37 <- get("geneCluster.37", go.env.de.novo)
cluster.119 <- get("geneCluster.119", go.env.de.novo)

write.table(data.frame(geneCluster.37 = cluster.37),
            file = "Data/geneCluster.37.csv",
            row.names = FALSE,
            sep = "\t")
write.table(data.frame(geneCluster.119 = cluster.119),
            file = "Data/geneCluster.119.csv",
            row.names = FALSE,
            sep = "\t")
