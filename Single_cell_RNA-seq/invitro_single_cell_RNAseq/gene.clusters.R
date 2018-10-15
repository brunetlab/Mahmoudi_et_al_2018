n.cores <- 20

library(Biobase)
library(scde)
setwd("~/Salah_Katja/SingleCells/Pagoda")

### Careful: Takes 2.5 hours when n.cores = 1!
### Evaluate overdispersion of de novo gene sets:

load("Data/varnorm/subtract.aspect_varnorm_knn.error.models.RData")
dir.create("Data/gene.clusters")

pdf("Data/gene.clusters/gene.clusters_subtract.aspect_varnorm_knn.error.models.pdf",
    width = 10, height = 5)
clpca <- pagoda.gene.clusters(varinfo.new,
                              n.cores = n.cores,
                              plot = TRUE)
dev.off()
save(clpca,
     file = "Data/gene.clusters/gene.clusters_subtract.aspect_varnorm_knn.error.models.RData")


