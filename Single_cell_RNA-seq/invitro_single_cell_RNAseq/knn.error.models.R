## scde works with module load destiny.
## Pagoda: Pathway and Gene Set Overdispersion Analysis
## I follow the tutorial of pagoda on: http://hms-dbmi.github.io/scde/pagoda.html

n.cores <- 1

library(Biobase)
library(scde)
setwd("/Volumes/Seagate\ Backup\ Plus\ Drive/Salah/Single_Cell_data/PAGODA")
load("./Data/eset.RData")

## Build error models (knn.error.models saves plots into current working dir):
## with groups as cell types:
dir.create("Data/knn.error.models")
setwd("Data/knn.error.models")
knn <- knn.error.models(exprs(eset),
                        verbose = 1,
                        n.cores = n.cores)
save(knn,
     file = "knn.error.models.RData")
