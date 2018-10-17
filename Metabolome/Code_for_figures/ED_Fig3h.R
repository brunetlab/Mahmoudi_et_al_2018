# Objective: Generate a clustering tree of metabolome data from young and old cells, using correlation and average
# Extended Data Figure3h
# Use processed data from Pre-processing_Metabolomics_Fibroblasts.R

###########################################################################################
# Required R packages
library('pvclust')    # v2.0-0
sessionInfo()
###########################################################################################
# Set working directory
setwd("Metabolome/")

# Upload data
load("/Metabolome/RAW_data/Metabolome_range-scaled_values_pos_neg.Rda")

# Clustering using pvclust, use correlation and average
e.clustering <- pvclust(data_nz.norm,nboot=1000,method.dist="correlation",method.hclust="average")

# color code
## young = orange
## good old (culture 29m1|29m2|29m6|29m7) = cyan blue
## bad old (culture 29m3|29m4|29m5|29m8) =  dark blue


pdf("ED_Fig_3h.pdf", width =4.3, height =5 )
plot(e.clustering,col.pv=c("white","white","white"))
dev.off()





