# Objective: Generate a clustering tree of metabolome data from young and old derived iPSC, using correlation and average
# Extended Data Figure 4f
# Use range scaled signal intensity data

###########################################################################################
# Need the following package
library('pvclust')    # v2.0-0
sessionInfo()
###########################################################################################
# Set working directory
setwd("Metabolome/")

# Load data
load("/Metabolome/RAW_data/Metabolome_iPSC_range-scaled_values_pos_neg.Rda")


# Clustering using pvclust, use correlation and average (as average is less sensitive to outliers)
e.clustering <- pvclust(data_nz.norm,nboot=10,method.dist="correlation",method.hclust="average")

# find the color code for orange and darkblue
# color code
## young = orange
## good old 1 (culture 29m6) = cyan blue
## good old 2 (culture 29m7) = light cyan blue
## bad old 1 (culture 29m2) =  dark blue


pdf("Clustering_Metabolome_range-scaled_pos_neg_iPSC.pdf", width =4.3, height =5 )
plot(e.clustering,col.pv=c("white","white","white"))
dev.off()




