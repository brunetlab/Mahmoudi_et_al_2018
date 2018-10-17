# Objective: Generate a PCA plot of PC1 and PC2 of metabolome data of young and old Fibroblasts
# Extended Data Figure 3c-d
# Use processed data from Pre-processing_Metabolomics_Fibroblasts.R

####################################################################################
# Required R packages

####################################################################################
# Set working directory
setwd("Metabolome/")

# Upload data
load("/Metabolome/RAW_data/Metabolome_range-scaled_values_pos_neg.Rda")

# PCA using prcomp function
e.pca <- prcomp(t(data_nz.norm), scale = T)

# get the precentage for each component
summary(e.pca)
# Importance of components:
#                           PC1     PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10    PC11    PC12     PC13     PC14     PC15      PC16
# Standard deviation     40.1108 26.9578 23.61767 20.36105 19.52266 17.25011 16.45501 16.22431 15.53229 15.29049 14.6511 14.6097 14.12149 14.09776 13.82115 4.322e-14
# Proportion of Variance  0.2676  0.1209  0.09276  0.06895  0.06339  0.04949  0.04503  0.04378  0.04012  0.03888  0.0357  0.0355  0.03316  0.03305  0.03177 0.000e+00
# Cumulative Proportion   0.2676  0.3884  0.48119  0.55014  0.61352  0.66301  0.70804  0.75182  0.79194  0.83082  0.8665  0.9020  0.93518  0.96823  1.00000 1.000e+00

# color code
## young = orange
## good old (culture 29m1|29m2|29m6|29m7) = cyan blue
## bad old (culture 29m3|29m4|29m5|29m8) =  dark blue

# Make PCA plot based PC1 and PC2
pdf("Extended_Data_Fig2d_PCA_PC1-2", width =4.3, height =5 )
  par(mar=c(4.1,4.1,6,1), xpd=TRUE)
  plot(e.pca$x[,c(1,2)],pch = 15,  xlab = 'PC1 (27%)', ylab = 'PC3 (12%)',cex.lab = 1.5, col = "white",cex=0.6, ylim = c(-80, 80), xlim=c(-100,57)) 
  points(e.pca$x[grep("3m", colnames(data_nz.norm)) ,c(1,2)],  pch=15, cex=3, 
         col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255))
  points(e.pca$x[grep("29m1|29m2|29m6|29m7", colnames(data_nz.norm)) ,c(1,2)],pch=15, cex=3, 
         col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255))
  points(e.pca$x[grep("29m3|29m4|29m5|29m8", colnames(data_nz.norm)) ,c(1,2)],pch=15, cex=3, 
         col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255))
  legend("topright",c("Young","Good Old","Bad Old"), inset=c(0,0),horiz=T,
   col=c(col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255),
         col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255),
         col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255)),
   pch=c(15,15,15), cex = 1 , pt.cex=2)
dev.off()

# Make PCA plot based PC1 and PC3
pdf("Extended_Data_Fig2d_PCA_PC1-2", width =4.3, height =5 )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
plot(e.pca$x[,c(1,3)],pch = 15,  xlab = 'PC1 (27%)', ylab = 'PC3 (9%)',cex.lab = 1.5, col = "white",cex=0.6, ylim = c(-80, 80), xlim=c(-100,57)) 
points(e.pca$x[grep("3m", colnames(data_nz.norm)) ,c(1,3)],  pch=15, cex=3, 
       col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255))
points(e.pca$x[grep("29m1|29m2|29m6|29m7", colnames(data_nz.norm)) ,c(1,3)],pch=15, cex=3, 
       col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255))
points(e.pca$x[grep("29m3|29m4|29m5|29m8", colnames(data_nz.norm)) ,c(1,3)],pch=15, cex=3, 
       col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255))
legend("topright",c("Young","Good Old","Bad Old"), inset=c(0,0),horiz=T,
       col=c(col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255),
             col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255),
             col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255)),
       pch=c(15,15,15), cex = 1 , pt.cex=2)
dev.off()
