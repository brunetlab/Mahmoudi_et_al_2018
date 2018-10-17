# Objective: Generate a PCA plot of PC1-PC2 PC1-PC3 of metabolome of young and old iPSC at passage 23
# ED_Fig4e
# Use processed data from Pre-processing_Metabolomics_iPSC.R
####################################################################################
# Required R packages

####################################################################################
# Set working directory
setwd("Metabolome/")

# Upload data
load("/Metabolome/RAW_data/Metabolome_iPSC_range-scaled_values_pos_neg.Rda")

# PCA using prcomp function
e.pca <- prcomp(t(data_nz.norm), scale = T)

# get the precentage for each component
summary(e.pca)
# Importance of components:
#                           PC1     PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10     PC11    PC12     PC13
# Standard deviation     39.0629 29.4660 22.67017 22.00682 20.76572 19.21205 18.80666 18.13016 17.61145 17.40286 16.69586 16.4658 6.09e-14
# Proportion of Variance  0.2527  0.1438  0.08512  0.08021  0.07142  0.06113  0.05858  0.05444  0.05137  0.05016  0.04617  0.0449 0.00e+00
# Cumulative Proportion   0.2527  0.3965  0.48163  0.56184  0.63326  0.69439  0.75296  0.80740  0.85877  0.90893  0.95510  1.0000 1.00e+00

# find the color code for orange and darkblue
# color code
## young = orange
## good old 1 (culture 29m6) = cyan blue
## good old 2 (culture 29m7) = light cyan blue
## bad old 1 (culture 29m2) =  dark blue


# Make PCA plot based PC1 and PC2
pdf("PCA_Metabolome_range-scaled_pos_neg_iPSC_PC1-2.pdf", width =4.3, height =5 )
  par(mar=c(4.1,4.1,6,1), xpd=TRUE)
  plot(e.pca$x[,c(1,2)],
       pch = 15,  
       xlab = 'PC1 (25%)', 
       ylab = 'PC2 (14%)',
       cex.lab = 1.5, 
       col = "white",
       cex=0.6 
       # ylim = c(-27, 52), 
       # xlim=c(-110,48)
       ) 
  points(e.pca$x[grep("3m", colnames(data_nz.norm)) ,c(1,2)],  pch=15, cex=3, 
         col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255))
  points(e.pca$x[grep("29m6", colnames(data_nz.norm)) ,c(1,2)],pch=15, cex=3, 
         col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255))
  points(e.pca$x[grep("29m7", colnames(data_nz.norm)) ,c(1,2)],pch=15, cex=3, 
         col=rgb(red=31,green=217,blue=255, alpha = 120, maxColorValue = 255))
  points(e.pca$x[grep("29m2", colnames(data_nz.norm)) ,c(1,2)],pch=15, cex=3, 
         col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255))
  legend("topright",c("Young","Good Old","Bad Old"), inset=c(0,0),horiz=T,
         col=c(col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255),
               col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255),
               col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255)),
         pch=c(15,15,15), cex = 1 , pt.cex=2)
dev.off()


pdf("PCA_Metabolome_range-scaled_pos_neg_iPSC_PC1-3.pdf", width =4.3, height =5 )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
plot(e.pca$x[,c(1,3)],
     pch = 15,  
     xlab = 'PC1 (25%)', 
     ylab = 'PC3 (9%)',
     cex.lab = 1.5, 
     col = "white",
     cex=0.6 
     # ylim = c(-27, 52), 
     # xlim=c(-110,48)
) 
points(e.pca$x[grep("3m", colnames(data_nz.norm)) ,c(1,3)],  pch=15, cex=3, 
       col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255))
points(e.pca$x[grep("29m6", colnames(data_nz.norm)) ,c(1,3)],pch=15, cex=3, 
       col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255))
points(e.pca$x[grep("29m7", colnames(data_nz.norm)) ,c(1,3)],pch=15, cex=3, 
       col=rgb(red=31,green=217,blue=255, alpha = 120, maxColorValue = 255))
points(e.pca$x[grep("29m2", colnames(data_nz.norm)) ,c(1,3)],pch=15, cex=3, 
       col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255))
legend("topright",c("Young","Good Old","Bad Old"), inset=c(0,0),horiz=T,
       col=c(col=rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255),
             col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255),
             col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255)),
       pch=c(15,15,15), cex = 1 , pt.cex=2)
dev.off()


