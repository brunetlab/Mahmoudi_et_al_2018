# Objective: Generate a PCA plot of PC1 and PC3 of metabolome good and bad old Fibroblasts at pasage 3
# ED_Fig3j
# Use processed data from Pre-processing_Metabolomics_Fibroblasts.R
#####################################################################################################################
# Set working directory
setwd("Metabolome/")

# Load data
load("/Metabolome/RAW_data/Metabolome_range-scaled_values_pos_neg.Rda")
colnames(data_nz.norm)

# Extract only old samples 
data_nz.norm.old <- data_nz.norm[,grep("29m", colnames(data_nz.norm))]
colnames(data_nz.norm.old)

# PCA using prcomp function
e.pca <- prcomp(t(data_nz.norm.old), scale = T)

# get the precentage for each component
summary(e.pca)
# Importance of components:
#                           PC1     PC2     PC3     PC4     PC5      PC6      PC7       PC8
# Standard deviation     43.8354 34.4230 27.2296 24.8482 23.6351 22.91856 21.53492 6.194e-14
# Proportion of Variance  0.3196  0.1971  0.1233  0.1027  0.0929  0.08735  0.07713 0.000e+00
# Cumulative Proportion   0.3196  0.5166  0.6399  0.7426  0.8355  0.92287  1.00000 1.000e+00

# color code
## good old (culture 29m1|29m2|29m6|29m7) = cyan blue
## bad old (culture 29m3|29m4|29m5|29m8) =  dark blue

# Make PCA plot based PC2 and PC3
pdf("ED_Fig3j.pdf", width =4.3, height =5 )
  par(mar=c(4.1,4.1,6,1), xpd=TRUE)
  plot(e.pca$x[,c(2,3)],
       pch = 15,  
       xlab = 'PC1 (20%)', 
       ylab = 'PC2 (12%)',
       cex.lab = 1.5, 
       col = "white",
       cex=0.6, 
       ylim = c(-45, 45),
       xlim=c(-70,50)
       )
  points(e.pca$x[grep("29m1|29m2|29m6|29m7", colnames(data_nz.norm)) ,c(2,3)],pch=15, cex=3, 
         col=rgb(red=31,green=217,blue=255, alpha = 204, maxColorValue = 255))
  points(e.pca$x[grep("29m3|29m4|29m5|29m8", colnames(data_nz.norm)) ,c(2,3)],pch=15, cex=3, 
         col=rgb(red=23,green=0,blue=134, alpha = 204, maxColorValue = 255))
  legend("topright",c("Good Old","Bad Old"), inset=c(0,-0.33),horiz=T,
         col=c("lightblue","darkblue"), text.col=c("lightblue","darkblue"), pch=c(17,17), cex = 1.5 , pt.cex=3)
dev.off()

