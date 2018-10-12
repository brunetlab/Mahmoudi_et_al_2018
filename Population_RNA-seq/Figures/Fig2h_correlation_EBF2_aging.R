# Objective: To generate a correlation plot with EBF2 expressin and age in human primary fibroblasts

# Data: Vst-values from DESeq2 analysis of healthy human fibroblastss
###########################################################################################
# Need the following packages
library("DESeq2")   # v1.20.0
sessionInfo()
######################################################################################################################################
###############################    removing low expressed genes  ----------------------------------------------------------
######################################################################################################################################
# Set woking directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Correlation_plot/')

# load DESeq2 analysis
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Human_healthy/DESeq2_vst_human_Fibs_healthy_all_genes_linear_model.Rda")

# specifiy batch and condtion, (for some reason had to add 0.001 to make them numeric)
Age = c(69, 1, 0, 30,	
        0, 43, 29,	
        87, 71, 89,
        1, 66, 29)

Ebf2 <- vstMat[grep("EBF2", row.names(vstMat)),]

pdf("Correlation_Age_Ebf2_expression.pdf", width =5, height =5, onefile=F )
plot(Age,Ebf2[1,], col = "grey", pch =15)
abline(lm(Ebf2[1,]~ Age), col= "red")
dev.off()
