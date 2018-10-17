# Objective: To determine significantly different metabolite features between young and old iPSC at passage 23
# Analysis will be used to generate Extended Data Figure 4g
# Use median signal intensities of pos and neg mode
####################################################################################
# Need the following packages
library(qvalue)    # qvalue_2.6.0   
sessionInfo()
####################################################################################
# set working directory
setwd("/Metabolome/")

# In a loop calculate mean 3m, mean 29m, mean across all samples, Log2FoldChange, Log2 fold change (Old Fib vs Young Fib), 
# determine significance between young and old using a non-parametric test (wilcoxon-test), 
# finally calculate FDR (BH) 

# Upload Neg mode data
load("~/Desktop/Github/Metabolome/RAW_data/Metabolome_iPSC_median_signal_intensity_NEG.Rda") 
# Determine the relevant columns
colnames(median_n.all)

median_n.all.stat <- as.data.frame(median_n.all)
for (i in 1:nrow(median_n.all)) {
  x <- as.numeric(median_n.all[i,grep("3m", colnames(median_n.all))])
  y <- as.numeric(median_n.all[i,grep("29m", colnames(median_n.all))])
  median_n.all.stat$mean.3m[i] <- mean(y)
  median_n.all.stat$mean.29m[i] <- mean(x)
  median_n.all.stat$baseMean[i] <- mean(c(x,y))
  median_n.all.stat$Log2FC[i] <- log2((mean(x)/mean(y)))
  wx <- wilcox.test(x,y, alternative="two.sided")
  median_n.all.stat$Wilcoxon[i] <- wx$p.value
}
median_n.all.stat$p.adjusted.BH <- p.adjust(median_n.all.stat$Wilcoxon, "BH")
qval_neg <- qvalue(median_n.all.stat$Wilcoxon)
median_n.all.stat$Q_value <- qval_neg$qvalues
# Determine how many features have BH < 0.05
length(which(median_n.all.stat$p.adjusted.BH < 0.05)) # 0
length(which(median_n.all.stat$Q_value < 0.05)) # 0


# Upload Pos mode data
load("~/Desktop/Github/Metabolome/RAW_data/Metabolome_iPSC_median_signal_intensity_POS.Rda")
colnames(median_p.all)

median_p.all.stat <- as.data.frame(median_p.all)
for (i in 1:nrow(median_p.all)) {
  x <- as.numeric(median_p.all[i,grep("3m", colnames(median_n.all))])
  y <- as.numeric(median_p.all[i,grep("29m", colnames(median_n.all))])
  median_p.all.stat$mean.3m[i] <- mean(y)
  median_p.all.stat$mean.29m[i] <- mean(x)
  median_p.all.stat$baseMean[i] <- mean(c(x,y))
  median_p.all.stat$Log2FC[i] <- log2((mean(x)/mean(y)))
  wx <- wilcox.test(x,y, alternative="two.sided")
  median_p.all.stat$Wilcoxon[i] <- wx$p.value
}
median_p.all.stat$p.adjusted.BH <- p.adjust(median_p.all.stat$Wilcoxon, "BH")
qval_pos <- qvalue(median_p.all.stat$Wilcoxon)
median_p.all.stat$Q_value <- qval_pos$qvalues
# Determine how many features have BH < 0.05
length(which(median_p.all$p.adjusted.BH < 0.05)) # 0
length(which(median_p.all$Q_value < 0.05)) # 0


# Merge the two tables together including all data
data_all <- rbind(median_p.all,median_p.all)
write.csv(data_all, file='Differential_analysis_iPSC_metabolome_pos_neg.csv')
save(data_all, file='Differential_analysis_iPSC_metabolome_pos_neg.Rda')
