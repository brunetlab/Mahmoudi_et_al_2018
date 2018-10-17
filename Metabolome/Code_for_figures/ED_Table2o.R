# Objective: To determine and generate a table of significantly differently abundant metabolite features
# Extended Data Table 2o
# Use processed data from Pre-processing_Metabolomics_Fibroblasts.R

####################################################################################
# Required R packages
library(qvalue) # qvalue_2.6.0 
sessionInfo()
####################################################################################
# set working directory
setwd("Metabolome/")

# In a loop calculate mean 3m, mean 29m, mean across all samples, Log2FoldChange	Log2 fold change (Old Fib vs Young Fib), 
# determine significance between young and old using a non-parametric test (wilcoxon-test), 
# finally calculate FDR (BH) 

# Upload Neg mode data
load("/Metabolome/RAW_data/Metabolome_median_signal_intensity_NEG.Rda") 
colnames(median_n.all)

median_n.all.stat <- as.data.frame(median_n.all)
for (i in 1:nrow(median_n.all)) {
  x <- as.numeric(median_n.all[i,grep("29m", colnames(median_n.all))])
  y <- as.numeric(median_n.all[i,grep("3m", colnames(median_n.all))])
  median_n.all.stat$mean.3m[i] <- mean(y)
  median_n.all.stat$mean.29m[i] <- mean(x)
  median_n.all.stat$baseMean[i] <- mean(c(x,y))
  median_n.all.stat$Log2FC[i] <- log2((mean(x)/mean(y)))
  wx <- wilcox.test(x,y, alternative="two.sided")
  median_n.all.stat$Wilcoxon[i] <- wx$p.value
}
median_n.all.stat$p.adjusted.BH <- p.adjust(median_n.all.stat$Wilcoxon, "BH")
q.n <- qvalue(median_n.all.stat$Wilcoxon)
median_n.all.stat$p.adjusted.Qvalues <- q.n$qvalues

# Determine how many features have BH < 0.05
length(which(median_n.all.stat$p.adjusted.BH < 0.05)) # 343
length(which(median_n.all.stat$p.adjusted.Qvalues < 0.05)) # 343

save(median_n.all.stat, file='Differential_analysis_fibroblast_metabolome_neg.Rda')


####################################################################################

# Upload Pos mode data
load("/Metabolome/RAW_data/Metabolome_median_signal_intensity_POS.Rda")
colnames(median_p.all)

median_p.all.stat <- as.data.frame(median_p.all)
for (i in 1:nrow(median_p.all)) {
  x <- as.numeric(median_p.all[i,grep("29m", colnames(median_n.all))])
  y <- as.numeric(median_p.all[i,grep("3m", colnames(median_n.all))])
  median_p.all.stat$mean.3m[i] <- mean(y)
  median_p.all.stat$mean.29m[i] <- mean(x)
  median_p.all.stat$baseMean[i] <- mean(c(x,y))
  median_p.all.stat$Log2FC[i] <- log2((mean(x)/mean(y)))
  wx <- wilcox.test(x,y, alternative="two.sided")
  median_p.all.stat$Wilcoxon[i] <- wx$p.value
}
median_p.all.stat$p.adjusted.BH <- p.adjust(median_p.all.stat$Wilcoxon, "BH")
q.p <- qvalue(median_p.all.stat$Wilcoxon)
median_p.all.stat$p.adjusted.Qvalues <- q.p$qvalues

# Determine how many features have BH < 0.05
length(which(median_p.all.stat$p.adjusted.BH < 0.05)) # 313
length(which(median_p.all.stat$p.adjusted.Qvalues < 0.05)) # 394

save(median_p.all.stat, file='Differential_analysis_fibroblasts_metabolome_pos.Rda')

# Merge the two tables together including all data
data_all <- rbind(median_p.all.stat,median_n.all.stat)
write.csv(data_all, file='Differential_analysis_fibroblasts_metabolome_pos_neg.csv')
save(data_all, file='Differential_analysis_fibroblasts_metabolome_pos_neg.Rda')



## Extended Data Table 20
# Extract Mode, mxmed, rtmed, baseMean,Log2FC,Wilcoxon, p.adjusted.BH info from merged table
colnames(data_all)
data_c = data_all[,c(1:3,22:26)]
colnames(data_c)
#  "Mode"               "mzmed"              "rtmed"              "baseMean"           "Log2FC"             "Wilcoxon"           "p.adjusted.BH"      "p.adjusted.Qvalues"

# Extract only the significant ones (FDR <0.05 and FC > 1.5)
data_c.sig <- data_c[which(data_c$p.adjusted.Qvalues < 0.05 & abs(data_c$Log2FC) > 0.58),] # 484
data_c.sig.up <- data_c[which(data_c$p.adjusted.Qvalues < 0.05 & data_c$Log2FC > 0.58),] # 295
data_c.sig.dn <- data_c[which(data_c$p.adjusted.Qvalues < 0.05 & data_c$Log2FC < -0.58),] # 189
write.csv(data_c.sig, file='Extended_Data_table2b.csv')




