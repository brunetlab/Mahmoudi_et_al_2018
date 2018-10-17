# Objective: To pre-process metabolomic data from young and old iPSC for PCA and clustering

# DATA: Upload table from xcms (pos and neg mode separately) provided by Xiyan
# Positive mode comes in triplicates and negative in duplicates, determine median for each feature based on the technical replicas, 
# and merge the tables 

# Pre-process by range-scaling (as recommended by Berg et al, BMC, 2006)

# Generated 3 tables:
# 1) Median Signal intensity NEG mode + median m/z (mzmed), median retention time (rtmed) and mode info (pos/neg) (will be used for differential analysis)
# 2) Median Signal intensity POS mode + median m/z (mzmed), median retention time (rtmed) and mode info (pos/neg) (will be used for differential analysis)
# 3) Range-scaled signal intensity pos and neg mode (will be used for clustering and PCA)

###########################################################################################
# set working directory
setwd("/Metabolome/")

# Process pos and neg mode separately and merge them into one table.

###########################################################################################
# Neg mode
# Upload neg mode data, 
data_n = read.csv("~/Desktop/Github/Metabolome/RAW_data/XCMS_neg_mode_iPSC.csv") 

# NEG MODE: choose columns for data and make an index list that will be used to merge the technical replicas  
colnames(data_n)
data_neg = data_n[,grep("3M|29M", colnames(data_n))]
data_neg = data_neg[,-grep("_neg", colnames(data_neg))]
colnames(data_neg)
row.names(data_neg) = data_n[, 2]

# Create an index file, note neg mode was ran in technical triplicates
index_n = as.factor(c(rep("iPSC_3m1A4",3),rep("iPSC_3m2C1",3),rep("iPSC_3m2D9",3),rep("iPSC_3m4E1",3),rep("iPSC_3m4E2",3),
                      rep("iPSC_29m2C1",3),rep("iPSC_29m2D8",3),rep("iPSC_29m2E4",3),rep("iPSC_29m6E1",3),rep("iPSC_29m6E4",3),rep("iPSC_29m6F7",3),rep("iPSC_29m7G7",3),rep("iPSC_29m7H9",3)))

# Calculate the median for each metabolite for each sample
# make a loop where you look at the data frame row-wise, split by label (index list), apply median on each list, 
# return median and transpose data frame 
median_n.t <- apply(data_neg, 1, function(x) {
  y <- split(x, index_n)
  z <- sapply(y, median)
  return(z)
})
median_n <- t(median_n.t)

# Create a table with ionization mode (pos/neg), median m/z (mzmed), median retention time (rtmed), followed by median signal intensities of all features
colnames(data_n)
median_n.all <- cbind(data_n[,c(62,9,12)],median_n)
colnames(median_n.all)

# write/savle table 
write.table(median_n.all, file="Metabolome_iPSC_median_signal_intensity_NEG.txt", quote =F, sep="\t")
save(median_n.all, file="Metabolome_iPSC_median_signal_intensity_NEG.Rda")

###########################################################################################
# Pos mode
# Upload Pos mode data
data_p = read.csv("~/Desktop/Github/Metabolome/RAW_data/XCMS_pos_mode_iPSC.csv") 

colnames(data_p)
data_pos = data_p[,grep("3M|29M", colnames(data_p))]
data_pos = data_pos[,-grep("_pos", colnames(data_pos))]
colnames(data_pos)
row.names(data_pos) = data_p[, 2]

# Create an index file, note POS mode was ran in technical triplicates except 3m5
index_p = as.factor(c(rep("iPSC_3m1A4",3),rep("iPSC_3m2C1",3),rep("iPSC_3m2D9",3),rep("iPSC_3m4E1",3),rep("iPSC_3m4E2",3),
                      rep("iPSC_29m2C1",3),rep("iPSC_29m2D8",3),rep("iPSC_29m2E4",3),rep("iPSC_29m6E1",3),rep("iPSC_29m6E4",3),rep("iPSC_29m6F7",3),rep("iPSC_29m7G7",3),rep("iPSC_29m7H9",3)))

                     
# Calculate the median for each metabolite for each sample
# make a loop where you look at the data frame row-wise, split by label (index list), apply median on each list, 
# return median and transpose data frame 

median_p.t <- apply(data_pos, 1, function(x) {
  y <- split(x, index_p)
  z <- sapply(y, median)
  return(z)
})
median_p <- t(median_p.t)

# Create a table with ionization mode (pos/neg), median m/z (mzmed), median retention time (rtmed), followed by median signal intensities of all features
colnames(data_p)
median_p.all <- cbind(data_p[,c(62,9,12)],median_p)
colnames(median_p.all)

# write/savle table 
write.table(median_p.all, file="Metabolome_iPSC_median_signal_intensity_POS.txt", quote =F, sep="\t")
save(median_p.all, file="Metabolome_iPSC_median_signal_intensity_POS.Rda")


################## Merge Pos and Neg mode, replace all 0 with a low number, log2 normalize and range scale
# Merge the samples
data_t = rbind(median_p[,1:13], median_n[,1:13])

#Pre-processing of data; range-scaling (mean centered and divided by the range of each variable) 
data_nz.mean <- apply(data_t,1,mean)
data_nz.centered <- data_t - data_nz.mean
data_nz.max <- apply(data_nz.centered,1,max)
data_nz.min <- apply(data_nz.centered,1,min)
data_nz.range <- data_nz.max - data_nz.min
data_nz.norm <- data_nz.centered/data_nz.range

# Write and save table with Log10 transformed, range scaled data
write.table(data_nz.norm, file="Metabolome_iPSC_range-scaled_values_pos_neg.txt", quote =F, sep="\t")
save(data_nz.norm, file="Metabolome_iPSC_range-scaled_values_pos_neg.Rda")
