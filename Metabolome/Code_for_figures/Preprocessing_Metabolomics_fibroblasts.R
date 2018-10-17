# Objective: Pre-processing of metabolomic data from young and old fibroblasts at passage 3 for differential abundence, PCA and clustering

# DATA: Load table from xcms (pos and neg mode separately)
# Positive mode comes in triplicates and negative in duplicates, determine median for each feature based on the technical replicas and merge the tables 
# Pre-process by range-scaling (as recommended by Berg et al, BMC, 2006)

# Generated 3 tables:
# 1) Median Signal intensity NEG mode + median m/z (mzmed), median retention time (rtmed) and mode info (pos/neg) (will be used for differential analysis)
# 2) Median Signal intensity POS mode + median m/z (mzmed), median retention time (rtmed) and mode info (pos/neg) (will be used for differential analysis)
# 3) Range-scaled signal intensity pos and neg mode (will be used for clustering and PCA)

###########################################################################################
# set working directory
setwd(/Metabolome/RAW_data/")

# Process pos and neg mode separately and merge them into one table.

###########################################################################################
# Neg mode
# Upload neg mode data, 
data_n = read.csv("/Metabolome/RAW_data/XCMS_neg_mode_fibroblasts.csv") 

# NEG MODE: choose columns for data and make an index list that will be used to merge the technical replicas  
colnames(data_n)
data_neg = data_n[,grep("3M|29M", colnames(data_n))]
data_neg = data_neg[,-grep("_neg", colnames(data_neg))]
colnames(data_neg)
row.names(data_neg) = data_n[, 2]

# Create an index file, note neg mode was ran in technical duplicates
index_n = as.factor(c(rep("3m1",2),rep("3m2",2),rep("3m3",2),rep("3m4",2),rep("3m5",2),rep("3m6",2),rep("3m7",2),rep("3m8",2),
                      rep("29m1",2),rep("29m2",2),rep("29m3",2),rep("29m4",2),rep("29m5",2),rep("29m6",2),rep("29m7",2),rep("29m8",2)))

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
median_n.all <- cbind(data_n[,c(55,9,12)],median_n)
colnames(median_n.all)

# write/save table 
write.table(median_n.all, file="Metabolome_median_signal_intensity_NEG.txt", quote =F, sep="\t")
save(median_n.all, file="Metabolome_median_signal_intensity_NEG.Rda")

###########################################################################################
# Pos mode
# Upload Pos mode data
data_p = read.csv("/Metabolome/RAW_data/XCMS_pos_mode_fibroblasts.csv") 

colnames(data_p)
data_pos = data_p[,grep("3M|29M", colnames(data_n))]
data_pos = data_pos[,-grep("_pos", colnames(data_pos))]
colnames(data_pos)
row.names(data_pos) = data_p[, 2]

# Create an index file, note POS mode was ran in technical triplicates except 3m5
index_p = as.factor(c(rep("3m1",2),rep("3m2",2),rep("3m3",2),rep("3m4",2),rep("3m5",2),rep("3m6",2),rep("3m7",2),rep("3m8",2),
                      rep("29m1",2),rep("29m2",2),rep("29m3",2),rep("29m4",2),rep("29m5",2),rep("29m6",2),rep("29m7",2),rep("29m8",2)))

                     
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
median_p.all <- cbind(data_p[,c(55,9,12)],median_p)
colnames(median_p.all)

# write/save table 
write.table(median_p.all, file="Metabolome_median_signal_intensity_POS.txt", quote =F, sep="\t")
save(median_p.all, file="Metabolome_median_signal_intensity_POS.Rda")


################## Merge Pos and Neg mode, replace all 0 with a lowest value in the dataset, Range scale
# Merge the samples
data_t = rbind(median_p, median_n)

# Determine the number of 0 values 
length(which(data_t == 0))  #156

# Detrmine the minimum value that is not 0 
min.value <- min(data_t[data_t>0]) # 68.05737

# Replacing all zero values with the lowest value in the dataset
data_nz <- data_t
for (i in 1:nrow(data_t)) {
  for (j in 1:ncol(data_t)) {
    if (data_nz[i,j] == 0) {
      data_nz[i,j] = min.value
    }
  }
}
min(data_nz)

#Pre-processing of data; range-scaling (mean centered and divided by the range of each variable) 
data_nz.mean <- apply(data_t,1,mean)
data_nz.centered <- data_t - data_nz.mean
data_nz.max <- apply(data_nz.centered,1,max)
data_nz.min <- apply(data_nz.centered,1,min)
data_nz.range <- data_nz.max - data_nz.min
data_nz.norm <- data_nz.centered/data_nz.range

# Write and save table with Log10 transformed, range scaled data
write.table(data_nz.norm, file="Metabolome_range-scaled_values_pos_neg.txt", quote =F, sep="\t")
save(data_nz.norm, file="Metabolome_range-scaled_values_pos_neg.Rda")
