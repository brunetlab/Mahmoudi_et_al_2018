# Objective: Generate a strip plot comparing the DEGs in fibroblast to DEGs in iPSCs
# Extended Data Figure 3d
# Use a merged table of DEG for FIb and iPSC, with information on cell type and color codes:
# 2 = grey, 3 = darkblue, 4 = orange, 5 = gold1, 6 = dodgerblue4


###########################################################################################
# Required package
library("DESeq2")
library("lattice")
###########################################################################################
# Set working directory
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Stripplot/Striplot_transcriptome/')

# Load data from Fibroblasts analysis
load('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro/DESeq2_DEG_all_genes.Rda')
my.fibs <- as.data.frame(resOrdered)
my.fibs$Cell.type <- "Fibroblast (pooled)"
# Reverse log2FoldChange so that + mean up in old
my.fibs$log2FoldChange <- -my.fibs$log2FoldChange
head(my.fibs)

# Substract fib non-degs and mark with 1
my.fibs.nondeg <- my.fibs[-which(my.fibs$padj < 0.05 & abs(my.fibs$log2FoldChange) > 0.58),]
my.fibs.nondeg$sig <- 1
head(my.fibs.nondeg)

# Substract fib DEGs that go down with age and mark with 2
my.fibs.deg.up  <- my.fibs[which(my.fibs$padj < 0.05 & my.fibs$log2FoldChange > 0.58),]  
my.fibs.deg.up$sig <- 2
head(my.fibs.deg.up)

# Substract fib DEGs that go up with age and mark with 3
my.fibs.deg.down <- my.fibs[which(my.fibs$padj < 0.05 & my.fibs$log2FoldChange < -0.58),]  
my.fibs.deg.down$sig <- 3
head(my.fibs.deg.down)

# Load data from iPSC analysis
load('~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/iPSC_in_vitro/DESeq2_iPSC_vst_values.Rda')
my.ipsc <- as.data.frame(resOrdered)
my.ipsc$Cell.type <- "iPSC (pooled)"
head(my.ipsc)

# Substract fib non-degs and mark with 1
my.ipsc.nondeg <- my.ipsc[-which(my.ipsc$padj < 0.05 & abs(my.ipsc$log2FoldChange) > 0),]
my.ipsc.nondeg$sig <- 1
head(my.ipsc.nondeg)

# Substract fib DEGs that go down with age and mark with 2
my.ipsc.deg.up  <- my.ipsc[which(my.ipsc$padj < 0.05 & my.ipsc$log2FoldChange > 0),]  
my.ipsc.deg.up$sig <- 2
head(my.ipsc.deg.up)

# Substract fib DEGs that go up with age and mark with 3
my.ipsc.deg.down <- my.ipsc[which(my.ipsc$padj < 0.05 & my.ipsc$log2FoldChange < -0),]  
my.ipsc.deg.down$sig <- 3
head(my.ipsc.deg.down)

# Merge all the data int one data frame 
my.data.all <- rbind(my.fibs.nondeg,
                     my.fibs.deg.up,
                     my.fibs.deg.down,
                     my.ipsc.nondeg,
                     my.ipsc.deg.up,
                     my.ipsc.deg.down)


# Color codes
col2rgb("lightgrey") # 211, 211, 211
col2rgb("darkblue") # 0, 0, 139
col2rgb("orange") # 255, 165, 0
col2rgb("gold1")   # 255, 215, 37
col2rgb("dodgerblue4") # 16, 78, 139


# Make the striplot
#pdf("EDF3e_Strip_plot_DEG_iPSC_FIB.pdf", width =4.3, height =5 )
png("EDF3e_Strip_plot_DEG_iPSC_FIB.png", width =600, height =800 )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
stripplot(log2FoldChange ~ Cell.type,
          data = my.data.all,
          jitter.data=TRUE,
          factor = 1,
          vertical = TRUE,
          ylab = "Log2 fold change",
          grid = "h",
          groups = my.data.all$sig,
          col = c(rgb(211, 211, 211, alpha = 204, maxColorValue = 255), #azure4
                  rgb(0, 0, 139, alpha = 204, maxColorValue = 255), # down with age fibs
                  rgb(255, 165, 0, alpha = 204, maxColorValue = 255) # up with age fibs
          ),
          pch = 17,
          cex=2
)
dev.off()

