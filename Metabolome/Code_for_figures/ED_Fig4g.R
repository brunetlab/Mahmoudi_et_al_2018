# Objective: Generate a strip plot comparing the DEGs in fibroblast to iPSCs
# Extended Data Figure 4g
# Use a merged table of DEG for FIb and iPSC, with information on cell type and color codes:
# 2 = grey, 3 = dakblue, 4 = orange, 5 = gold1, 6 = dodgerblue4

###########################################################################################
# Required R package
library(lattice)
###########################################################################################
# Set working directory
setwd("Metablome/")

# Upload a merged table of DEG for FIb and iPSC, cell type and col.sign
my.data <- read.csv('Differential_analysis_metabolome_pos_neg_FIB_iPSC.csv', stringsAsFactors = F)

# Color codes
col2rgb("lightgrey") # 211, 211, 211
col2rgb("darkblue") # 0, 0, 139
col2rgb("orange") # 255, 165, 0
#col2rgb("gold1")   # 255, 215, 37
#col2rgb("dodgerblue4") # 16, 78, 139


# Make the striplot
pdf("ED_Fig4g.pdf", width =4.3, height =5 )
par(mar=c(4.1,4.1,6,1), xpd=TRUE)
stripplot(log2FoldChange ~ Cell.type,
          data = my.data,
          jitter.data=TRUE,
          factor = 1,
          vertical = TRUE,
          ylab = "Log2 fold change (old/young)",
          grid = "h",
          groups = col.sign,
          # darkblue, gold1
          col = c(
            rgb(211, 211, 211, alpha = 204, maxColorValue = 255), #lightgrey
            rgb(255, 165, 0, alpha = 204, maxColorValue = 255), # up with age fibs
            rgb(0, 0, 139, alpha = 204, maxColorValue = 255) # down with age fibs
          ),
          pch = 15
)
dev.off()