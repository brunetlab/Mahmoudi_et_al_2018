# Objective: Generate plots for gene expression of genes related to Cytokine-cytokine receptor signaling/TNF signaling/Fibroblast activation between young and old across different subgroups
# Figure 4x
# Use data generated on cluster using Step3_pagoda.. :
## young_old_norm_counts_gene5500_mito10perc.Rdata (gene counts)

################################################################################################################
# packages required
library(ggplot2)
sessionInfo()
################################################################################################################
# upload files required 
setwd("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/Transcriptome_signatures_youngold_PAGODA_clusters/Expected_data/")
load("~/Dropbox/Code_checking_SM/For_Matt/Essentials/Expression_data_PAGODA/young_old_norm_counts_gene5500_mito10perc.Rdata")
load("~/Dropbox/Code_checking_SM/For_Matt/PAGODA/PAGODA_Visualization/Visualization_KEGG_disease_aging/Expected_data/cluster.cells.top.aspects_KEGG_DISEASE_AGING_varnorm_knn.error.models.repeat_set_seed.RData")

# Upload kegg pathways 
gene.sets <- read.table("~/Dropbox/Code_checking_SM/For_Matt/Essentials/Pathways/KEGG_2016_disease_aging_activation_lower_case_for_PAGODA_090518.txt", sep="\t", header = F, quote="",row.names = 1)
gene.set.data <- as.data.frame(t(gene.sets))
gene.expression <- as.data.frame(young_old_norm_counts_gene5500_mito10perc)

# Input color information to be used by PAGODA
age.info <- c(rep("A_young",length(grep("Young", colnames(young_old_norm_counts_gene5500_mito10perc)))),
              rep("B_old",length(grep("Old", colnames(young_old_norm_counts_gene5500_mito10perc)))))

# Define the different subpopulations (6 main ones based on PAGODA), and rename so that order is the same as on PAGODA heatmap
sg <- as.character(cutree(hc, k=6))
sg[sg == 1] <- "f1"
sg[sg == 2] <- "c2"
sg[sg == 3] <- "a3"
sg[sg == 4] <- "d4"
sg[sg == 5] <- "b5"
sg[sg == 6] <- "e6"

# Define color for jitter (dots) for plots
am = col=rgb(red=0,green=0,blue=0, alpha = 204, maxColorValue = 255)

################################################################################################################
# Extract TNF signaling interaction genes and determine median expression across all genes
genes.tnf <- gene.set.data$`TNF signaling pathway`
genes.tnf <- na.omit(genes.tnf)

# Extract relevant genes
expression.tnf <- gene.expression[which(toupper(row.names(gene.expression)) %in% toupper(genes.tnf)),]
row.names(expression.tnf)

# Determine SUM expression of genes of interest
expression.tnf.sum <- as.data.frame(colSums(expression.tnf))

# input age and subgroup
expression.tnf.sum$Age <- age.info
expression.tnf.sum$Subgroup <- sg

# double check everything looks as expected
head(expression.tnf.sum)
tail(expression.tnf.sum)

# ggplot2 plot  
p <- ggplot(expression.tnf.sum, 
            aes(x=sg, 
                y=expression.tnf.sum$`colSums(expression.tnf)`, 
                fill=Age)) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.7, aes(colour=am)) +
  geom_violin(draw_quantiles = c(0.5), width=0.7, size =0.8) +
  scale_fill_manual(values=c("orange", "darkblue"))

max.exp <- ceiling(max(expression.tnf.sum$`colSums(expression.tnf)`))

my.plot <- p + 
  ggtitle("TNF signaling pathway") +
  xlab("Subgroup") +
  ylab("Expression") +
  theme_classic() +
  scale_y_continuous(limits = c(-1, max.exp))

my.plot

# save plot
ggsave(filename="TNF_signaling_pathway_young_old.pdf", plot=my.plot, width = 3.5, height = 2)

# Test statistical significance in different subgroups between young and old - use non-parametric wilcoxon test
# subset different subgroups and perform stats
sg1.tnf <- expression.tnf.sum[grep("f1", expression.tnf.sum$Subgroup),]
sg1.tnf <- wilcox.test(sg1.tnf$`colSums(expression.tnf)` ~ sg1.tnf$Age)
# W = 38315, p-value = 0.7205

sg2.tnf <- expression.tnf.sum[grep("c2", expression.tnf.sum$Subgroup),]
sg2.tnf <- wilcox.test(sg2.tnf$`colSums(expression.tnf)` ~ sg2.tnf$Age)
# W = 19713, p-value = 0.708

sg3.tnf <- expression.tnf.sum[grep("a3", expression.tnf.sum$Subgroup),]
sg3.tnf <- wilcox.test(sg3.tnf$`colSums(expression.tnf)` ~ sg3.tnf$Age)
# W = 28761, p-value = 0.01851

sg4.tnf <- expression.tnf.sum[grep("d4", expression.tnf.sum$Subgroup),]
sg4.tnf <- wilcox.test(sg4.tnf$`colSums(expression.tnf)` ~ sg4.tnf$Age)
# W = 118580, p-value = 0.4219

sg5.tnf <- expression.tnf.sum[grep("b5", expression.tnf.sum$Subgroup),]
sg5.tnf <- wilcox.test(sg5.tnf$`colSums(expression.tnf)` ~ sg5.tnf$Age)
# W = 23873, p-value = 0.4412

sg6.tnf <- expression.tnf.sum[grep("e6", expression.tnf.sum$Subgroup),]
sg6.tnf <- wilcox.test(sg6.tnf$`colSums(expression.tnf)` ~ sg6.tnf$Age)
#W = 4026, p-value = 0.4127

#Reorder based on graph
p.adjust(c(sg3.tnf$p.value, sg5.tnf$p.value, sg2.tnf$p.value, sg4.tnf$p.value, sg6.tnf$p.value, sg1.tnf$p.value), method = "BH")
# 0.1110461 0.6617291 0.7204793 0.6617291 0.6617291 0.7204793

################################################################################################################
################################################################################################################
# Extract Fibroblast activation genes and determine median expression across all genes
genes.FA <- gene.set.data$`Fibroblast activation`
genes.FA <- na.omit(genes.FA)

# Extract relevant genes 
expression.FA <- gene.expression[which(toupper(row.names(gene.expression)) %in% toupper(genes.FA)),]
expression.FA <- na.omit(expression.FA)
row.names(expression.FA)

# Determine SUM expression of genes of interest
expression.FA.sum <- as.data.frame(colSums(expression.FA))

# input age and subgroup
expression.FA.sum$Age <- age.info
expression.FA.sum$Subgroup <- sg

# double check everything looks as expected
head(expression.FA.sum)
tail(expression.FA.sum)

# ggplot2 plot 
p <- ggplot(expression.FA.sum, 
            aes(x=sg, 
                y=expression.FA.sum$`colSums(expression.FA)`, 
                fill=Age)) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.7, aes(colour=am)) +
  geom_violin(draw_quantiles = c(0.5), width=0.7, size =0.8) +
  scale_fill_manual(values=c("orange", "darkblue"))


max.exp <- ceiling(max(expression.FA.sum$`colSums(expression.FA)`))

my.plot <- p + 
  ggtitle("Fibroblast activation") +
  xlab("Subgroup") +
  ylab("Expression") +
  theme_classic() +
  scale_y_continuous(limits = c(-1, max.exp))

my.plot

# save plot
ggsave(filename="Fibroblast_activation_young_old.pdf", plot=my.plot, width = 3.5, height = 2)

# Test statistical significance in different subgroups between young and old - use non-parametric wilcoxon test
# subset different subgroups and perform stats
sg1.FA <- expression.FA.sum[grep("f1", expression.FA.sum$Subgroup),]
sg1.FA <- wilcox.test(sg1.FA$`colSums(expression.FA)` ~ sg1.FA$Age)
# W = 35742, p-value = 0.0886

sg2.FA <- expression.FA.sum[grep("c2", expression.FA.sum$Subgroup),]
sg2.FA <- wilcox.test(sg2.FA$`colSums(expression.FA)` ~ sg2.FA$Age)
# W = 18089, p-value = 0.2857

sg3.FA <- expression.FA.sum[grep("a3", expression.FA.sum$Subgroup),]
sg3.FA <- wilcox.test(sg3.FA$`colSums(expression.FA)` ~ sg3.FA$Age)
# W = 22500, p-value = 0.03287

sg4.FA <- expression.FA.sum[grep("d4", expression.FA.sum$Subgroup),]
sg4.FA <- wilcox.test(sg4.FA$`colSums(expression.FA)` ~ sg4.FA$Age)
# W = 88224, p-value = 9.914e-10

sg5.FA <- expression.FA.sum[grep("b5", expression.FA.sum$Subgroup),]
sg5.FA <- wilcox.test(sg5.FA$`colSums(expression.FA)` ~ sg5.FA$Age)
# W = 22507, p-value = 0.7758

sg6.FA <- expression.FA.sum[grep("e6", expression.FA.sum$Subgroup),]
sg6.FA <- wilcox.test(sg6.FA$`colSums(expression.FA)` ~ sg6.FA$Age)
# W = 3429, p-value = 0.3455

p.adjust(c(sg3.FA$p.value, sg5.FA$p.value, sg2.FA$p.value, sg4.FA$p.value, sg6.FA$p.value, sg1.FA$p.value), method = "BH")
#9.859737e-02 7.757718e-01 4.145501e-01 5.948481e-09 4.145501e-01 1.771965e-01


################################################################################################################
################################################################################################################
# Extract Fibroblast activation genes and determine median expression across all genes
genes.CK <- gene.set.data$`Cytokine-cytokine receptor interaction`
genes.CK <- na.omit(genes.CK)

# Extract relevant genes 
expression.CK <- gene.expression[which(toupper(row.names(gene.expression)) %in% toupper(genes.CK)),]
expression.CK <- na.omit(expression.CK)
row.names(expression.CK)

# Determine SUM expression of genes of interest
expression.CK.sum <- as.data.frame(colSums(expression.CK))

# input age and subgroup
expression.CK.sum$Age <- age.info
expression.CK.sum$Subgroup <- sg

# double check everything looks as expected
head(expression.CK.sum)
tail(expression.CK.sum)

# ggplot2 plot 
p <- ggplot(expression.CK.sum, 
            aes(x=sg, 
                y=expression.CK.sum$`colSums(expression.CK)`, 
                fill=Age)) + 
  geom_jitter(shape=16, position=position_jitter(0.3), alpha=.7, aes(colour=am)) +
  geom_violin(draw_quantiles = c(0.5), width=0.7, size =0.8) +
  scale_fill_manual(values=c("orange", "darkblue"))


max.exp <- ceiling(max(expression.CK.sum$`colSums(expression.CK)`))

my.plot <- p + 
  ggtitle("Cytokine-cytokine receptor interaction") +
  xlab("Subgroup") +
  ylab("Expression") +
  theme_classic() +
  scale_y_continuous(limits = c(-1, max.exp))

my.plot

# save plot
ggsave(filename="Cytokine-cytokine_receptor_interaction_young_old.pdf", plot=my.plot, width = 3.5, height = 2)

# Test statistical significance in different subgroups between young and old - use non-parametric wilcoxon test
# subset different subgroups and perform stats
sg1.CK <- expression.CK.sum[grep("f1", expression.CK.sum$Subgroup),]
sg1.CK <- wilcox.test(sg1.CK$`colSums(expression.CK)` ~ sg1.CK$Age)
# W = 35309, p-value = 0.05371

sg2.CK <- expression.CK.sum[grep("c2", expression.CK.sum$Subgroup),]
sg2.CK <- wilcox.test(sg2.CK$`colSums(expression.CK)` ~ sg2.CK$Age)
# W = 21328, p-value = 0.07039

sg3.CK <- expression.CK.sum[grep("a3", expression.CK.sum$Subgroup),]
sg3.CK <- wilcox.test(sg3.CK$`colSums(expression.CK)` ~ sg3.CK$Age)
# W = 26575, p-value = 0.4309

sg4.CK <- expression.CK.sum[grep("d4", expression.CK.sum$Subgroup),]
sg4.CK <- wilcox.test(sg4.CK$`colSums(expression.CK)` ~ sg4.CK$Age)
# W = 117120, p-value = 0.638

sg5.CK <- expression.CK.sum[grep("b5", expression.CK.sum$Subgroup),]
sg5.CK <- wilcox.test(sg5.CK$`colSums(expression.CK)` ~ sg5.CK$Age)
# W = 21935, p-value = 0.4672

sg6.CK <- expression.CK.sum[grep("e6", expression.CK.sum$Subgroup),]
sg6.CK <- wilcox.test(sg6.CK$`colSums(expression.CK)` ~ sg6.CK$Age)
# W = 3705, p-value = 0.8988

p.adjust(c(sg3.CK$p.value, sg5.CK$p.value, sg2.CK$p.value, sg4.CK$p.value, sg6.CK$p.value, sg1.CK$p.value), method = "BH")
# 0.7008732 0.7008732 0.2111667 0.7655836 0.8988047 0.2111667
