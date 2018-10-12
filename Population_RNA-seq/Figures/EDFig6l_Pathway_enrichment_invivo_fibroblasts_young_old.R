# Objective: Make graph for transcriptional signatures associated with age in vivo fibroblasts using Katja's script. 
# Will use the following pathways terms (donloaded from Enrichr): KEGG 2016 + activation
# Generate a graph visualizing the results

###################################################################################################
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Pathway_Enrichment/Pathway_enrichment_young_old_fibs_in_vivo/')

load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vivo/DESeq2_DEG_young_old_in_vivo_all_genes.Rda")
# Upload data
my.data <- as.data.frame(resOrdered)
# Make gene names in CAPTIAL letter
row.names(my.data) <- toupper(row.names(my.data))
# remove ows with NAs
my.data <- my.data[!is.na(my.data$stat),]

# FIRST! Upload function 
source("~/Dropbox/RNAseq/R-scripts/KEGG_pathway_enrichment.R")

#########################       Kegg Pathways
# Run function as stated, start with pathways that change with age (genes both up and down) using abs = T 
set.seed(12345)
kegg.results.change.w.age <- KEGG.pathway.enrichment(gene.sets.file = "~/Dropbox/Code_checking_SM/For_Brittany/Essentials/Pathway_genesets/KEGG_2016_disease_activation_UPPER_case.txt",
                                                     gene.sets.list = NULL,
                                                     test.results = my.data,
                                                     statistic = "stat",
                                                     fold.change = "log2FoldChange",
                                                     num.samples = 10000,
                                                     abs = FALSE)
sum(kegg.results.change.w.age$pathway.results$adj.p.val < 0.05)
# 25

# subset only the significant ones (cut-ff 0.01 as so many significant ones)
kegg.sig <- kegg.results.change.w.age$pathway.results[which(kegg.results.change.w.age$pathway.results$adj.p.val < 0.01),]
kegg.sig

# 3 order graph by fold change
kegg.sig.order <- kegg.sig[order(kegg.sig$fold.change, decreasing = T),]


# Make a barplot 
pdf("KEGG_pathway_enrichment.change.w.age_in_vivo_basal_180805.pdf", width = 8, height = 15)
par(mar=c(4.1, 25, 55.1, 7.2), xpd = TRUE)
x <- -kegg.sig.order$fold.change
col <- character(length = length(x))
col[x < 0] <- rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255)
col[x > 0] <- rgb(red=0,green=0,blue=139, alpha = 204, maxColorValue = 255)
b <- barplot(x,
             horiz = TRUE,
             names.arg = rownames(kegg.sig.order),
             xlab = "Average log2 fold change (Old / Young)",
             border = NA,
             las = 1,
             col = col,
             axes = FALSE,
             space = 0.05)
stars <- character(length = nrow(kegg.sig.order))
pvals <- kegg.sig.order$adj.p.val
stars[ pvals > 0.05 ] <- ""
stars[ pvals <= 0.05 ] <- "*"
stars[ pvals <= 0.01 ] <- "**"
stars[ pvals <= 0.001 ] <- "***"
pos <- integer(length = nrow(kegg.sig.order))
pos[ x < 0 ] <- 4
pos[ x > 0 ] <- 2
text(x = x,
     y = b[,1] - 0.3,
     stars,
     pos = pos)
axis(1,
     at = seq(-0.5, 0.5, by = 0.5))
dev.off()
