# Objective: Make barplot that visualized the the transcriptional pathways that are associated with RE based on linear regresional model and top vs bottom5 RE
# Determine the pathways that are comon between the folwoin two analyses:
## Extended_Data_Table_x_pathway_enrichment_linear_modeling_RE
## Extended_Data_Table_x_pathway_enrichment_Top5_vs_bottom5_RE

###################################################################################################
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Figures/Pathway_Enrichment/Fig2g_Pathways_associated_with_RE/')

# load data 
linear.modeling <- read.csv("~/Dropbox/Code_checking_SM/For_Brittany/Tables/Extended_Data_Table_x_pathway_enrichment_linear_modeling_RE/kegg.results.pathways.linear_modeling_young_old_RE_180805.csv",
                            stringsAsFactors = F, row.names = 1)

linear.modeling.sig <-  linear.modeling[which(linear.modeling$adj.p.val < 0.05),]

tb5 <- read.csv("~/Dropbox/Code_checking_SM/For_Brittany/Tables/Extended_Data_Table_x_pathway_enrichment_Top5_vs_bottom5_RE/kegg.results.pathways.TopBottom5_RE_young_old_fibs_180805.csv",
                stringsAsFactors = F, row.names = 1)
tb5.sig <-  tb5[which(tb5$adj.p.val < 0.05),]

# Identify common pathways 
common.pathways <- linear.modeling.sig[which(row.names(linear.modeling.sig) %in% row.names(tb5.sig)),]

#Order pathways based on fold change
common.pathways.order <- common.pathways[order(common.pathways$fold.change, decreasing = F),]
common.pathways.order

# Make a barplot 
pdf("KEGG_pathway_enrichment.associated_with_RE_180805.pdf", width = 8, height = 15)
par(mar=c(4.1, 25, 60.1, 7.2), xpd = TRUE)
x <- common.pathways.order$fold.change
col <- character(length = length(x))
col[x < 0] <- rgb(red=255,green=165,blue=0, alpha = 204, maxColorValue = 255)
col[x > 0] <- rgb(red=0,green=0,blue=139, alpha = 204, maxColorValue = 255)
b <- barplot(x,
             horiz = TRUE,
             names.arg = rownames(common.pathways.order),
             xlab = "Average log2 fold change (Good to Bad reprogramming)",
             border = NA,
             las = 1,
             col = col,
             axes = FALSE,
             space = 0.1)
stars <- character(length = nrow(common.pathways.order))
pvals <- common.pathways.order$adj.p.val
stars[ pvals > 0.05 ] <- ""
stars[ pvals <= 0.05 ] <- "*"
stars[ pvals <= 0.01 ] <- "**"
stars[ pvals <= 0.001 ] <- "***"
pos <- integer(length = nrow(common.pathways.order))
pos[ x < 0 ] <- 4
pos[ x > 0 ] <- 2
text(x = x,
     y = b[,1] - 0.3,
     stars,
     pos = pos)
axis(1,
     at = seq(-0.2, 0.2, by = 0.05))
dev.off()
