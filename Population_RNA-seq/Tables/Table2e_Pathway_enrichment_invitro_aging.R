# Objective: Look for gene set enrichment in genes that change in old upon wounding using Katja's script. 
# Will use the following pathways terms (downloaded from Enrichr): 
##### KEGG 2016 + activation

# DATA: DESeq2 data from DESeq2_analysis_DEGs_VST_values_Linear_modeling_RE
# Generates a table with all the results
###################################################################################################
setwd('~/Dropbox/Code_checking_SM/For_Brittany/Tables/Katja_Extended_Data_Table_young_old_invitro/')

# Load data
load("~/Dropbox/Code_checking_SM/For_Brittany/Core_Analysis/Fibroblast_in_vitro/DESeq2_DEG_all_genes.Rda")
my.data <- as.data.frame(resOrdered)
# Make gene names in CAPTIAL letter
row.names(my.data) <- toupper(row.names(my.data))
# remove ows with NAs
my.data <- my.data[!is.na(my.data$stat),]

# FIRST! Upload function 
source("~/Dropbox/Code_checking_SM/For_Brittany/Essentials/R-scripts/KEGG_pathway_enrichment3.R")

#########################       Kegg Pathways
# Run function as stated
set.seed(12345)
kegg.results.young_old_invitro <- KEGG.pathway.enrichment(gene.sets.file = "~/Dropbox/Code_checking_SM/For_Brittany/Essentials/Pathway_genesets/KEGG_2016_disease_activation_UPPER_case.txt",
                                                     gene.sets.list = NULL,
                                                     test.results = my.data,
                                                     statistic = "stat",
                                                     fold.change = "log2FoldChange",
                                                     num.samples = 10000,
                                                     abs = FALSE)
sum(kegg.results.young_old_invitro$pathway.results$adj.p.val < 0.05)
# 59
save(kegg.results.young_old_invitro,
      file = "kegg.results.pathways.kegg.results.young_old_invitro.RData", row.names = T)
write.csv(kegg.results.young_old_invitro$pathway.results,
            file= "kegg.results.pathways.kegg.results.young_old_invitro.csv")

