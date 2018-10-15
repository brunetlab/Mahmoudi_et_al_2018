Pagoda analysis Katja

1. Take expression data (vst.values) from Salah
2. Make it into an expression set (script: import.R in the script folder)
3. Run KNN models on the cluster (script: knn.error.models.R) can also run on your own computer. The output of knn.error.models.R will then be used for the next step
4. Run the pagoda_KEGG_Aging_Activation.R until clpca is loaded.
5. To identify novel gene cluster, run gene.clusters.R.
6. Go back to pagoda_KEGG_Aging_Activation.R to finish the last of the script
7. Visualize the data using pagoda_visualization_KEGG_Aging_Activation.R. This will give you a heat map of the overdispersed pathways.
8. Identify pathways (de novo gene clusters) that may drive the signatures and may not be interesting in this case (e.g. cell cycle). Note! Do not need to make KNN models oridentify de novo gene clusters again. 
9. Account for these pathways using pagoda_subtractCellCycle_KEGG_Aging_Activation.R. This will provide a new varinfo, which you then use for the rest of the analysis.
10. To visualize the final data use pagoda_visualization_subtractCellCycle_KEGG_Aging_Activation.R.



To create heatmaps of single genes that uses the pagoda clustering use: heatmaps_subtractedCellCycle_KEGG_Aging_Activation.R

