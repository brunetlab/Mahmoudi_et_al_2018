setwd('/Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_height_fibro')

# 2016-03-02 - do differential height on Salah's data

library('DiffBind')
library("rtracklayer")

################################################################################
####################   Diffbind Analysis : H3K4me3 height  ####################
################################################################################

DF.aging <- dba(sampleSheet="Fibro_FIXSEQ_exp_K4me3_height1.csv",skipLines=1,attributes=c(DBA_ID,DBA_CONDITION))
DF.aging <- dba.count(DF.aging)


pdf("DF_aging_H3K4me3_height_clust_heatmap_Stitched.pdf")
plot(DF.aging, colScheme="Reds")
dev.off()

DF.aging <- dba.contrast(DF.aging, categories=DBA_CONDITION,minMembers=2)
DF.aging <- dba.analyze(DF.aging,method=DBA_ALL_METHODS)
DF.aging
# 
# 6 Contrasts:
#   Group1 Members1 Group2 Members2 DB.edgeR DB.DESeq DB.DESeq2
#1     3m        2    29m        3     6635      637      6979


DF.aging.3mvs29m <- dba.report(DF.aging,,method=DBA_EDGER, th=0.05, bLoss = TRUE, bGain = TRUE)
export(DF.aging.3mvs29m[DF.aging.3mvs29m$Fold < 0], "DiffBind_H3K4me3peaksDF_aging_3mvs29m_DBA_EDGER_FDR0.05_GAINED.bed")
export(DF.aging.3mvs29m[DF.aging.3mvs29m$Fold > 0], "DiffBind_H3K4me3peaksDF_aging_3mvs29m_DBA_EDGER_FDR0.05_LOST.bed")


pdf("MA_plot_DF_3m_vs_29m_H3K4me3peaksDiffbind_DBA_EDGER.pdf")
dba.plotMA(DF.aging, method=DBA_EDGER, th=0.05)
dev.off()

pdf("Diffbind_PCA_Total_DF_H3K4me3Peaks.pdf")
dba.plotPCA(DF.aging,DBA_CONDITION,method=DBA_EDGER,cex.pt=0.5)
dev.off()

pdf("Diffbind_Heatmap_DF_H3K4me3Peaks_allsamples_sig3v29m.pdf")
dba.plotHeatmap(DF.aging,mask=DF.aging$masks$All,method=DBA_EDGER,
                correlations=FALSE,scale="row", colScheme="RdYlBu", th = 0.05,
                ColAttributes = NULL,Colv=NULL,distMethod="correlation", margin = 15) 
dev.off()

library('pheatmap')
pdf("Manual_pheatmap_DF_H3K4me3Peaks_allsamples_sig3v29m.pdf")
pheatmap(DF.aging$allvectors[as.numeric(names(DF.aging.3mvs29m)),4:8],scale="row",cluster_cols = F,show_rownames = F,
         main="H3K4me3 intensity changed in 3 vs 29m",
         cellwidth = 50, cellheight = 0.05,border_color = NA,clustering_distance_rows = "correlation") 
dev.off()



#############
# S2 norm needed -> export !!!!

write.table(DF.aging$allvectors,file="DF_Diffbind_matrix_concentration.txt",quote=F,row.names=F,sep="\t")

pdf("concentration_boxplot_DF_aging.pdf")
boxplot(DF.aging$allvectors[,4:8],log='y',las=2)
dev.off()
