setwd('/Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth_old_method/')
source('get_cov_hists.r')

# 2016-07-28
# Old method for Salah
density.per.peak.3m1 <- read.csv('3m1_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS_cov_regions.bed',header=F,sep="\t")
density.per.peak.3m2 <- read.csv('3m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS_cov_regions.bed',header=F,sep="\t")
density.per.peak.29m2 <- read.csv('29m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS_cov_regions.bed',header=F,sep="\t")
density.per.peak.29m6 <- read.csv('29m6_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS_cov_regions.bed',header=F,sep="\t")
density.per.peak.29m7 <- read.csv('29m7_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS_cov_regions.bed',header=F,sep="\t")

my.hist.3m1 <- density.per.peak.3m1[which(density.per.peak.3m1$V1 == 'all'),]
my.hist.3m2 <- density.per.peak.3m2[which(density.per.peak.3m2$V1 == 'all'),]
my.hist.29m2 <- density.per.peak.29m2[which(density.per.peak.29m2$V1 == 'all'),]
my.hist.29m6 <- density.per.peak.29m6[which(density.per.peak.29m6$V1 == 'all'),]
my.hist.29m7 <- density.per.peak.29m7[which(density.per.peak.29m7$V1 == 'all'),]

pdf("2016-07-28_Hists_DF_Pre_adj_on_mergedCalling.pdf")
plot(0:39,my.hist.3m1$V5[1:40],type='l',col="darkblue",ylim=c(0,0.3),xlab="Coverage (X)",ylab="Histograms of H3K4me3 reference peak coverage")
lines(0:39,my.hist.3m2$V5[1:40],type='l',col="cadetblue")
lines(0:39,my.hist.29m2$V5[1:40],type='l',col="lightgreen")
lines(0:39,my.hist.29m6$V5[1:40],type='l',col="darkgreen")
lines(0:39,my.hist.29m7$V5[1:40],type='l',col="green")
legend("bottomright",c("3 mths (rep 1)","3 mths (rep 2)","29 mths (rep 2)","29 mths (rep 6)","29 mths (rep 7)"),
       col=c("darkblue","cadetblue","lightgreen","darkgreen","green"),bty='n',pch='_',cex = 0.75 ,pt.cex=1.5)
dev.off()


##########

DF3m1.hist <- calibrate_hists("3m1_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS",".",my.hist.3m2$V5,"DF_3m2") # "25 % downsampling, index : 5"
DF3m2.hist <- calibrate_hists("3m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS",".",my.hist.3m2$V5,"DF_3m2") # "90 % downsampling, index : 18"  "100 % downsampling, index : 20"
DF29m2.hist <- calibrate_hists("29m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READ",".",my.hist.3m2$V5,"DF_3m2") # "20 % downsampling, index : 4"
DF29m6.hist <- calibrate_hists("29m6_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS",".",my.hist.3m2$V5,"DF_3m2") # "35 % downsampling, index : 7"
DF29m7.hist <- calibrate_hists("29m7_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READ",".",my.hist.3m2$V5,"DF_3m2") # "45 % downsampling, index : 9"

### graphical determination!!!!!! (KS not working here), and that was what was done in BD paper
pdf("2016-07-28_Hists_DF_Post_adj_on_mergedCalling.pdf")
plot(0:39,DF3m1.hist[[11]]$V5[1:40],type='l',col="darkblue",ylim=c(0,0.4),xlab="Coverage (X)",ylab="Ehist of H3K4me3 reference peak coverage")
lines(0:39,DF3m2.hist[[20]]$V5[1:40],type='l',col="cadetblue", lwd=2)
lines(0:39,DF29m2.hist[[12]]$V5[1:40],type='l',col="lightgreen")
lines(0:39,DF29m6.hist[[13]]$V5[1:40],type='l',col="darkgreen")
lines(0:39,DF29m7.hist[[15]]$V5[1:40],type='l',col="green")
legend("topright",c("3 mths (rep 1) 55% DS","3 mths (rep 2) 100% DS","29 mths (rep 2) 60% DS","29 mths (rep 6) 65% DS","29 mths (rep 7) 75% DS"),
       col=c("darkblue","cadetblue","lightgreen","darkgreen","green"),bty='n',pch='_',cex = 0.75 ,pt.cex=1.5)
dev.off()

  