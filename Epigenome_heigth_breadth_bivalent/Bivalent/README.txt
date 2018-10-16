# 2016-06-01

do bivalent analysis using meta peaks for analysis (overlap of changed peaks with age)

ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed
ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K27me3_height_fibro/ALL_AGES_MERGED_DermalFibro_H3K27me3_broad_peaks.bed 


intersectBed -wa -a ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed -b ALL_AGES_MERGED_DermalFibro_H3K27me3_broad_peaks.bed > ALL_AGES_MERGED_DermalFibro_Bivalent_domains.bed

annotatePeaks.pl ALL_AGES_MERGED_DermalFibro_Bivalent_domains.bed mm9 > HOMER_ALL_AGES_MERGED_DermalFibro_Bivalent_domains.xls
