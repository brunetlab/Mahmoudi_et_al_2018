# 2016-05-19

# run for Salah, potentially included in his paper

cat /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m2_ear_fibroblasts_H3K27me3_1_trimmed.FIXSEQ_CLEANED_READS.bed /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m1_ear_fibroblasts_H3K27me3_1_trimmed.FIXSEQ_CLEANED_READS.bed \
 /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m7_ear_fibroblasts_H3K27me3_1_trimmed.FIXSEQ_CLEANED_READS.bed /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m6_ear_fibroblasts_H3K27me3_1_trimmed.FIXSEQ_CLEANED_READS.bed /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m2_ear_fibroblasts_H3K27me3_1_trimmed.FIXSEQ_CLEANED_READS.bed  \
	> ALL_AGES_MERGED_DermalFibro_H3K27me3.FIXSEQ_CLEANED_READS.bed


macs2 -t ALL_AGES_MERGED_DermalFibro_H3K27me3.FIXSEQ_CLEANED_READS.bed -c ../H3K4me3_breadth/ALL_AGES_MERGED_DermalFibro_INPUT.FIXSEQ_CLEANED_READS.bed -n ALL_AGES_MERGED_DermalFibro_H3K27me3 -f 'BED' -g 'mm' --broad --keep-dup=all

# 2016-05-23

cat DF_H3K27me3_Diffbind_matrix_concentration.txt | cut -f 1,2,3 > DF_H3K27me3_Diffbind_peaks_coord.bed
annotatePeaks.pl DF_H3K27me3_Diffbind_peaks_coord.bed mm9 > HOMER_DF_H3K27me3_Diffbind_peaks_coord.xls

