#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed 
#
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m1_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m7_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m6_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m1_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/3m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m6_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/BAM_clean/29m7_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS.bed 
#
#breadth_ds_coverage_BED.sh 3m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed
#breadth_ds_coverage_BED.sh 3m1_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed
#breadth_ds_coverage_BED.sh 29m7_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed
#breadth_ds_coverage_BED.sh 29m6_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed
#breadth_ds_coverage_BED.sh 29m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS.bed ALL_AGES_MERGED_DermalFibro_H3K4me3_peaks.bed
#
#
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/29m7_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/29m6_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/29m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/3m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed
#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/3m1_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed

#ln -s /Volumes/MyBook_5/Salah_Fibro_data/ChIP-seq/H3K4me3_breadth/ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed 

#macs2 -t 29m7_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS75percent.bed -c 29m7_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed -n 29m7_ear_fibroblasts_H3K4me3_ADJ -f 'BED' -g 'mm' --broad --keep-dup=all
#macs2 -t 29m6_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS65percent.bed -c 29m6_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed -n 29m6_ear_fibroblasts_H3K4me3_ADJ -f 'BED' -g 'mm' --broad --keep-dup=all
#macs2 -t 29m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS60percent.bed -c 29m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed -n 29m2_ear_fibroblasts_H3K4me3_ADJ -f 'BED' -g 'mm' --broad --keep-dup=all
#macs2 -t 3m2_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS100percent.bed -c 3m2_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed -n 3m2_ear_fibroblasts_H3K4me3_ADJ -f 'BED' -g 'mm' --broad --keep-dup=all
#macs2 -t 3m1_ear_fibroblasts_H3K4me3_combined_trimmed.FIXSEQ_CLEANED_READS55percent.bed -c 3m1_ear_fibroblasts_INPUT_combined_trimmed.FIXSEQ_CLEANED_READS_DS_15203520.bed -n 3m1_ear_fibroblasts_H3K4me3_ADJ -f 'BED' -g 'mm' --broad --keep-dup=all
#
#get_broad_from_bed.pl 1 ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed 3m1_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 3m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 29m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 29m6_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 29m7_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed
#
#ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed 95th percentile: 6228
#3m1_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 95th percentile: 5495
#3m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 95th percentile: 5395
#29m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 95th percentile: 5637
#29m6_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 95th percentile: 5963
#29m7_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed 95th percentile: 6007
#

### IDR like calls for young and old broad peaks
coverageBed -a broad_domains_29m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed -b broad_domains_ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed \
	| awk -v min=0.5 '{if ($NF > min) print $0}' \
	| coverageBed -a broad_domains_29m6_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed -b - \
	| awk -v min=0.5 '{if ($NF > min) print $0}' \
	| coverageBed -a broad_domains_29m7_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed -b - \
	| awk -v min=0.5 '{if ($NF > min) print $0}' \
	| intersectBed -wa -u -f 1 -a broad_domains_ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed -b - | cut -f 1-4 \
	 > broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed

coverageBed -a broad_domains_3m2_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed -b broad_domains_ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed \
	| awk -v min=0.5 '{if ($NF > min) print $0}' \
	| coverageBed -a broad_domains_3m1_ear_fibroblasts_H3K4me3_ADJ_broad_peaks.bed -b - \
	| awk -v min=0.5 '{if ($NF > min) print $0}' \
	| intersectBed -wa -u -f 1 -a broad_domains_ALL_AGES_MERGED_DermalFibro_H3K4me3_broad_peaks.bed -b - | cut -f 1-4 \
	 > broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed

wc -l broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed # 591
wc -l broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed # 614
intersectBed -a broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed -b broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed | wc -l
#473

intersectBed -a broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed -b broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed \
	> Unchanged_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed

intersectBed -v -a broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed -b broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed \
	> 29mONLY_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed
	
intersectBed -v -b broad_domains_29m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed -a broad_domains_3m_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed \
	> 3mONLY_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed
	
118 29mONLY_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed
473 Unchanged_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed
141 3mONLY_broad_domains_ear_fibroblasts_H3K4me3_IDR_broad_peaks.bed
