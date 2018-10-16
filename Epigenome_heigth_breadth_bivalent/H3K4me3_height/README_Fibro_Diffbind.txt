
# How to annotate the Diffbind table with peaks?
cat DF_Diffbind_matrix_concentration.txt | cut -f 1,2,3 > DF_Diffbind_peaks_coord.bed
annotatePeaks.pl DF_Diffbind_peaks_coord.bed mm9 > HOMER_DF_Diffbind_peaks_coord.xls

Will do the DESeq2 linear modeling?
