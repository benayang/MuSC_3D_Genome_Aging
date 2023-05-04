prefix='/nas/homes/benyang/HiC'

# TOBIAS BINDetect \
# --motifs $prefix/16_TOBIAS/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
# --signals $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
# --peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
# --output-peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_young_diff_loop_anchors.bed \
# --genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
# --cond-names young_merged_ATAC aged_merged_ATAC \
# --outdir $prefix/16_TOBIAS/young_diff_loops \
# --prefix bindetect_young_diff_loop_anchors \
# --cores 50

# TOBIAS BINDetect \
# --motifs $prefix/16_TOBIAS/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
# --signals $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
# --peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
# --output-peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_aged_diff_loop_anchors.bed \
# --genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
# --cond-names young_merged_ATAC aged_merged_ATAC \
# --outdir $prefix/16_TOBIAS/aged_diff_loops \
# --prefix bindetect_aged_diff_loop_anchors \
# --cores 50

# TOBIAS BINDetect \
# --motifs $prefix/16_TOBIAS/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
# --signals $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
# --peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
# --genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
# --cond-names young_merged_ATAC aged_merged_ATAC \
# --outdir $prefix/16_TOBIAS/genomewide \
# --prefix bindetect_genomewide \
# --cores 50

# TOBIAS BINDetect \
# --motifs $prefix/16_TOBIAS/JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
# --signals $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
# --peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_annotated.bed \
# --peak-header $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_annotated_header.txt \
# --output-peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
# --genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
# --cond-names young_merged_ATAC aged_merged_ATAC \
# --outdir $prefix/16_TOBIAS/merged_peakset_annotated \
# --prefix bindetect_merged_peakset \
# --cores 50

TOBIAS BINDetect \
--motifs $prefix/16_TOBIAS/expressed_JASPAR2022_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
--signals $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
--peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_annotated_parsed.bed \
--peak-header $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_annotated_header.txt \
--output-peaks $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
--genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
--cond-names young_merged_ATAC aged_merged_ATAC \
--outdir $prefix/16_TOBIAS/merged_peakset_annotated_expressed_TF \
--prefix bindetect_merged_peakset_expressed \
--cores 50
