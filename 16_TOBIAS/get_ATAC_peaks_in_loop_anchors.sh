prefix='/nas/homes/benyang/HiC'

# get merged peakset across ATAC peaks
cut -f1-3 $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak | awk '{printf "%s\t%s\n", $0,"Aged"}' | \
uniq | sort -k1,1 -k2,2n > $prefix/16_TOBIAS/atac.d0.Aged.bed

cut -f1-3 $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak | awk '{printf "%s\t%s\n", $0,"Young"}' | \
uniq | sort -k1,1 -k2,2n > $prefix/16_TOBIAS/atac.d0.Young.bed

cat $prefix/16_TOBIAS/atac.d0.Aged.bed $prefix/16_TOBIAS/atac.d0.Young.bed | \
uniq | sort -k1,1 -k2,2n | \
bedtools merge -i stdin -c 4 -o distinct > $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed

# get individual ATAC peaks in loop anchors
bedtools intersect -wa \
-a $prefix/16_TOBIAS/atac.d0.Aged.bed \
-b $prefix/05_Arrowhead_HICCUPS/aged_merged_loop_anchors.bed | \
uniq | sort -k1,1 -k2,2n > $prefix/16_TOBIAS/aged_ATAC_peaks_in_loop_anchors.bed

bedtools intersect -wa \
-a $prefix/16_TOBIAS/atac.d0.Young.bed \
-b $prefix/05_Arrowhead_HICCUPS/young_merged_loop_anchors.bed | \
uniq | sort -k1,1 -k2,2n > $prefix/16_TOBIAS/young_ATAC_peaks_in_loop_anchors.bed

# get merged loop anchor lists
cat $prefix/05_Arrowhead_HICCUPS/aged_merged_loop_anchors.bed \
$prefix/05_Arrowhead_HICCUPS/young_merged_loop_anchors.bed | \
uniq | sort -k1,1 -k2,2n | \
bedtools merge -i stdin > $prefix/16_TOBIAS/young_aged_loop_anchors_merged.bed

bedtools intersect -wa \
-a $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
-b $prefix/16_TOBIAS/young_aged_loop_anchors_merged.bed | \
uniq | sort -k1,1 -k2,2n > $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed

# get merged non diff loop anchor lists
cat $prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_non_diff_loop_anchors.bed \
$prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_non_diff_loop_anchors.bed | \
uniq | sort -k1,1 -k2,2n | \
bedtools merge -i stdin -c 4 -o distinct > $prefix/16_TOBIAS/young_aged_non_diff_loop_anchors_merged.bed

bedtools intersect -wa \
-a $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
-b $prefix/16_TOBIAS/young_aged_non_diff_loop_anchors_merged.bed > \
$prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_non_diff_loop_anchors.bed

bedtools intersect -wa \
-a $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
-b $prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_diff_loop_anchors.bed > \
$prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_aged_diff_loop_anchors.bed

bedtools intersect -wa \
-a $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
-b $prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_diff_loop_anchors.bed > \
$prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_young_diff_loop_anchors.bed