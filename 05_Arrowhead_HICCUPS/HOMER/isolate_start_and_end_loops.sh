prefix=/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS

cut -f 1-3 $prefix/aged_merged_hiccups/aged_merged_loops_noHeader.bedpe | awk '{OFS="\t"}{print $s,NR}' > aged_merged_loops_start.bed
cut -f 4-6 $prefix/aged_merged_hiccups/aged_merged_loops_noHeader.bedpe | awk '{OFS="\t"}{print $s,NR}' > aged_merged_loops_end.bed

cut -f 1-3 $prefix/young_merged_hiccups/young_merged_loops_noHeader.bedpe | awk '{OFS="\t"}{print $s,NR}' > young_merged_loops_start.bed
cut -f 4-6 $prefix/young_merged_hiccups/young_merged_loops_noHeader.bedpe | awk '{OFS="\t"}{print $s,NR}' > young_merged_loops_end.bed

cat aged_merged_loops_start.bed aged_merged_loops_end.bed | cut -f 1-3 | sort -k1,1 -k2,2n | uniq > aged_merged_loop_anchors.bed
cat young_merged_loops_start.bed young_merged_loops_end.bed | cut -f 1-3 | sort -k1,1 -k2,2n | uniq > young_merged_loop_anchors.bed



tail -n +2 $prefix/hiccups_diff/aged_differential_loops2.bedpe | cut -f 1-3 | awk '{OFS="\t"}{print $s,NR}' > aged_merged_diff_loops_start.bed
tail -n +2 $prefix/hiccups_diff/aged_differential_loops2.bedpe | cut -f 4-6 | awk '{OFS="\t"}{print $s,NR}' > aged_merged_diff_loops_end.bed

tail -n +2 $prefix/hiccups_diff/young_differential_loops1.bedpe | cut -f 1-3 | awk '{OFS="\t"}{print $s,NR}' > young_merged_diff_loops_start.bed
tail -n +2 $prefix/hiccups_diff/young_differential_loops1.bedpe | cut -f 4-6 | awk '{OFS="\t"}{print $s,NR}' > young_merged_diff_loops_end.bed

cat aged_merged_diff_loops_start.bed aged_merged_diff_loops_end.bed | cut -f 1-3 | sort -k1,1 -k2,2n | uniq > aged_merged_diff_loop_anchors.bed
cat young_merged_diff_loops_start.bed young_merged_diff_loops_end.bed | cut -f 1-3 | sort -k1,1 -k2,2n | uniq > young_merged_diff_loop_anchors.bed