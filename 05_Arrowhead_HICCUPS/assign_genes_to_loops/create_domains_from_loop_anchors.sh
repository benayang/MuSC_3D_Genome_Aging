mainDir="/nas/homes/benyang/HiC"

prefix=/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS

cut -f 1,2,6 $prefix/aged_merged_hiccups/merged_loops_noHeader.bedpe | sort -k1,1 -k2,2n | uniq > aged_merged_loop_domains.bed
cut -f 1,2,6 $prefix/young_merged_hiccups/merged_loops_noHeader.bedpe | sort -k1,1 -k2,2n | uniq > young_merged_loop_domains.bed

tail -n +3 $prefix/hiccups_diff/young_differential_loops1.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n | uniq > young_diff_loop_domains.bed
tail -n +3 $prefix/hiccups_diff/aged_differential_loops2.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n | uniq > aged_diff_loop_domains.bed 