prefix=/nas/homes/benyang/HiC

cut -f1,2,6 $prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_differential_loops2.bedpe | tail -n +3 | uniq > $prefix/05_Arrowhead_HICCUPS/hiccups_diff/aged_diff_loop_domains.bed
cut -f1,2,6 $prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_differential_loops1.bedpe | tail -n +3 | uniq > $prefix/05_Arrowhead_HICCUPS/hiccups_diff/young_diff_loop_domains.bed