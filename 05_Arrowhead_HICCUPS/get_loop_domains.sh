prefix="/nas/homes/benyang/HiC"

awk '{OFS="\t"} {if(NR>2) print}' $prefix/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n | uniq > $prefix/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_domains.bed

awk '{OFS="\t"} {if(NR>2) print}' $prefix/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops.bedpe | cut -f 1,2,6 | sort -k1,1 -k2,2n | uniq > $prefix/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_domains.bed