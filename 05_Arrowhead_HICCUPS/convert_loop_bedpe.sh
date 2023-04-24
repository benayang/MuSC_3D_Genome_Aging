awk '{OFS="\t"} NR>2 {print $1,$2,$3,$4,$5,$6,$12/$14}' aged_merged_hiccups/merged_loops.bedpe > aged_merged_hiccups/merged_loops_links.bedpe
awk '{OFS="\t"} NR>2 {print $1,$2,$3,$4,$5,$6,$12/$14}' young_merged_hiccups/merged_loops.bedpe > young_merged_hiccups/merged_loops_links.bedpe

awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' young_merged_hiccups/merged_loops.bedpe | tail -n +3 > young_merged_hiccups/merged_loops_coordOnly.bedpe
awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' aged_merged_hiccups/merged_loops.bedpe | tail -n +3 > aged_merged_hiccups/merged_loops_coordOnly.bedpe
