prefix=/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS

tail -n +3 aged_differential_loops2.bedpe | \
cut -f 1-3 >> aged_diff_loop_anchors.bed 

tail -n +3 aged_differential_loops2.bedpe | \
cut -f 4-6 >> aged_diff_loop_anchors.bed

sort -k1,1 -k2,2n aged_diff_loop_anchors.bed > aged_diff_loop_anchors.sorted.bed

tail -n +3 young_differential_loops1.bedpe | \
cut -f 1-3 >> young_diff_loop_anchors.bed 

tail -n +3 young_differential_loops1.bedpe | \
cut -f 4-6 >> young_diff_loop_anchors.bed

sort -k1,1 -k2,2n young_diff_loop_anchors.bed > young_diff_loop_anchors.sorted.bed