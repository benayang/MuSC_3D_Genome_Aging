prefix=/nas/homes/benyang/HiC

#### Motifs in all loop anchors
#findMotifsGenome.pl $prefix/05_Arrowhead_HICCUPS/HOMER/young_merged_loop_anchors.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_merged_loop_anchors -p 55
#findMotifsGenome.pl $prefix/05_Arrowhead_HICCUPS/HOMER/aged_merged_loop_anchors.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_merged_loop_anchors -p 55

#### Motifs in all chanors of loops within TADs
# bedtools intersect -wa -a $prefix/05_Arrowhead_HICCUPS/HOMER/young_merged_loop_anchors.bed \
# -b $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
# sort -k1,1 -k2,2n | uniq > young_merged_loop_anchors_in_TAD_domain.bed

# bedtools intersect -wa -a $prefix/05_Arrowhead_HICCUPS/HOMER/aged_merged_loop_anchors.bed \
# -b $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed | \
# sort -k1,1 -k2,2n | uniq > aged_merged_loop_anchors_in_TAD_domain.bed

# findMotifsGenome.pl young_merged_loop_anchors_in_TAD_domain.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_merged_loop_anchors_in_TADs -p 55
# findMotifsGenome.pl aged_merged_loop_anchors_in_TAD_domain.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_merged_loop_anchors_in_TADs -p 55

#### Motifs in A compartment
# bedtools intersect -wa -a $prefix/05_Arrowhead_HICCUPS/HOMER/young_merged_loop_anchors.bed \
# -b $prefix/04_FANC/compartmentExpression/compartmentBed/young.A.bed | \
# sort -k1,1 -k2,2n | uniq > young_merged_loop_anchors_in_A_compartment.bed

# bedtools intersect -wa -a $prefix/05_Arrowhead_HICCUPS/HOMER/aged_merged_loop_anchors.bed \
# -b $prefix/04_FANC/compartmentExpression/compartmentBed/aged.A.bed | \
# sort -k1,1 -k2,2n | uniq > aged_merged_loop_anchors_in_A_compartment.bed

# findMotifsGenome.pl young_merged_loop_anchors_in_A_compartment.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_merged_loop_anchors_in_AComp -p 55
# findMotifsGenome.pl aged_merged_loop_anchors_in_A_compartment.bed mm10 $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_merged_loop_anchors_in_AComp -p 55

#### Motifs in ATAC sites within loop anchors across all ATAC sites
# findMotifsGenome.pl \
# $prefix/16_TOBIAS/young_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_ATAC_loop_anchors \
# -bg $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
# -p 55

# findMotifsGenome.pl \
# $prefix/16_TOBIAS/aged_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_ATAC_loop_anchors \
# -bg $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
# -p 55

#### Motifs in ATAC sites within loop anchors against genome background
# findMotifsGenome.pl \
# $prefix/16_TOBIAS/young_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_ATAC_loop_anchors_genome_bg \
# -p 55

# findMotifsGenome.pl \
# $prefix/16_TOBIAS/aged_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_ATAC_loop_anchors_genome_bg \
# -p 55

#### Motifs in ATAC sites within loop anchors against age-specific ATAC background
# findMotifsGenome.pl \
# $prefix/16_TOBIAS/young_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_ATAC_loop_anchors_ATAC_bg \
# -bg $prefix/16_TOBIAS/atac.d0.Young.bed \
# -p 55

# findMotifsGenome.pl \
# $prefix/16_TOBIAS/aged_ATAC_peaks_in_loop_anchors.bed \
# mm10 \
# $prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_ATAC_loop_anchors_ATAC_bg \
# -bg $prefix/16_TOBIAS/atac.d0.Aged.bed \
# -p 55

#### Motifs in ATAC sites within loop anchors against all loop anchor ATAC sites
findMotifsGenome.pl \
$prefix/16_TOBIAS/young_ATAC_peaks_in_loop_anchors.bed \
mm10 \
$prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/young_ATAC_loop_anchors_ATAC_loops_bg \
-bg $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
-p 30

findMotifsGenome.pl \
$prefix/16_TOBIAS/aged_ATAC_peaks_in_loop_anchors.bed \
mm10 \
$prefix/05_Arrowhead_HICCUPS/HOMER/call_homer/aged_ATAC_loop_anchors_ATAC_loops_bg \
-bg $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks_merged_loop_anchors.bed \
-p 30
