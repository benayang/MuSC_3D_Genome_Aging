prefix=/nas/homes/benyang/HiC

grep A $prefix/04_FANC/without_KR_normalization/aged.ab_100kb.bed | \
bedtools intersect -wa -u -a $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b stdin | sort -k1,1 -k2,2n > $prefix/08_HiCExplorer/TAD_boundary_compartment_overlap/aged_TADboundary_AComp_overlap.bed

grep B $prefix/04_FANC/without_KR_normalization/aged.ab_100kb.bed | \
bedtools intersect -wa -u -a $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b stdin | sort -k1,1 -k2,2n > $prefix/08_HiCExplorer/TAD_boundary_compartment_overlap/aged_TADboundary_BComp_overlap.bed

grep A $prefix/04_FANC/without_KR_normalization/young.ab_100kb.bed | \
bedtools intersect -wa -u -a $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b stdin | sort -k1,1 -k2,2n > $prefix/08_HiCExplorer/TAD_boundary_compartment_overlap/young_TADboundary_AComp_overlap.bed

grep B $prefix/04_FANC/without_KR_normalization/young.ab_100kb.bed | \
bedtools intersect -wa -u -a $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
-b stdin | sort -k1,1 -k2,2n > $prefix/08_HiCExplorer/TAD_boundary_compartment_overlap/young_TADboundary_BComp_overlap.bed