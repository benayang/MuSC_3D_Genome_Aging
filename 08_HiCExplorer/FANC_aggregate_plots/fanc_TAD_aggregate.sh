prefix=/nas/homes/benyang/HiC

#### Plot domains
# fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.40kb.hic@40kb@KR \
# $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_domains_40kb.agg \
# -p $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_domains_40kb_oe_large.agg.png \
# -m $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_domains_40kb_oe_large.agg.txt \
# -e -l -r 1.0

# fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@40kb@KR \
# $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_domains_40kb.agg \
# -p $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_domains_40kb_oe_large.agg.png \
# -m $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_domains_40kb_oe_large.agg.txt \
# -e -l -r 1.0

#### Plot Boundaries
# fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.40kb.hic@40kb@KR \
# $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_boundaries_oe_40kb.agg \
# -w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_boundaries_40kb_oe_large.agg.png \
# -m $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_TAD_boundaries_40kb_oe_large.agg.txt \
# -e -l

# fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@40kb@KR \
# $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_boundaries_oe_40kb.agg \
# -w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_boundaries_40kb_oe_large.agg.png \
# -m $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_TAD_boundaries_40kb_oe_large.agg.txt \
# -e -l

#### Plot gained/lost Boundaries
fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.40kb.hic@40kb@KR \
$prefix/08_HiCExplorer/unique_aged_TAD_boundary.bed \
$prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_boundaries_oe_40kb.agg \
-w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_boundaries_40kb_oe_large.agg.png \
-m $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_boundaries_40kb_oe_large.agg.txt \
-e -l

fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@40kb@KR \
$prefix/08_HiCExplorer/unique_young_TAD_boundary.bed \
$prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_boundaries_oe_40kb.agg \
-w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_boundaries_40kb_oe_large.agg.png \
-m $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_boundaries_40kb_oe_large.agg.txt \
-e -l

fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.40kb.hic@40kb@KR \
$prefix/08_HiCExplorer/shared_TAD_boundary.bed \
$prefix/08_HiCExplorer/FANC_aggregate_plots/aged_shared_TAD_boundaries_oe_40kb.agg \
-w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_shared_TAD_boundaries_40kb_oe_large.agg.png \
-m $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_shared_TAD_boundaries_40kb_oe_large.agg.txt \
-e -l

fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@40kb@KR \
$prefix/08_HiCExplorer/shared_TAD_boundary.bed \
$prefix/08_HiCExplorer/FANC_aggregate_plots/young_shared_TAD_boundaries_oe_40kb.agg \
-w 5mb -p $prefix/08_HiCExplorer/FANC_aggregate_plots/young_shared_TAD_boundaries_40kb_oe_large.agg.png \
-m $prefix/08_HiCExplorer/FANC_aggregate_plots/young_shared_TAD_boundaries_40kb_oe_large.agg.txt \
-e -l

#### Plot loops
# fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.hic@50kb@KR \
# $prefix/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_coordOnly.bedpe \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_loops_5kb.agg \
# -p $prefix/08_HiCExplorer/FANC_aggregate_plots/aged_merged_loops_5kb_oe_large.agg.png \
# --loops

# fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@5kb@KR \
# $prefix/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_coordOnly.bedpe \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_loops_5kb.agg \
# -p $prefix/08_HiCExplorer/FANC_aggregate_plots/young_merged_loops_5kb_oe_large.agg.png \
# --loops