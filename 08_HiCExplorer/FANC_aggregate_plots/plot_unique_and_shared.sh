prefix=/nas/homes/benyang/HiC

# fanc aggregate $prefix/02_HIC/aged.merged/aged.merged.40kb.hic@40kb@KR \
# $prefix/08_HiCExplorer/TAD\ expression/unique_and_shared_genes/unique_aged_TAD_domain.bed \
# $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_domains_40kb.agg \
# -p $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_domains_40kb_oe_large.agg.png \
# -m $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_aged_TAD_domains_40kb_oe_large.agg.txt \
# -e -l -r 1.0

fanc aggregate $prefix/02_HIC/young.merged/young.merged.hic@40kb@KR \
$prefix/08_HiCExplorer/TAD\ expression/unique_and_shared_genes/unique_young_TAD_domain.bed \
$prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_domains_40kb.agg \
-p $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_domains_40kb_oe_large.agg.png \
-m $prefix/08_HiCExplorer/FANC_aggregate_plots/unique_young_TAD_domains_40kb_oe_large.agg.txt \
-e -l -r 1.0