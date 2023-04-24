prefix=/nas/homes/benyang/HiC

# awk '{OFS="\t"}{print $1,$2,$3,$5}' "$prefix/04_FANC/without KR normalization/young.merged_100kb_noKR.ev" > young.merged_100kb_noKR.bedgraph
# awk '{OFS="\t"}{print $1,$2,$3,$5}' "$prefix/04_FANC/without KR normalization/aged.merged_100kb_noKR.ev" > aged.merged_100kb_noKR.bedgraph

#hicCompartmentalization -m $prefix/08_HiCExplorer/aged.merged_100kb_KR_obsexp.cool \
#--pca ./aged.merged_100kb_noKR.filtered.bedgraph \
hicCompartmentalization -m $prefix/08_HiCExplorer/aged.merged_100kb_KR_obsexp.cool \
--pca $prefix/08_HiCExplorer/hicCompartmentalization/hicPCA/aged_pca1.filtered.bedgraph \
-o aged_merged_100kb_KR_filtered_compartmentalization_noOffset_hicExplorerPCA.png \
--quantile 50 \
--outliers 2.5 \
--outputMatrix aged_merged_100kb_KR_filtered_compartmentalization_noOffset_hicExplorerPCA.npz

# hicCompartmentalization -m $prefix/08_HiCExplorer/young.merged_100kb_KR_obsexp.cool \
# --pca ./young.merged_100kb_noKR.filtered.bedgraph \
hicCompartmentalization -m $prefix/08_HiCExplorer/young.merged_100kb_KR_obsexp.cool \
--pca $prefix/08_HiCExplorer/hicCompartmentalization/hicPCA/young_pca1.filtered.bedgraph \
-o young_merged_100kb_KR_filtered_compartmentalization_noOffset_hicExplorerPCA.png \
--quantile 50 \
--outliers 2.5 \
--outputMatrix young_merged_100kb_KR_filtered_compartmentalization_noOffset_hicExplorerPCA.npz