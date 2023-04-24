prefix=/nas/homes/benyang/HiC

hicInterIntraTAD -m $prefix/08_HiCExplorer/hicDifferentialTAD/aged.merged_40000_minNorm_KR.cool \
--tadDomains $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileName aged.merged_hicInterIntraTAD.txt \
--outFileNameRatioPlot aged.merged_hicInterIntraTAD.png \
--dpi 300 \
--threads 55

hicInterIntraTAD -m $prefix/08_HiCExplorer/hicDifferentialTAD/young.merged_40000_minNorm_KR.cool \
--tadDomains $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
--outFileName young.merged_hicInterIntraTAD.txt \
--outFileNameRatioPlot young.merged_hicInterIntraTAD.png \
--dpi 300 \
--threads 55

# hicInterIntraTAD -m $prefix/08_HiCExplorer/aged.merged_40kb_KR.cool \
# --tadDomains $prefix/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_domains.bed \
# --outFileName aged.merged_loops_hicInterIntraTAD.txt \
# --outFileNameRatioPlot aged.merged_loops_hicInterIntraTAD.png \
# --dpi 300 \
# --threads 55

# hicInterIntraTAD -m $prefix/08_HiCExplorer/young.merged_40kb_KR.cool \
# --tadDomains $prefix/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_domains.bed \
# --outFileName young.merged_loops_hicInterIntraTAD.txt \
# --outFileNameRatioPlot young.merged_loops_hicInterIntraTAD.png \
# --dpi 300 \
# --threads 55