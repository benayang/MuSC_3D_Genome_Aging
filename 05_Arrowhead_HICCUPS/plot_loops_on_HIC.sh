prefix=/nas/homes/benyang/HiC

hicPlotMatrix --matrix $prefix/08_HiCExplorer/aged.merged_5kb_KR.cool \
--log1p \
--outFileName aged_merged_5kb_loops_region2.png \
--region chr7:66500000-69000000 \
--loops aged_merged_hiccups/merged_loops.bedpe \
--dpi 300
#--region chr3:78000000-81000000 \

hicPlotMatrix --matrix $prefix/08_HiCExplorer/young.merged_5kb_KR.cool \
--log1p \
--outFileName young_merged_5kb_loops_region2.png \
--region chr7:66500000-69000000 \
--loops young_merged_hiccups/merged_loops.bedpe \
--dpi 300