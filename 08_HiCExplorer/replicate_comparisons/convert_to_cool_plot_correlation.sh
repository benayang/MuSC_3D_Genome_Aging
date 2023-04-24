#hicConvertFormat -m 4504-KS-1.hic \
#-o 4504-KS-1_250kb.cool \
#--inputFormat hic \
#--outputFormat cool \
#--resolutions 250000

hicCorrelate -m ./4140-KS-2_250kb_250000.cool \
./4140-KS-3_250kb_250000.cool \
./4140-KS-1_250kb_250000.cool \
./4504-KS-1_250kb_250000.cool \
--threads 50 \
--log1p \
--method=pearson \
--labels aged.r1 aged.r2 young.r1 young.r2 \
--outFileNameHeatmap young_aged_heatmap_KR_pearson \
--plotFileFormat png \
--outFileNameScatter young_aged_scatterplot_KR_pearson \
--plotFileFormat png
