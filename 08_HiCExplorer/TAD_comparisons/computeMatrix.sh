# computeMatrix reference-point \
# -S "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# --samplesLabel Young Aged \
# -R aggregate_TAD_boundaries.bed \
# --outFileName aggregate_TAD_boundaries_scores.gz \
# --outFileNameMatrix aggregate_TAD_boundaries_scores_mat.txt \
# --outFileSortedRegions aggregate_TAD_boundaries_sortedRegions.bed \
# --referencePoint center \
# --beforeRegionStartLength 200000 \
# --afterRegionStartLength 200000 \
# --binSize 40000 \
# --sortRegions keep \
# -bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
# -p 45

# plotHeatmap \
# -m aggregate_TAD_boundaries_scores.gz \
# --plotType se \
# --dpi 300 \
# --kmeans 3 \
# --silhouette \
# --outFileName aggregate_TAD_boundaries_scores.png

# computeMatrix reference-point \
# -S "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# --samplesLabel Young Aged \
# -R "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
# --outFileName young_TAD_boundaries_scores.gz \
# --outFileNameMatrix young_TAD_boundaries_scores_mat.txt \
# --outFileSortedRegions young_TAD_boundaries_sortedRegions.bed \
# --referencePoint center \
# --beforeRegionStartLength 400000 \
# --afterRegionStartLength 400000 \
# --binSize 40000 \
# --sortRegions keep \
# -bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
# -p 45

# plotHeatmap \
# -m young_TAD_boundaries_scores.gz \
# --plotType se \
# --dpi 300 \
# --kmeans 3 \
# --silhouette \
# --outFileName young_TAD_boundaries_scores.png

# computeMatrix reference-point \
# -S "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# --samplesLabel Young Aged \
# -R "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
# --outFileName aged_TAD_boundaries_scores.gz \
# --outFileNameMatrix aged_TAD_boundaries_scores_mat.txt \
# --outFileSortedRegions aged_TAD_boundaries_sortedRegions.bed \
# --referencePoint center \
# --beforeRegionStartLength 400000 \
# --afterRegionStartLength 400000 \
# --binSize 40000 \
# --sortRegions keep \
# -bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
# -p 45

# plotHeatmap \
# -m aged_TAD_boundaries_scores.gz \
# --plotType se \
# --dpi 300 \
# --kmeans 3 \
# --silhouette \
# --outFileName aged_TAD_boundaries_scores.png

# computeMatrix reference-point \
# -S "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bw" \
# --samplesLabel Young Aged \
# -R "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
# "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
# "aggregate_TAD_boundaries.bed" \
# --outFileName grouped_TAD_boundaries_scores.gz \
# --outFileNameMatrix grouped_TAD_boundaries_scores_mat.txt \
# --outFileSortedRegions grouped_TAD_boundaries_sortedRegions.bed \
# --referencePoint center \
# --beforeRegionStartLength 400000 \
# --afterRegionStartLength 400000 \
# --binSize 40000 \
# --sortRegions keep \
# -bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
# -p 45

# plotHeatmap \
# -m grouped_TAD_boundaries_scores.gz \
# --plotType se \
# --dpi 300 \
# --outFileName grouped_TAD_boundaries_scores.png

computeMatrix reference-point \
-S aged_vs_young_TAD_score_log2ratio.bw \
--samplesLabel Aged_vs_Young \
-R "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
"/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed" \
"aggregate_TAD_boundaries.bed" \
--outFileName log2_grouped_TAD_boundaries_scores.gz \
--outFileNameMatrix log2_grouped_TAD_boundaries_scores_mat.txt \
--outFileSortedRegions log2_grouped_TAD_boundaries_sortedRegions.bed \
--referencePoint center \
--beforeRegionStartLength 400000 \
--afterRegionStartLength 400000 \
--binSize 40000 \
--sortRegions keep \
-bl "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
-p 45

plotHeatmap \
-m log2_grouped_TAD_boundaries_scores.gz \
--plotType se \
--dpi 300 \
--outFileName log2_grouped_TAD_boundaries_scores.png