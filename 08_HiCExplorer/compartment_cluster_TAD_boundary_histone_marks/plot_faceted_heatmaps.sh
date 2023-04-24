prefix=/nas/homes/benyang/HiC

for age in young aged
do
    #for feature in atac H3K4me3 H4K20me1 H3K27me3

    ##### A/B 
    # plotHeatmap -m $prefix/08_HiCExplorer/TAD_histone_marks/${age}_AB.kmeans4.mat.gz \
    # --whatToShow "heatmap and colorbar" \
    # --outFileName ${age}.AB.bin10kb.kmeans4.hmpOnly.png \
    # --zMin -0.05 --zMax 0.05 \
    # --heatmapHeight 10 \
    # --heatmapWidth 3.5 \
    # --xAxisLabel "Boundary Distance (bp)" \
    # --samplesLabel "A/B" \
    # --regionsLabel "1" "2" "3" "4" \
    # --colorMap RdYlBu_r RdYlBu \
    # --dpi 300

    ##### H4K20me1 
    # plotHeatmap -m ../${age}.H4K20me1.bin10kb.mat.gz \
    # --whatToShow "heatmap and colorbar" \
    # --outFileName ${age}.H4K20me1.bin10kb.kmeans4.hmpOnly.png \
    # --zMin 0 --zMax 800\
    # --heatmapHeight 10 \
    # --heatmapWidth 3.5 \
    # --xAxisLabel "Boundary Distance (bp)" \
    # --samplesLabel "H4K20me1" \
    # --regionsLabel "1" "2" "3" "4" \
    # --colorMap RdYlBu_r RdYlBu \
    # --dpi 300

    ##### ATAC 
    plotHeatmap -m ../${age}.atac.bin10kb.mat.gz \
    --whatToShow "heatmap and colorbar" \
    --outFileName ${age}.ATAC.bin10kb.kmeans4.hmpOnly.png \
    --zMin 0 --zMax 80 \
    --heatmapHeight 10 \
    --heatmapWidth 3.5 \
    --xAxisLabel "Boundary Distance (bp)" \
    --samplesLabel "ATAC" \
    --regionsLabel "1" "2" "3" "4" \
    --colorMap RdYlBu_r RdYlBu \
    --dpi 300

    ##### H3K4me3 
    plotHeatmap -m ../${age}.H3K4me3.bin10kb.mat.gz \
    --whatToShow "heatmap and colorbar" \
    --outFileName ${age}.H3K4me3.bin10kb.kmeans4.hmpOnly.png \
    --zMin 0 --zMax 80 \
    --heatmapHeight 10 \
    --heatmapWidth 3.5 \
    --xAxisLabel "Boundary Distance (bp)" \
    --samplesLabel "H3K4me3" \
    --regionsLabel "1" "2" "3" "4" \
    --colorMap RdYlBu_r RdYlBu \
    --dpi 300

    ##### H3K27me3 
    plotHeatmap -m ../${age}.H3K27me3.bin10kb.mat.gz \
    --whatToShow "heatmap and colorbar" \
    --outFileName ${age}.H3K27me3.bin10kb.kmeans4.hmpOnly.png \
    --zMin 0 --zMax 40 \
    --heatmapHeight 10 \
    --heatmapWidth 3.5 \
    --xAxisLabel "Boundary Distance (bp)" \
    --samplesLabel "H3K27me3" \
    --regionsLabel "1" "2" "3" "4" \
    --colorMap RdYlBu_r RdYlBu \
    --dpi 300

done