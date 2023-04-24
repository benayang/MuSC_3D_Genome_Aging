prefix=/nas/homes/benyang/HiC

outdir=$prefix/08_HiCExplorer/kmeans_TAD_boundaries/track_coverage_plots

############# Plot histone markers for each cluster
for feature in atac H3K4me3 H4K20me1 H3K27me3
do

    # computeMatrix reference-point \
    # -R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster1.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster2.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster4.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster5.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster6.bedgraph \
    # -S $prefix/HistoneTracks/get_signal_bigwigs/$feature.count.rpkm.track.nodup.Aged.bw \
    # --outFileName $outdir/aged.$feature.bin5kb.mat.gz \
    # --referencePoint center -b 250000 -a 250000 \
    # --binSize 5000 \
    # -p 55

    plotHeatmap -m $outdir/aged.$feature.bin5kb.mat.gz \
    --outFileName $outdir/aged.$feature.bin5kb.png \
    --plotType se \
    --perGroup \
    --xAxisLabel "distance (bp)" \
    --samplesLabel $feature \
    --regionsLabel "Cluster 1" "Cluster 2" "Cluster 4" "Cluster 5" "Cluster 6" \
    --dpi 300

    # computeMatrix reference-point \
    # -R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster1.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster2.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster3.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster5.bedgraph \
    # $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster6.bedgraph \
    # -S $prefix/HistoneTracks/get_signal_bigwigs/$feature.count.rpkm.track.nodup.Aged.bw \
    # --outFileName $outdir/young.$feature.bin5kb.mat.gz \
    # --referencePoint center -b 250000 -a 250000 \
    # --binSize 5000 \
    # -p 55

    plotHeatmap -m $outdir/young.$feature.bin5kb.mat.gz \
    --outFileName $outdir/young.$feature.bin5kb.png \
    --plotType se \
    --perGroup \
    --xAxisLabel "distance (bp)" \
    --samplesLabel $feature \
    --regionsLabel "Cluster 1" "Cluster 2" "Cluster 3" "Cluster 5" "Cluster 6" \
    --dpi 300

done
