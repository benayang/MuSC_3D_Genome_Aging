prefix=/nas/homes/benyang/HiC

outdir=$prefix/08_HiCExplorer/kmeans_TAD_boundaries/track_coverage_plots

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster4.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster6.bedgraph \
-S /nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
--outFileName $outdir/aged.ATAC.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster4.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster6.bedgraph \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName $outdir/aged.H3K27me3.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster4.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/aged_cluster6.bedgraph \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName $outdir/aged.H3K4me3.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster3.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster6.bedgraph \
-S /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
--outFileName $outdir/young.ATAC.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster3.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster6.bedgraph \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName $outdir/young.H3K27me3.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster1.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster2.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster3.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster5.bedgraph \
$prefix/08_HiCExplorer/kmeans_TAD_boundaries/young_cluster6.bedgraph \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName $outdir/young.H3K4me3.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 1000 \
-p 50

for f in ATAC H3K27me3 H3K4me3
do

plotHeatmap -m $outdir/aged.$f.foldchange.mat.gz \
--outFileName $outdir/aged.$f.foldchange.histone.png \
--perGroup \
--plotType se \
--xAxisLabel "distance (bp)" \
--samplesLabel $f \
--dpi 300

plotHeatmap -m $outdir/young.$f.foldchange.mat.gz \
--outFileName $outdir/young.$f.foldchange.histone.png \
--perGroup \
--plotType se \
--xAxisLabel "distance (bp)" \
--samplesLabel $f \
--dpi 300

done