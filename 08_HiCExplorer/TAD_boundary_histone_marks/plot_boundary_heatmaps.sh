prefix=/nas/homes/benyang/HiC
trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

for f in "A" "B" "A_to_B" "B_to_A" "StaticA" "StaticB"
do

########## RPKM heatmaps
# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADBoundary.bed \
# -S $trackPrefix/H4K20me1.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/H3K27me3.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
# --outFileName aged.$f.mat.gz \
# --referencePoint center -b 60000 -a 60000 \
# --binSize 100 \
# -p 50

# plotHeatmap -m aged.$f.mat.gz \
# --outFileName aged.$f.histone.png \
# --xAxisLabel "distance (bp)" \
# --samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
# --dpi 300

# computeMatrix reference-point \
# -R $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADBoundary.bed \
# -S $trackPrefix/H4K20me1.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/atac.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/H3K27me3.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/H3K4me3.count.rpkm.track.nodup.Young.bw \
# --outFileName young.$f.mat.gz \
# --referencePoint center -b 60000 -a 60000 \
# --binSize 100 \
# -p 50

# plotHeatmap -m young.$f.mat.gz \
# --outFileName young.$f.histone.png \
# --xAxisLabel "distance (bp)" \
# --samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
# --dpi 300

################## fold change marks
computeMatrix reference-point \
-R $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADBoundary.bed \
-S /nas/homes/annashch/H4K20me1/caper_out/chip/7aab8e70-c16d-481f-ab6d-ebd0e3683fae/call-macs2_signal_track_pooled/execution/rep.pooled.fc.signal.bigwig \
/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName aged.$f.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 100 \
-p 50

plotHeatmap -m aged.$f.foldchange.mat.gz \
--outFileName aged.$f.foldchange.histone.png \
--xAxisLabel "distance (bp)" \
--samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
--dpi 300

computeMatrix reference-point \
-R $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADBoundary.bed \
-S /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
/nas/homes/annashch/H4K20me1/caper_out/chip/72cab657-4520-4eee-b95a-d87eac1cecbf/call-macs2_signal_track_pooled/execution/rep.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
--outFileName young.$f.foldchange.mat.gz \
--referencePoint center -b 60000 -a 60000 \
--binSize 100 \
-p 50

plotHeatmap -m young.$f.foldchange.mat.gz \
--outFileName young.$f.foldchange.histone.png \
--xAxisLabel "distance (bp)" \
--samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
--dpi 300

done
