prefix=/nas/homes/benyang/HiC
trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

computeMatrix scale-regions \
-R $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Aged_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
/nas/homes/annashch/H4K20me1/caper_out/chip/7aab8e70-c16d-481f-ab6d-ebd0e3683fae/call-macs2_signal_track_pooled/execution/rep.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/d0_Aged.RNA.bigwig \
--regionBodyLength 15000 \
--beforeRegionStartLength 15000 --afterRegionStartLength 15000 \
--outFileName aged.domain.flank.foldchange.mat.gz \
--binSize 1000 \
-p 55

plotHeatmap -m aged.domain.flank.foldchange.mat.gz \
--outFileName aged.domain.flank.foldchange.png \
--outFileNameMatrix aged.domain.flank.foldchange.txt \
--xAxisLabel "distance (bp)" \
--samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 RNA \
--dpi 300

computeMatrix scale-regions \
-R $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-S $prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/get_signal_bigwigs/macs2_signal/Young_H3K27me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig \
/nas/homes/annashch/H4K20me1/caper_out/chip/72cab657-4520-4eee-b95a-d87eac1cecbf/call-macs2_signal_track_pooled/execution/rep.pooled.fc.signal.bigwig \
$prefix/HistoneTracks/d0_Young.RNA.bigwig \
--regionBodyLength 15000 \
--beforeRegionStartLength 15000 --afterRegionStartLength 15000 \
--outFileName young.domain.flank.foldchange.mat.gz \
--binSize 1000 \
-p 55

plotHeatmap -m young.domain.flank.foldchange.mat.gz \
--outFileName young.domain.flank.foldchange.png \
--outFileNameMatrix young.domain.flank.foldchange.txt \
--xAxisLabel "distance (bp)" \
--samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 RNA \
--dpi 300

####################### RPKM tracks
# computeMatrix scale-regions \
# -R $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# -S $trackPrefix/H4K20me1.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/H3K27me3.count.rpkm.track.nodup.Aged.bw \
# $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
# --regionBodyLength 15000 \
# --beforeRegionStartLength 15000 --afterRegionStartLength 15000 \
# --outFileName aged.domain.flank.mat.gz \
# --binSize 1000 \
# -p 50

# plotHeatmap -m aged.domain.flank.mat.gz \
# --outFileName aged.domain.flank.png \
# --outFileNameMatrix aged.domain.flank.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
# --dpi 300

# computeMatrix scale-regions \
# -R $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# -S $trackPrefix/H4K20me1.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/atac.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/H3K27me3.count.rpkm.track.nodup.Young.bw \
# $trackPrefix/H3K4me3.count.rpkm.track.nodup.Young.bw \
# --regionBodyLength 15000 \
# --beforeRegionStartLength 15000 --afterRegionStartLength 15000 \
# --outFileName young.domain.flank.mat.gz \
# --binSize 1000 \
# -p 50

# plotHeatmap -m young.domain.flank.mat.gz \
# --outFileName young.domain.flank.png \
# --outFileNameMatrix young.domain.flank.txt \
# --xAxisLabel "distance (bp)" \
# --samplesLabel H4K20me1 ATAC H3K27me3 H3K4me3 \
# --dpi 300
