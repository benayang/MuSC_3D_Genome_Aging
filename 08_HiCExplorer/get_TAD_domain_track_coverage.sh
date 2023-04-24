prefix=/nas/homes/benyang/HiC

trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

###### Fold-Change
# multiBigwigSummary BED-file -p 50 \
# --bwfiles "/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
# $trackPrefix/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
# --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
# --labels ATAC H3K4me3 \
# --BED $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# -out aged_domain_scores.npz --outRawCounts aged_domain_scores.tab

# multiBigwigSummary BED-file -p 50 \
# --bwfiles "/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
# $trackPrefix/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
# --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
# --labels ATAC H3K4me3 \
# --BED $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
# -out young_domain_scores.npz --outRawCounts young_domain_scores.tab

###### RPKM
multiBigwigSummary BED-file -p 50 \
--bwfiles $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
$trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
--blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
--labels ATAC H3K4me3 \
--BED $prefix/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-out aged_domain_scores_RPKM.npz --outRawCounts aged_domain_scores_RPKM.tab

multiBigwigSummary BED-file -p 50 \
--bwfiles $trackPrefix/atac.count.rpkm.track.nodup.Young.bw \
$trackPrefix/H3K4me3.count.rpkm.track.nodup.Young.bw \
--blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
--labels ATAC H3K4me3 \
--BED $prefix/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed \
-out young_domain_scores_RPKM.npz --outRawCounts young_domain_scores_RPKM.tab