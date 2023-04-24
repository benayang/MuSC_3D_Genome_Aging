prefix=/nas/homes/benyang/HiC

trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs
outdir=$prefix/05_Arrowhead_HICCUPS/hiccups_diff/loop_domain_track_coverage

for f in stable aged young
do

    multiBigwigSummary BED-file -p 50 \
    --bwfiles "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
    "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/05_Arrowhead_HICCUPS/olap_${f}_loop_domains.bed \
    -out $outdir/olap_analysis/aged_scores_per_olap_${f}_loop.npz --outRawCounts $outdir/olap_analysis/aged_scores_per_olap_${f}_loop.tab

    multiBigwigSummary BED-file -p 50 \
    --bwfiles "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
    "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/05_Arrowhead_HICCUPS/olap_${f}_loop_domains.bed \
    -out $outdir/olap_analysis/young_scores_per_olap_${f}_loop.npz --outRawCounts $outdir/olap_analysis/young_scores_per_olap_${f}_loop.tab

done

# for f in young_diff aged_diff young_non_diff aged_non_diff
# do

#     multiBigwigSummary BED-file -p 50 \
#     --bwfiles "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
#     "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
#     --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
#     --labels ATAC H3K4me3 \
#     --BED $prefix/05_Arrowhead_HICCUPS/hiccups_diff/${f}_loop_domains.bed \
#     -out $outdir/diff_analysis/aged_scores_per_${f}_loop.npz --outRawCounts $outdir/diff_analysis/aged_scores_per_${f}_loop.tab

#     multiBigwigSummary BED-file -p 50 \
#     --bwfiles "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
#     "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
#     --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
#     --labels ATAC H3K4me3 \
#     --BED $prefix/05_Arrowhead_HICCUPS/hiccups_diff/${f}_loop_domains.bed \
#     -out $outdir/diff_analysis/young_scores_per_${f}_loop.npz --outRawCounts $outdir/diff_analysis/young_scores_per_${f}_loop.tab

# done

# for f in young_diff aged_diff young_non_diff aged_non_diff
# do

    # multiBigwigSummary BED-file -p 50 \
    # --bwfiles "/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
    # $trackPrefix/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
    # --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    # --labels ATAC H3K4me3 \
    # --BED $prefix/05_Arrowhead_HICCUPS/hiccups_diff/${f}_loop_domains.bed \
    # -out $outdir/aged_scores_per_${f}_loop.npz --outRawCounts $outdir/aged_scores_per_${f}_loop.tab

    # multiBigwigSummary BED-file -p 50 \
    # --bwfiles "/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
    # $trackPrefix/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
    # --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    # --labels ATAC H3K4me3 \
    # --BED $prefix/05_Arrowhead_HICCUPS/hiccups_diff/${f}_loop_domains.bed \
    # -out $outdir/young_scores_per_${f}_loop.npz --outRawCounts $outdir/young_scores_per_${f}_loop.tab

# done