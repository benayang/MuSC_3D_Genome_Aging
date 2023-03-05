prefix=/nas/homes/benyang/HiC

trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

for f in "A_to_B" "B_to_A" "StaticA" "StaticB" "aged.A" "aged.B" "young.A" "young.B" "static"
do

    multiBigwigSummary BED-file -p 50 \
    --bwfiles "/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-macs2_signal_track_pooled/execution/MuSC_d0_Old_10K_R2_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
    $trackPrefix/macs2_signal/Aged_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$f.bed \
    -out aged_scores_per_$f.npz --outRawCounts aged_scores_per_$f.tab

    multiBigwigSummary BED-file -p 50 \
    --bwfiles "/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-macs2_signal_track_pooled/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.tn5.pooled.fc.signal.bigwig" \
    $trackPrefix/macs2_signal/Young_H3K4me3/rep.pooled_x_ctl.pooled.fc.signal.bigwig \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$f.bed \
    -out young_scores_per_$f.npz --outRawCounts young_scores_per_$f.tab

done
