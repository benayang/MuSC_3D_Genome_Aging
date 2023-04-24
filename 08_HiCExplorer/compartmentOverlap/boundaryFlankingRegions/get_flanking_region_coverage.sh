prefix=/nas/homes/benyang/HiC
trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

for f in "A" "B" "A_to_B" "B_to_A" "StaticA" "StaticB"
do

    multiBigwigSummary BED-file -p 55 \
    --bwfiles $trackPrefix/H4K20me1.count.rpkm.track.nodup.Aged.bw \
    $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
    $trackPrefix/H3K27me3.count.rpkm.track.nodup.Aged.bw \
    $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
    --labels H4K20me1 ATAC H3K27me3 H3K4me3 \
    --BED $prefix/08_HiCExplorer/compartmentOverlap/boundary_flanking_regions/aged.$f.TADBoundary.flank.bed \
    -out get_flank_coverage/aged_binned_scores_per_$f.npz --outRawCounts get_flank_coverage/aged_binned_scores_per_$f.tab \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz

    multiBigwigSummary BED-file -p 55 \
    --bwfiles $trackPrefix/H4K20me1.count.rpkm.track.nodup.Young.bw \
    $trackPrefix/atac.count.rpkm.track.nodup.Young.bw \
    $trackPrefix/H3K27me3.count.rpkm.track.nodup.Young.bw \
    $trackPrefix/H3K4me3.count.rpkm.track.nodup.Young.bw \
    --labels H4K20me1 ATAC H3K27me3 H3K4me3 \
    --BED $prefix/08_HiCExplorer/compartmentOverlap/boundary_flanking_regions/young.$f.TADBoundary.flank.bed \
    -out get_flank_coverage/young_binned_scores_per_$f.npz --outRawCounts get_flank_coverage/young_binned_scores_per_$f.tab \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz


done