prefix=/nas/homes/benyang/HiC

trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

for f in "A_to_B" "B_to_A" "StaticA" "StaticB" "aged.A" "aged.B" "young.A" "young.B"
do

    multiBigwigSummary BED-file -p 50 \
    --bwfiles $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
    $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$f.bed \
    -out aged_scores_per_$f.npz --outRawCounts aged_scores_per_$f.tab

    multiBigwigSummary BED-file -p 50 \
    --bwfiles $trackPrefix/atac.count.rpkm.track.nodup.Young.bw \
    $trackPrefix/H3K4me3.count.rpkm.track.nodup.Young.bw \
    --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
    --labels ATAC H3K4me3 \
    --BED $prefix/04_FANC/compartmentExpression/compartmentBed/$f.bed \
    -out young_scores_per_$f.npz --outRawCounts young_scores_per_$f.tab

done