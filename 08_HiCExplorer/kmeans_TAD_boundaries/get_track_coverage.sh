prefix=/nas/homes/benyang/HiC

outdir=$prefix/08_HiCExplorer/kmeans_TAD_boundaries/track_coverage
trackPrefix=$prefix/HistoneTracks/get_signal_bigwigs

for age in young aged
do
    for c in `seq 1 1 6`
    do
        
        # Boundary proper
        # multiBigwigSummary BED-file -p 50 \
        # --bwfiles $trackPrefix/H4K20me1.count.rpkm.track.nodup.Aged.bw \
        # $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
        # $trackPrefix/H3K27me3.count.rpkm.track.nodup.Aged.bw \
        # $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
        # --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
        # --labels H4K20me1 ATAC H3K27me3 H3K4me3 \
        # --BED $prefix/08_HiCExplorer/kmeans_TAD_boundaries/${age}_cluster$c.bedgraph \
        # -out $outdir/${age}_scores_per_cluster$c.npz --outRawCounts $outdir/${age}_scores_per_cluster$c.tab
        # Boundary flanking
        bedtools flank -i $prefix/08_HiCExplorer/kmeans_TAD_boundaries/${age}_cluster$c.bedgraph \
        -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 40000 | sort -k1,1 -k2,2n > $outdir/${age}_cluster${c}_flank.bed
        
        multiBigwigSummary BED-file -p 50 \
        --bwfiles $trackPrefix/H4K20me1.count.rpkm.track.nodup.Aged.bw \
        $trackPrefix/atac.count.rpkm.track.nodup.Aged.bw \
        $trackPrefix/H3K27me3.count.rpkm.track.nodup.Aged.bw \
        $trackPrefix/H3K4me3.count.rpkm.track.nodup.Aged.bw \
        --blackListFileName /nas/homes/benyang/JC_H3K27me3/encode_genome_data/ENCFF547MET.bed.gz \
        --labels H4K20me1 ATAC H3K27me3 H3K4me3 \
        --BED $outdir/${age}_cluster${c}_flank.bed \
        -out $outdir/${age}_scores_per_cluster${c}_flank.npz --outRawCounts $outdir/${age}_scores_per_cluster${c}_flank.tab


    done
done