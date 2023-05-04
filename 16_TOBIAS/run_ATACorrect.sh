prefix='/nas/homes/benyang/HiC'

TOBIAS ATACorrect \
--bam $prefix/HistoneTracks/get_signal_bigwigs/atac.merged.nodup.Young.bam \
--genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
--peaks $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak \
--blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
--outdir $prefix/16_TOBIAS \
--prefix young_merged_ATAC \
--cores 50

TOBIAS ATACorrect \
--bam $prefix/HistoneTracks/get_signal_bigwigs/atac.merged.nodup.Aged.bam \
--genome "/nas/homes/benyang/Genome_References/mm10_no_alt_analysis_set_ENCODE.fasta" \
--peaks $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak \
--blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" \
--outdir $prefix/16_TOBIAS \
--prefix aged_merged_ATAC \
--cores 50