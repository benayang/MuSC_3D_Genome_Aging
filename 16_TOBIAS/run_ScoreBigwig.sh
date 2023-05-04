prefix='/nas/homes/benyang/HiC'

# Individual peak sets
# TOBIAS ScoreBigwig \
# --signal $prefix/16_TOBIAS/aged_merged_ATAC_corrected.bw \
# --output $prefix/16_TOBIAS/aged_merged_ATAC_bigwig_score.bw \
# --regions $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak \
# --cores 50

# TOBIAS ScoreBigwig \
# --signal $prefix/16_TOBIAS/young_merged_ATAC_corrected.bw \
# --output $prefix/16_TOBIAS/young_merged_ATAC_bigwig_score.bw \
# --regions $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak \
# --cores 50

# Merged peak sets for differential detection via BINDetect
TOBIAS ScoreBigwig \
--signal $prefix/16_TOBIAS/aged_merged_ATAC_corrected.bw \
--output $prefix/16_TOBIAS/aged_merged_peakset_ATAC_bigwig_score.bw \
--regions $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
--cores 50

TOBIAS ScoreBigwig \
--signal $prefix/16_TOBIAS/young_merged_ATAC_corrected.bw \
--output $prefix/16_TOBIAS/young_merged_peakset_ATAC_bigwig_score.bw \
--regions $prefix/16_TOBIAS/young_aged_merged_ATAC_peaks.bed \
--cores 50