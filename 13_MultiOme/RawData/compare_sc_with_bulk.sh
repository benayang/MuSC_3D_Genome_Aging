# Young
multiBigwigSummary bins \
-b all_reps_combined_young_MuSC_atac_possorted.rpkm.bw \
/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw \
--labels single-cell bulk \
-p 45 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-o young_bigwig_summary_results.npz

plotCorrelation --corData young_bigwig_summary_results.npz \
-c pearson \
-p scatterplot \
-o young_bigwig_summary.png \
--log1p \
--removeOutliers \
--plotFileFormat png

# Aged
multiBigwigSummary bins \
-b all_reps_aged_MuSC_atac_possorted.rpkm.bw \
/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw \
--labels single-cell bulk \
-p 45 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-o aged_bigwig_summary_results.npz

plotCorrelation --corData aged_bigwig_summary_results.npz \
-c pearson \
-p scatterplot \
-o aged_bigwig_summary.png \
--log1p \
--removeOutliers \
--plotFileFormat png