file_dir='/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR/bigwigs'

# multiBamSummary \
# --bamfiles $file_dir/Aged_MuSC.bam \
# /nas/homes/annashch/H4K20me1/atac.merged.Old.bam \
# --outFileName $file_dir/aged_MuSC_atac_bam.npz \
# --binSize 1000 \
# --labels scATAC ATAC \
# -p 45 \
# -v --samFlagExclude 780 --minMappingQuality 30 \
# -bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed

# multiBamSummary \
# --bamfiles $file_dir/combined_young_MuSC.bam \
# /nas/homes/annashch/H4K20me1/atac.merged.Young.bam \
# --outFileName $file_dir/young_MuSC_atac_bam.npz \
# --binSize 1000 \
# --labels scATAC ATAC \
# -p 45 \
# -v --samFlagExclude 780 --minMappingQuality 30 \
# -bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed

plotCorrelation \
--corData $file_dir/aged_MuSC_atac_bam.npz \
--corMethod pearson \
--log1p \
--skipZeros \
--whatToPlot scatterplot \
--plotFile $file_dir/aged_MuSC_atac_bam_pearson_log1p.png \
--labels Single-cell Bulk \
--plotFileFormat png \
--plotHeight 10 --plotWidth 10

plotCorrelation \
--corData $file_dir/young_MuSC_atac_bam.npz \
--corMethod pearson \
--log1p \
--skipZeros \
--whatToPlot scatterplot \
--plotFile $file_dir/young_MuSC_atac_bam_pearson_log1p.png \
--labels Single-cell Bulk \
--plotFileFormat png \
--plotHeight 10 --plotWidth 10