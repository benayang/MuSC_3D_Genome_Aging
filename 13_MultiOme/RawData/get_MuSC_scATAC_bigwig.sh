prefix='/nas/homes/benyang/HiC/13_MultiOme'

########## Aged
sinto filterbarcodes \
-b $prefix/RawData/Aged/10x_analysis_5422-JL/Sample_5422-JL-1/atac_possorted_bam.bam \
-c $prefix/RawData/all_reps_aged_MuSC_cell_barcodes.tsv \
-p 50 \
--outdir $prefix/RawData

samtools index $prefix/RawData/all_reps_aged_MuSC.bam

# read length is 50bp
bamCoverage \
-b $prefix/RawData/all_reps_aged_MuSC.bam \
-o $prefix/RawData/all_reps_aged_MuSC_atac_possorted.rpkm.bw \
-of bigwig -bs 1 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-p 50 --normalizeUsing RPKM --effectiveGenomeSize 2308125349

########## Young
sinto filterbarcodes \
-b $prefix/RawData/Young/10x_analysis_5644-JL/Sample_5644-JL-1/atac_possorted_bam.bam \
-c $prefix/RawData/all_reps_young_MuSC_cell_barcodes.tsv \
-p 50 \
--outdir $prefix/RawData

samtools index $prefix/RawData/all_reps_young_MuSC.bam

# read length is 50bp
bamCoverage \
-b $prefix/RawData/all_reps_young_MuSC.bam \
-o $prefix/RawData/all_reps_young_MuSC_atac_possorted.rpkm.bw \
-of bigwig -bs 1 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-p 50 --normalizeUsing RPKM --effectiveGenomeSize 2308125349

########## Young_v2
sinto filterbarcodes \
-b /nas/datasets/6098-JL/10x_analysis_6098-JL/Sample_6098-JL-1/atac_possorted_bam.bam \
-c $prefix/RawData/all_reps_young_v2_MuSC_cell_barcodes.tsv \
-p 50 \
--outdir $prefix/RawData

samtools index $prefix/RawData/all_reps_young_v2_MuSC.bam

# read length is 50bp
bamCoverage \
-b $prefix/RawData/all_reps_young_v2_MuSC.bam \
-o $prefix/RawData/all_reps_young_v2_MuSC_atac_possorted.rpkm.bw \
-of bigwig -bs 1 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-p 50 --normalizeUsing RPKM --effectiveGenomeSize 2308125349

######### Combine young bam files
samtools merge -@ 45 all_reps_combined_young_MuSC.bam "/nas/homes/benyang/HiC/13_MultiOme/RawData/all_reps_young_MuSC.bam" "/nas/homes/benyang/HiC/13_MultiOme/RawData/all_reps_young_v2_MuSC.bam"

bamCoverage \
-b $prefix/RawData/all_reps_combined_young_MuSC.bam \
-o $prefix/RawData/all_reps_combined_young_MuSC_atac_possorted.rpkm.bw \
-of bigwig -bs 1 \
-bl /nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed \
-p 50 --normalizeUsing RPKM --effectiveGenomeSize 2308125349