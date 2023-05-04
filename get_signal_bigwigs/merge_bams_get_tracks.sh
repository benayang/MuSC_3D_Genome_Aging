#!/bin/bash
prefix=/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs

######### ATAC
#aged_ATAC_prefix=/nas/homes/annashch/Age_ATAC/caper_out/atac/b4245fba-248f-4201-bd66-907bb8e97a85/call-filter
#young_ATAC_prefix=/nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-filter

# samtools merge --threads 30 $prefix/atac.merged.nodup.Aged.bam \
# $aged_ATAC_prefix/shard-0/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R2_r1.merged.nodup.bam \
# $aged_ATAC_prefix/shard-1/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R3_r1.merged.nodup.bam \
# $aged_ATAC_prefix/shard-2/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R1_r1.merged.nodup.bam \
# $aged_ATAC_prefix/shard-3/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_21K_R1_r1.merged.nodup.bam \
# $aged_ATAC_prefix/shard-4/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_AGGCAGAA_S3_R1_001.merged.nodup.bam \
# $aged_ATAC_prefix/shard-5/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_R1_129772_TAGGCATG_S6_R1_001.merged.nodup.bam \
# $aged_ATAC_prefix/shard-6/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_R2_129772_CTCTCTAC_S7_R1_001.merged.nodup.bam

# samtools index $prefix/atac.merged.nodup.Aged.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/atac.merged.nodup.Aged.bam -o $prefix/atac.count.rpkm.track.nodup.Aged.bw 

# samtools merge --threads 30 $prefix/atac.merged.nodup.Young.bam \
# $young_ATAC_prefix/shard-0/execution/MuSC_d0_Y_10K_R1_r1.merged.nodup.bam \
# $young_ATAC_prefix/shard-1/execution/MuSC_d0_Y_10K_R2_r1.merged.nodup.bam \
# $young_ATAC_prefix/shard-2/execution/MuSC_d0_Y_10K_R2_r1.merged.nodup.bam \
# $young_ATAC_prefix/shard-3/execution/127397_CGTACTAG-GTGTAGAT_S1_R1_001.merged.nodup.bam

# samtools index $prefix/atac.merged.nodup.Young.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/atac.merged.nodup.Young.bam -o $prefix/atac.count.rpkm.track.nodup.Young.bw 

############ H4K20me1

# young_H4K20me1_prefix=/nas/homes/annashch/H4K20me1/caper_out/chip/72cab657-4520-4eee-b95a-d87eac1cecbf/call-filter
# aged_H4K20me1_prefix=/nas/homes/annashch/H4K20me1/caper_out/chip/7aab8e70-c16d-481f-ab6d-ebd0e3683fae/call-filter

# samtools merge --threads 30 $prefix/H4K20me1.merged.nodup.Young.bam \
# $young_H4K20me1_prefix/shard-0/execution/Young-H4K20me1-1_R1_001_val_1.nodup.bam \
# $young_H4K20me1_prefix/shard-1/execution/Young-H4K20me1-2_R1_001_val_1.nodup.bam \
# $young_H4K20me1_prefix/shard-2/execution/Young-H4K20me1-3_R1_001_val_1.nodup.bam

# samtools index $prefix/H4K20me1.merged.nodup.Young.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H4K20me1.merged.nodup.Young.bam -o $prefix/H4K20me1.count.rpkm.track.nodup.Young.bw 

# samtools merge --threads 30 $prefix/H4K20me1.merged.nodup.Aged.bam \
# $aged_H4K20me1_prefix/shard-0/execution/Old-H4K20me1-1_R1_001_val_1.nodup.bam \
# $aged_H4K20me1_prefix/shard-1/execution/Old-H4K20me1-2_R1_001_val_1.nodup.bam \
# $aged_H4K20me1_prefix/shard-2/execution/Old-H4K20me1-3_R1_001_val_1.nodup.bam

# samtools index $prefix/H4K20me1.merged.nodup.Aged.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H4K20me1.merged.nodup.Aged.bam -o $prefix/H4K20me1.count.rpkm.track.nodup.Aged.bw 

############ H3K27me3

# young_H3K27me3_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/croo/croo_Young_H3K27me3/align
# aged_H3K27me3_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/croo/croo_Aged_H3K27me3/align

# samtools merge --threads 30 $prefix/H3K27me3.merged.nodup.Young.bam \
# $young_H3K27me3_prefix/rep1/Young_H3K27me3_Rep1_SRR867402.fastq.srt.nodup.bam \
# $young_H3K27me3_prefix/rep2/Young_H3K27me3_Rep2_SRR867403.fastq.srt.nodup.bam

# samtools index $prefix/H3K27me3.merged.nodup.Young.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H3K27me3.merged.nodup.Young.bam -o $prefix/H3K27me3.count.rpkm.track.nodup.Young.bw 

# samtools merge --threads 30 $prefix/H3K27me3.merged.nodup.Aged.bam \
# $aged_H3K27me3_prefix/rep1/Aged_H3K27me3_Rep1_SRR867416.fastq.srt.nodup.bam \
# $aged_H3K27me3_prefix/rep2/Aged_H3K27me3_Rep2_SRR867417.fastq.srt.nodup.bam

# samtools index $prefix/H3K27me3.merged.nodup.Aged.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H3K27me3.merged.nodup.Aged.bam -o $prefix/H3K27me3.count.rpkm.track.nodup.Aged.bw 

############ H3K4me3

# young_H3K4me3_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/croo/croo_Young_H3K4me3/align
# aged_H3K4me3_prefix=/nas/homes/benyang/HiC/HistoneTracks/Tom_Rando/croo/croo_Aged_H3K4me3/align

# samtools merge --threads 30 $prefix/H3K4me3.merged.nodup.Young.bam \
# $young_H3K4me3_prefix/rep1/Young_H3K4me3_Rep1_SRR867400.fastq.srt.nodup.bam \
# $young_H3K4me3_prefix/rep2/Young_H3K4me3_Rep2_SRR867401.fastq.srt.nodup.bam

# samtools index $prefix/H3K4me3.merged.nodup.Young.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H3K4me3.merged.nodup.Young.bam -o $prefix/H3K4me3.count.rpkm.track.nodup.Young.bw 

# samtools merge --threads 30 $prefix/H3K4me3.merged.nodup.Aged.bam \
# $aged_H3K4me3_prefix/rep1/Aged_H3K4me3_Rep1_SRR867414.fastq.srt.nodup.bam \
# $aged_H3K4me3_prefix/rep2/Aged_H3K4me3_Rep2_SRR867415.fastq.srt.nodup.bam

# samtools index $prefix/H3K4me3.merged.nodup.Aged.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H3K4me3.merged.nodup.Aged.bam -o $prefix/H3K4me3.count.rpkm.track.nodup.Aged.bw 

########## H3K27ac

# H3K27ac_prefix=/nas/homes/benyang/HiC/HistoneTracks/H3K27ac/chip/eafea56b-4971-4718-aaa9-b37b24febaba/call-filter

# samtools merge --threads 30 $prefix/H3K27ac.merged.nodup.bam \
# $H3K27ac_prefix/shard-0/execution/T0_H3K27ac_Rep1_SRR5984122.fastq.srt.nodup.bam \
# $H3K27ac_prefix/shard-1/execution/T0_H3K27ac_Rep2_SRR5984123.fastq.srt.nodup.bam

# samtools index $prefix/H3K27ac.merged.nodup.bam

# bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/H3K27ac.merged.nodup.bam -o $prefix/H3K27ac.count.rpkm.track.nodup.bw 

########## H3K27ac_GSE122867

bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 \
-b /nas/homes/benyang/HiC/HistoneTracks/H3K27ac_GSE122867/align/TA_2mo_H3K27ac_SRR8239713.filt.final.bam \
-o $prefix/H3K27ac_GSE122867.count.rpkm.track.nodup.2mo.bw

bamCoverage -p30 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 \
-b /nas/homes/benyang/HiC/HistoneTracks/H3K27ac_GSE122867/align/TA_20mo_H3K27ac_SRR8239715.filt.final.bam \
-o $prefix/H3K27ac_GSE122867.count.rpkm.track.nodup.20mo.bw