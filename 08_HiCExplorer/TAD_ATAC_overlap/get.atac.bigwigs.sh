#!/bin/bash
prefix='/nas/homes/benyang/HiC/08_HiCExplorer/TAD_ATAC_overlap/'

samtools merge --threads 50 $prefix/atac.merged.Aged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-1/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R3_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-3/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_21K_R1_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-6/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_R2_129772_CTCTCTAC_S7_R1_001.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-0/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R2_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-5/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_R1_129772_TAGGCATG_S6_R1_001.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-2/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/MuSC_d0_Old_10K_R1_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Aged/call-bowtie2/shard-4/execution/glob-3bcbe4e7489c90f75e0523ac6f3a9385/D0_Old_MuSC_AGGCAGAA_S3_R1_001.merged.bam

samtools index $prefix/atac.merged.Aged.bam

bamCoverage -p50 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/atac.merged.Aged.bam -o $prefix/atac.count.rpkm.track.Aged.bw 

samtools merge --threads 50 $prefix/atac.merged.Young.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-bowtie2/shard-0/execution/MuSC_d0_Y_10K_R1_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-bowtie2/shard-1/execution/MuSC_d0_Y_10K_R2_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-bowtie2/shard-2/execution/MuSC_d0_Y_10K_R2_r1.merged.bam /nas/homes/annashch/Age_ATAC/outputs/d0_Young/cromwell-executions/atac/e373e737-4927-458b-b2cb-023ed08dd76e/call-bowtie2/shard-3/execution/127397_CGTACTAG-GTGTAGAT_S1_R1_001.merged.bam 

samtools index $prefix/atac.merged.Young.bam

bamCoverage -p50 -v --normalizeUsing RPKM --binSize 1 --samFlagExclude 780 --minMappingQuality 30 -b $prefix/atac.merged.Young.bam -o $prefix/atac.count.rpkm.track.Young.bw 
