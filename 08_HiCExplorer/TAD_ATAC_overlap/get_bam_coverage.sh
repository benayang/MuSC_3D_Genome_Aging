mainDir="/nas/homes/benyang/HiC/"

young_promoter_boundary=$mainDir+'08_HiCExplorer/TAD expression/young.merged.TAD.promoter.boundaries.bed'
aged_promoter_boundary=$mainDir+'08_HiCExplorer/TAD expression/aged.merged.TAD.promoter.boundaries.bed'
young_non_promoter_boundary=$mainDir+'08_HiCExplorer/TAD expression/young.merged.TAD.nonpromoter.boundaries.bed'
aged_non_promoter_boundary=$mainDir+'08_HiCExplorer/TAD expression/aged.merged.TAD.nonpromoter.boundaries.bed'

young_ATAC='/media/drive1/annashch/nobel_lab_projects/age_V2/tracks/bams/merged_bams/d0_Young.merged.bam'
aged_ATAC='/media/drive1/annashch/nobel_lab_projects/age_V2/tracks/bams/merged_bams/d0_Aged.merged.bam'

# bamCoverage -b $young_ATAC \
# --binSize 1
# --region $young_promoter_boundary \
# -p 50 \
# --normalizeUsing RPKM \
# --minMappingQuality 30 \
# --samFlagExclude 780 \
# -v \
# -o young_ATAC_promoter_boundary.bw

