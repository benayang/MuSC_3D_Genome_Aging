#!/usr/bin/env bash

outdir="/nas/homes/benyang/HiC/04_FANC/compartmentExpression"
young_ATAC="/nas/homes/benyang/HiC/04_FANC/compartmentExpression/ATAC_overlap/atac.d0.Young.optimal.narrowPeak.gz"
aged_ATAC="/nas/homes/benyang/HiC/04_FANC/compartmentExpression/ATAC_overlap/atac.d0.Aged.optimal.narrowPeak.gz"
young_A="$outdir/compartmentBed/young.A.bed"
young_B="$outdir/compartmentBed/young.B.bed"
aged_A="$outdir/compartmentBed/aged.A.bed"
aged_B="$outdir/compartmentBed/aged.B.bed"
A_to_B="$outdir/A_to_B.bed"
B_to_A="$outdir/B_to_A.bed"
static="$outdir/static.bed"
mm10="/media/labuser/BenPassport/GenomeRefs/sizes.mm10"
tss="/media/labuser/BenPassport/Datasets/HiC/fastqs_4140-KS/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.bed"

bedtools slop -i $tss -g $mm10 -b 1000 > $outdir/mm10.tss.1kb.slop.bed

# A to B
#bedtools intersect -a $A_to_B -b $young_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/young_ATAC_A_to_B.1kb.tss.bed

#bedtools intersect -a $A_to_B -b $aged_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/aged_ATAC_A_to_B.1kb.tss.bed
# B to A
#bedtools intersect -a $B_to_A -b $young_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/young_ATAC_B_to_A.1kb.tss.bed

#bedtools intersect -a $B_to_A -b $aged_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/aged_ATAC_B_to_A.1kb.tss.bed
# static
#bedtools intersect -a $static -b $young_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/young_ATAC_static.1kb.tss.bed

#bedtools intersect -a $static -b $aged_ATAC | \
#	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/aged_ATAC_static.1kb.tss.bed
# A and B
bedtools intersect -a $young_A -b $young_ATAC | \
	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/young_ATAC_A.1kb.tss.bed

bedtools intersect -a $young_B -b $young_ATAC | \
	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/young_ATAC_B.1kb.tss.bed

bedtools intersect -a $aged_A -b $aged_ATAC | \
	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/aged_ATAC_A.1kb.tss.bed

bedtools intersect -a $aged_B -b $aged_ATAC | \
	bedtools intersect -wa -a $outdir/mm10.tss.1kb.slop.bed -b stdin > $outdir/ATAC_overlap/aged_ATAC_B.1kb.tss.bed

for f in young_ATAC_A_to_B.1kb.tss.bed aged_ATAC_A_to_B.1kb.tss.bed young_ATAC_B_to_A.1kb.tss.bed aged_ATAC_B_to_A.1kb.tss.bed young_ATAC_static.1kb.tss.bed aged_ATAC_static.1kb.tss.bed young_ATAC_A.1kb.tss.bed young_ATAC_B.1kb.tss.bed aged_ATAC_A.1kb.tss.bed aged_ATAC_B.1kb.tss.bed
do
cat $outdir/ATAC_overlap/$f | wc -l
done
