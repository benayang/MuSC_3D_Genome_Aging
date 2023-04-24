#!/usr/bin/env bash

outdir="/nas/homes/benyang/HiC/04_FANC/compartmentExpression/ATAC_overlap/compartments"
maindir="/nas/homes/benyang/HiC/04_FANC/compartmentExpression"
young_ATAC="/nas/homes/benyang/HiC/04_FANC/compartmentExpression/ATAC_overlap/atac.d0.Young.optimal.narrowPeak.gz"
aged_ATAC="/nas/homes/benyang/HiC/04_FANC/compartmentExpression/ATAC_overlap/atac.d0.Aged.optimal.narrowPeak.gz"
young_A="$maindir/compartmentBed/young.A.bed"
young_B="$maindir/compartmentBed/young.B.bed"
aged_A="$maindir/compartmentBed/aged.A.bed"
aged_B="$maindir/compartmentBed/aged.B.bed"
A_to_B="$maindir/A_to_B.bed"
B_to_A="$maindir/B_to_A.bed"
static="$maindir/static.bed"

for f in $A_to_B $B_to_A $static $young_A $young_B $aged_A $aged_B
do
fname=$(basename $f | rev | cut -c 5- | rev)		
bedtools intersect -a $f -b $young_ATAC > $outdir/young_ATAC_$fname.bed
bedtools intersect -a $f -b $aged_ATAC > $outdir/aged_ATAC_$fname.bed

cat $outdir/young_ATAC_$fname.bed | wc -l 
cat $outdir/aged_ATAC_$fname.bed | wc -l

done
