mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"

cat "$mainDir/04_FANC/without KR normalization/aged.ab_100kb.bed" | \
sort -k1,1 -k2,2n | \
bedtools closest -a stdin \
-b "$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed" -d > "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/tss/aged.ab.gene.bfilt.distance.bed"

cat "$mainDir/04_FANC/without KR normalization/young.ab_100kb.bed" | \
sort -k1,1 -k2,2n | \
bedtools closest -a stdin \
-b "$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed" -d > "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/tss/young.ab.gene.bfilt.distance.bed"

cat "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/A_to_B.bed" | \
sort -k1,1 -k2,2n | \
bedtools closest -a stdin \
-b "$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed" -d > "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/tss/A_to_B.gene.bfilt.distance.bed"

cat "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/B_to_A.bed" | \
sort -k1,1 -k2,2n | \
bedtools closest -a stdin \
-b "$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed" -d > "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/tss/B_to_A.gene.bfilt.distance.bed"

cat "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/static.bed" | \
sort -k1,1 -k2,2n | \
bedtools closest -a stdin \
-b "$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed" -d > "$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/tss/static.gene.bfilt.distance.bed"