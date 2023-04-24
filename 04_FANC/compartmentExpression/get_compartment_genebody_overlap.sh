mainDir="/home/benjy/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
bedDir="$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb"
sortedDir="$bedDir/sorted"
genesDir="$bedDir/genebodies"
young_A="$bedDir/young.A.bed"
young_B="$bedDir/young.B.bed"
aged_A="$bedDir/aged.A.bed"
aged_B="$bedDir/aged.B.bed"
A_to_B="$bedDir/A_to_B.bed"
B_to_A="$bedDir/B_to_A.bed"
static="$bedDir/static.bed"
mm10="$mainDir/sizes.mm10"
genebodies="$mainDir/get_tss/genebodies.gencode.vM25.basic.annotation.filtered.uniq.bed"

sort -k1,1 -k2,2n "$young_A" > "$sortedDir/young.A.sorted.bed"
sort -k1,1 -k2,2n "$young_B" > "$sortedDir/young.B.sorted.bed"
sort -k1,1 -k2,2n "$aged_A" > "$sortedDir/aged.A.sorted.bed"
sort -k1,1 -k2,2n "$aged_B" > "$sortedDir/aged.B.sorted.bed"

sort -k1,1 -k2,2n "$A_to_B" > "$sortedDir/A_to_B.sorted.bed"
sort -k1,1 -k2,2n "$B_to_A" > "$sortedDir/B_to_A.sorted.bed"
sort -k1,1 -k2,2n "$static" > "$sortedDir/static.sorted.bed"


bedtools intersect -wa -a "$genebodies" -b "$sortedDir/young.A.sorted.bed" > "$genesDir/young.A.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/young.B.sorted.bed" > "$genesDir/young.B.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/aged.A.sorted.bed" > "$genesDir/aged.A.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/aged.B.sorted.bed" > "$genesDir/aged.B.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/A_to_B.sorted.bed" > "$genesDir/A_to_B.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/B_to_A.sorted.bed" > "$genesDir/B_to_A.genebodies.bed"

bedtools intersect -wa -a "$genebodies" -b "$sortedDir/static.sorted.bed" > "$genesDir/static.genebodies.bed"
