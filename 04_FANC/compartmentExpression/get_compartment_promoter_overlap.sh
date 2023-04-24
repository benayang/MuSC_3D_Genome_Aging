mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
bedDir="$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb"
# bedDir="$bedDir/sorted"
tssDir="$bedDir/tss"
# young_A="$bedDir/young.A.bed"
# young_B="$bedDir/young.B.bed"
# aged_A="$bedDir/aged.A.bed"
# aged_B="$bedDir/aged.B.bed"
# A_to_B="$bedDir/A_to_B.bed"
# B_to_A="$bedDir/B_to_A.bed"
# static="$bedDir/static.bed"
mm10="$mainDir/sizes.mm10"
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.sorted.bed"

# sort -k1,1 -k2,2n "$young_A" > "$bedDir/young.A.bed"
# sort -k1,1 -k2,2n "$young_B" > "$bedDir/young.B.bed"
# sort -k1,1 -k2,2n "$aged_A" > "$bedDir/aged.A.bed"
# sort -k1,1 -k2,2n "$aged_B" > "$bedDir/aged.B.bed"

# sort -k1,1 -k2,2n "$A_to_B" > "$bedDir/A_to_B.bed"
# sort -k1,1 -k2,2n "$B_to_A" > "$bedDir/B_to_A.bed"
# sort -k1,1 -k2,2n "$static" > "$bedDir/static.bed"

bedtools slop -b 1000 -i "$tss" -g "$mm10" | \
bedtools intersect -wa -a stdin -b "$bedDir/young.A.bed" > "$tssDir/young.A.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/young.B.bed" > "$tssDir/young.B.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/aged.A.bed" > "$tssDir/aged.A.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/aged.B.bed" > "$tssDir/aged.B.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/A_to_B.bed" > "$tssDir/A_to_B.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/B_to_A.bed" > "$tssDir/B_to_A.1kb.promoter.bed"

bedtools slop -i "$tss" -g "$mm10" -b 1000 | \
bedtools intersect -wa -a stdin -b "$bedDir/static.bed" > "$tssDir/static.1kb.promoter.bed"
