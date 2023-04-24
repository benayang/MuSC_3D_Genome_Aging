mainDir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C"
H4K20me1="$mainDir/11_HistoneModsYoung-H4K20me1-trimmed-pooled_005.stringent.bed"
genebodyDir="$mainDir/04_FANC/compartmentExpression/compartmentBed/100kb/genebodies"
outDir="$mainDir/08_HiCExplorer/TAD histone overlap"

for f in A_to_B B_to_A static young.A young.B
do
bedtools intersect -wa \
-a "$genebodyDir/${f}.genebodies.bed" \
-b "$mainDir/11_HistoneMods/GSE129749_yqsc_h3k27me3_peaks.broadPeak.gz" > "$outDir/$f.TAD.GSE129749.H3K27me3.genebodies.bed"

cut -f 4 "$outDir/$f.TAD.GSE129749.H3K27me3.genebodies.bed" | sort | uniq > "$outDir/$f.TAD.GSE129749.H3K27me3.genebodies.genenames.txt"

done