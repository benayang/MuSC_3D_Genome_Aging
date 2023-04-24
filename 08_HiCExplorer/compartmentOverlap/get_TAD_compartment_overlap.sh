mainDir="/nas/homes/benyang/HiC"
mm10="/nas/homes/benyang/Genome_References/sizes.mm10"
tss="$mainDir/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed"
suffix="min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"
boundarysuffix="min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed"
outDir="$mainDir/08_HiCExplorer/compartmentOverlap"

# get boundary-compartment overlap
for f in "young.A" "young.B" "aged.A" "aged.B"
do 

age=${f%.*}

bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/$age.merged/40kb/$age.merged_$boundarysuffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/$f.TADBoundary.bed"

done

for f in "A_to_B" "B_to_A" "StaticA" "StaticB"
do

bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.$f.TADBoundary.bed"

bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.$f.TADBoundary.bed"

done
# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/young.A.bed" \
# -b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" > \
# "$outDir/young.A.TADBoundary.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/young.B.bed" \
# -b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$boundarysuffix" > \
# "$outDir/young.B.TADBoundary.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/aged.A.bed" \
# -b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" > \
# "$outDir/aged.A.TADBoundary.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/aged.B.bed" \
# -b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$boundarysuffix" > \
# "$outDir/aged.B.TADBoundary.bed"

# get domain-compartment overlap
for f in "young.A" "young.B" "aged.A" "aged.B"
do 

age=${f%.*}

bedtools intersect -wa -f 0.5 \
-a "$mainDir/08_HiCExplorer/$age.merged/40kb/$age.merged_$suffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/$f.TADDomain.bed"

done

for f in "A_to_B" "B_to_A" "StaticA" "StaticB"
do

bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$suffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/aged.$f.TADDomain.bed"

bedtools intersect -wa \
-a "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$suffix" \
-b "$mainDir/04_FANC/compartmentExpression/compartmentBed/$f.bed" | \
sort -k1,1 -k2,2n | uniq > "$outDir/young.$f.TADDomain.bed"

done

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/young.A.bed" \
# -b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$suffix" > \
# "$outDir/young.A.TADDomain.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/young.B.bed" \
# -b "$mainDir/08_HiCExplorer/young.merged/40kb/young.merged_$suffix" > \
# "$outDir/young.B.TADDomain.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/aged.A.bed" \
# -b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$suffix" > \
# "$outDir/aged.A.TADDomain.bed"

# bedtools intersect -wa \
# -a "$mainDir/04_FANC/compartmentExpression/compartmentBed/aged.B.bed" \
# -b "$mainDir/08_HiCExplorer/aged.merged/40kb/aged.merged_$suffix" > \
# "$outDir/aged.B.TADDomain.bed"