prefix=/nas/homes/benyang/HiC

for f in "A" "B" "A_to_B" "B_to_A" "StaticA" "StaticB"
do

bedtools flank -i $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADBoundary.bed \
-g /nas/homes/benyang/Genome_References/sizes.mm10 -b 40000 > aged.$f.TADBoundary.flank.bed

bedtools flank -i $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADBoundary.bed \
-g /nas/homes/benyang/Genome_References/sizes.mm10 -b 40000 > young.$f.TADBoundary.flank.bed

done