prefix=/nas/homes/benyang/HiC

for f in "A" "B" "A_to_B" "B_to_A" "StaticA" "StaticB"
do

# split TAD domain into 10 bins and add flanking regions of 10 bins with 10kb/bin

while read line 
do
    echo "$line" | bedtools flank -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | head -n 1 | \
    bedtools makewindows -b stdin -n 10 -i srcwinnum >> aged.$f.TADDomain.extended.bed

    echo "$line" | bedtools makewindows -b stdin -n 10 -i srcwinnum >> aged.$f.TADDomain.extended.bed

    echo "$line" | bedtools flank -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | tail -n +2 | \
    bedtools makewindows -b stdin -n 10 -i srcwinnum >> aged.$f.TADDomain.extended.bed
done < $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADDomain.bed

while read line 
do
    echo "$line" | bedtools flank -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | head -n 1 | \
    bedtools makewindows -b stdin -n 10 -i srcwinnum >> young.$f.TADDomain.extended.bed

    echo "$line" | bedtools makewindows -b stdin -n 10 -i srcwinnum >> young.$f.TADDomain.extended.bed

    echo "$line" | bedtools flank -i stdin -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | tail -n +2 | \
    bedtools makewindows -b stdin -n 10 -i srcwinnum >> young.$f.TADDomain.extended.bed
done < $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADDomain.bed

# while read line 
# do
#     bedtools flank -i `echo "$line"` -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | \
#     cat - `echo "$line"` | \
#     bedtools makewindows -b stdin -n 10 -i srcwinnum | \
#     sort -k1,1 -k2,2n >> young.$f.TADDomain.extended.bed
# done < $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADDomain.bed

# bedtools flank -i $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADDomain.bed \
# -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | \
# cat - $prefix/08_HiCExplorer/compartmentOverlap/aged.$f.TADDomain.bed | \
# bedtools makewindows -b stdin -n 10 -i srcwinnum | \
# sort -k1,1 -k2,2n | uniq > aged.$f.TADDomain.extended.bed

# bedtools flank -i $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADDomain.bed \
# -g /nas/homes/benyang/Genome_References/sizes.mm10 -b 100000 | \
# cat - $prefix/08_HiCExplorer/compartmentOverlap/young.$f.TADDomain.bed | \
# bedtools makewindows -b stdin -n 10 -i srcwinnum | \
# sort -k1,1 -k2,2n | uniq > young.$f.TADDomain.extended.bed


done

# awk -v num_win=5 '{len = $3-$2; size = len/num_win; \
# for (start=$2; start<$3; start+=size) \
# print $1"\t"start"\t"start+size}' < test