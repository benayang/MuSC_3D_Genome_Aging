prefix=/nas/homes/benyang/HiC

# for f in end start
# do

# bedtools intersect -wo -a $prefix/05_Arrowhead_HICCUPS/HOMER/aged_merged_loops_$f.bed -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > aged_${f}_CTCF.bed
# bedtools intersect -wo -a $prefix/05_Arrowhead_HICCUPS/HOMER/young_merged_loops_$f.bed -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed > young_${f}_CTCF.bed

# done

# create random intervals of the same number and average width of CTCF anchors
for f in end start
do
bedtools random -seed 6135 -g /nas/homes/benyang/Genome_References/sizes.mm10 \
-n `cat aged_${f}_CTCF.bed | wc -l` \
-l `awk '{sum+=($3-$2)}END{print sum/NR}' aged_${f}_CTCF.bed` > aged_${f}_random.bed

bedtools random -seed 6135 -g /nas/homes/benyang/Genome_References/sizes.mm10 \
-n `cat young_${f}_CTCF.bed | wc -l` \
-l `awk '{sum+=($3-$2)}END{print sum/NR}' young_${f}_CTCF.bed` > young_${f}_random.bed
done

# find CTCF occurrences in the random intervals
for f in end start
do
bedtools intersect -wo -a $prefix/05_Arrowhead_HICCUPS/HOMER/aged_${f}_random.bed -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed | \
sort -k1,1 -k2,2n | uniq > aged_${f}_CTCF_random.bed
bedtools intersect -wo -a $prefix/05_Arrowhead_HICCUPS/HOMER/young_${f}_random.bed -b $prefix/08_HiCExplorer/homer/ctcf.scan.out.bed | \
sort -k1,1 -k2,2n | uniq > young_${f}_CTCF_random.bed
done