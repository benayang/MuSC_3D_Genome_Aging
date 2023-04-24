#chrom=($(seq 1 19 1))
#chrom+=("X")
hic="young.merged"
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X
do
    echo $c
    java -jar "/home/benjy/c/Users/benjy/Documents/juicer_tools_1.22.01.jar" dump observed KR \
    "/home/benjy/c/Users/benjy/Downloads/${hic}.hic" "chr${c}" "chr${c}" BP 40000 \
    "/home/benjy/c/Users/benjy/Downloads/dumpedMats/${hic}_obs_KR_chr${c}_40kb.mat"
done