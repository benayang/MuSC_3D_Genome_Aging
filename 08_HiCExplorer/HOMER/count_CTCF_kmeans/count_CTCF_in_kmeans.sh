prefix=/nas/homes/benyang/HiC/08_HiCExplorer

for c in `seq 1 1 6`
do

bedtools intersect -c -a $prefix/kmeans_TAD_boundaries/aged_cluster${c}.bedgraph -b $prefix/homer/ctcf.scan.out.bed > aged_cluster${c}_CTCF.bed
bedtools intersect -c -a $prefix/kmeans_TAD_boundaries/young_cluster${c}.bedgraph -b $prefix/homer/ctcf.scan.out.bed > young_cluster${c}_CTCF.bed

bedtools intersect -wo -a $prefix/kmeans_TAD_boundaries/aged_cluster${c}.bedgraph -b $prefix/homer/ctcf.scan.out.bed > aged_cluster${c}_CTCF_strand.bed
bedtools intersect -wo -a $prefix/kmeans_TAD_boundaries/young_cluster${c}.bedgraph -b $prefix/homer/ctcf.scan.out.bed > young_cluster${c}_CTCF_strand.bed

done