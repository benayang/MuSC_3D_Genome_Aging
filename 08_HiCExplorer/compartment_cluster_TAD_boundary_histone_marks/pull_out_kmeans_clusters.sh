for c in `seq 1 1 4`
do 

grep "cluster_$c" young_AB_kmeans4_sortedRegions.txt | sort -k1,1 -k2,2n | uniq > kmeans_clusters/young_AB_kmeans4_cluster$c.txt
grep "cluster_$c" aged_AB_kmeans4_sortedRegions.txt | sort -k1,1 -k2,2n | uniq > kmeans_clusters/aged_AB_kmeans4_cluster$c.txt

done