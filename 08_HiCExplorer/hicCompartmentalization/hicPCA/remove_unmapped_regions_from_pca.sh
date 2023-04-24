for f in young_pca1 young_pca2 aged_pca1 aged_pca2
do

awk '{if($4 != -0.000000000000 && $4 != 0.000000000000) print $0}' ${f}.bedgraph > ${f}.filtered.bedgraph

done