prefix='/nas/homes/benyang/HiC'

hicPCA -m $prefix/08_HiCExplorer/aged.merged_100kb_KR.cool \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX \
--whichEigenvectors 1 2 \
--outputFileName aged_pca1.bedgraph aged_pca2.bedgraph \
--format bedgraph \
--extraTrack "/nas/homes/benyang/HiC/08_HiCExplorer/hicCompartmentalization/hicPCA/mm10_Gencode_vM25_genetrack.sorted.bed"

hicPCA -m $prefix/08_HiCExplorer/young.merged_100kb_KR.cool \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX \
--whichEigenvectors 1 2 \
--outputFileName young_pca1.bedgraph young_pca2.bedgraph \
--format bedgraph \
--extraTrack "/nas/homes/benyang/HiC/08_HiCExplorer/hicCompartmentalization/hicPCA/mm10_Gencode_vM25_genetrack.sorted.bed"