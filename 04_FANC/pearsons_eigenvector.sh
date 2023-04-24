#! /bin/bash

juicerJar="/media/labuser/BenPassport/Datasets/Hi-C Shallow Sequencing/Deep Sequencing/juicer_tools_1.22.01.jar"
norm="KR"
hic="/media/data/fastqs_4140-KS/aged_merged_paired_alignment_dedup.fixed.hic"
#hic="/media/labuser/BenPassport/Datasets/Hi-C Shallow Sequencing/Deep Sequencing/fastqs_4140-KS/02_HIC/4140-KS-1_CGATGTAT-AGATCTCG_S101/"
binsize="100000"
bpfrag="BP"
outPrefix="aged.merged"
outDir="/media/labuser/BenPassport/Datasets/Hi-C Shallow Sequencing/Deep Sequencing/fastqs_4140-KS/03_ABCompartments/$outPrefix/500kb"
#outDir="/media/labuser/BenPassport/Datasets/Hi-C Shallow Sequencing/Deep Sequencing/fastqs_4140-KS/03_ABCompartments/$outPrefix/50kb"

echo $outDir
mkdir -p "$outDir"

export _JAVA_OPTIONS="-Xms256M -Xmx20480M"

for c in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY
do
	echo $c '------------------'
	java -jar "$juicerJar" eigenvector -p $norm "$hic" $c $bpfrag $binsize "$outDir/${outPrefix}_${c}_KR_${binsize}${bpfrag}_eigen.txt"
	java -jar "$juicerJar" pearsons -p $norm "$hic" $c $bpfrag $binsize "$outDir/${outPrefix}_${c}_KR_${binsize}${bpfrag}_pearsons.txt"
done
