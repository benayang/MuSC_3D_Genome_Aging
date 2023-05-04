#! /bin/bash

blacklist="/media/labuser/BenPassport/GenomeRefs/mm10.blacklist.bed"	

#python get_tss.py

bedtools intersect -v -a tss.gencode.vM25.basic.annotation.bed -b $blacklist > tss.gencode.vM25.basic.annotation.filtered.bed
	
cat tss.gencode.vM25.basic.annotation.filtered.bed | \
	awk -F "\t" 'BEGIN{OFS="\t";} {print $1,$2,$3;}' | \
	sort -k1,1 -k2,2n | uniq > tss.gencode.vM25.basic.annotation.filtered.nogenename.bed

sort -k1,1 -k2,2n tss.gencode.vM25.basic.annotation.filtered.bed | uniq > tss.gencode.vM25.basic.annotation.filtered.uniq.bed
