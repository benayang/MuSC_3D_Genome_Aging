prefix='/nas/homes/benyang/HiC'

bedtools slop -b 1000 -i $prefix/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.knownGenes.bed -g /nas/homes/benyang/Genome_References/sizes.mm10 > tss_1kb.bed

for x in Shared Shifted Split Merged Indeterminate
do

    grep -w $x $prefix/08_HiCExplorer/young.TAD.domain.classified.bed | bedtools intersect -a ./tss_1kb.bed -b stdin | uniq > ${x}_TAD_genes.txt

done