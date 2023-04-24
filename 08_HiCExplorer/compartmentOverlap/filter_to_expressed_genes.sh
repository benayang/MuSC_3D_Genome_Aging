for f in `ls | grep promoter.knownGenes.bed`
do

fname=`basename -s .bed $f`
grep -wf expressed_genes_at_boundaries.txt $f > $fname.expressedGenes.bed

done