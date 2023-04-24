main='/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS'

cut -f4 $main/assign_genes_to_loops/young.merged.loop_anchor.knownGenes.bed | uniq > $main/HOMER/call_homer/young_anchor_genes.txt
cut -f4 $main/assign_genes_to_loops/aged.merged.loop_anchor.knownGenes.bed | uniq > $main/HOMER/call_homer/aged_anchor_genes.txt

young_genes=$main/HOMER/call_homer/young_anchor_genes.txt
aged_genes=$main/HOMER/call_homer/aged_anchor_genes.txt

young_outdir=$main/HOMER/call_homer/young_merged_loop_anchors_promoters
aged_outdir=$main/HOMER/call_homer/aged_merged_loop_anchors_promoters

findMotifs.pl $young_genes mouse $young_outdir -p 4 -start -400 -end 100 

findMotifs.pl $aged_genes mouse $aged_outdir -p 4 -start -400 -end 100 