#!/usr/bin/env bash

pjdir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/04_FANC/compartmentExpression/compartmentBed/100kb/tss"

#for f in young.A young.B aged.A aged.B
for f in A_to_B B_to_A static
do
	#cut -f 4 "$pjdir/$f.1kb.promoter.bed" | sort | uniq > "$pjdir/$f.1kb.promoter.genenames.txt"
	cut -f 4 "$pjdir/$f.tss.near.ab.bed" | sort | uniq > "$pjdir/$f.500bp.promoter.genenames.txt"
done

for f in aged young
do
	#cut -f 4 "$pjdir/$f.1kb.promoter.bed" | sort | uniq > "$pjdir/$f.1kb.promoter.genenames.txt"
	cut -f 4 "$pjdir/$f.tss.near.ab.bed" | sort | uniq > "$pjdir/$f.500bp.promoter.genenames.txt"
done

# grep -w -f $pjdir/Aged.genenames.txt $pjdir/Young.genenames.txt > $pjdir/common.genenames.txt
# grep -v -w -f $pjdir/Aged.genenames.txt $pjdir/Young.genenames.txt > $pjdir/unique.Young.genenames.txt 
# grep -v -w -f $pjdir/Young.genenames.txt $pjdir/Aged.genenames.txt > $pjdir/unique.Aged.genenames.txt
