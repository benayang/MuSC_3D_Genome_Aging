#!/usr/bin/env bash

pjdir="/mnt/c/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/08_HiCExplorer/TAD expression"

for f in young.merged aged.merged
do
    cut -f 4 "$pjdir/$f.TAD.genebodies.bed" | sort | uniq > "$pjdir/$f.TAD.genenames.txt"
    cut -f 4 "$pjdir/$f.boundary.genebodies.bed" | sort | uniq > "$pjdir/$f.boundary.genenames.txt"
done
