for f in "A" "B" "A_to_B" "B_to_A" "StaticA" "StaticB"
do

echo -n "'ID'" $(cut -f4 aged.$f.TADBoundary.flank.bed) | tr " " "\n" | \
paste get_flank_coverage/aged_binned_scores_per_$f.tab - > get_flank_coverage/aged_binned_scores_per_${f}_ID.tab

echo -n "'ID'" $(cut -f4 young.$f.TADBoundary.flank.bed) | tr " " "\n" | \
paste get_flank_coverage/young_binned_scores_per_$f.tab - > get_flank_coverage/young_binned_scores_per_${f}_ID.tab

# echo -n "'ID'" $(cut -f4 ../aged.$f.TADBoundary.bed) | tr " " "\n" | \
# paste ../boundary_track_coverage/aged_scores_per_$f.tab - > ../boundary_track_coverage/aged_scores_per_${f}_ID.tab

# echo -n "'ID'" $(cut -f4 ../young.$f.TADBoundary.bed) | tr " " "\n" | \
# paste ../boundary_track_coverage/young_scores_per_$f.tab - > ../boundary_track_coverage/young_scores_per_${f}_ID.tab

done

# awk '{OFS="\t"} NR>1 {arr1[$8]+=$4; arr2[$8]+=$5; arr3[$8]+=$6; arr4[$8]+=$7; N[$8]++} END \
# { for (a in arr1) {avg1[a] = arr1[a]/N[a]} \
# for (a in arr2) {avg2[a] = arr2[a]/N[a]} \
# for (a in arr3) {avg3[a] = arr3[a]/N[a]} \
# for (a in arr4) {avg4[a] = arr4[a]/N[a]} \
# print $1,$2,$3,avg1}' get_flank_coverage/aged_binned_scores_per_A_ID.tab > get_averaged_flank/aged_binned_scores_per_A_ID_avg.tab