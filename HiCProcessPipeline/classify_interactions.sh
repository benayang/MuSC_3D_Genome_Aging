#! /bin/bash

echo Running
total=$(awk '{sum+=$5}END{print sum}' "$1")
inter_total=$(awk '{if ($1 == $3) sum+=$5}END{print sum}' "$1")
awk -v t="$total" '{if ($1 != $3)sum+=$5} END {print "total inter-chromosomal: " sum " (" 100*sum/t "%)"}' "$1"
awk -v t="$total" '{if ($1 == $3)sum+=$5} END {print "total intra-chromosomal: " sum " (" 100*sum/t "%)"}' "$1"
awk -v t="$total" -v it="$inter_total" '{if ($1 == $3 && $4-$2<=20000)sum+=$5} END {print "short-range: " sum " (" 100*sum/t "% / " 100*sum/it "%)"}' "$1"
awk -v t="$total" -v it="$inter_total" '{if ($1 == $3 && $4-$2>20000)sum+=$5} END {print "long-range: " sum " (" 100*sum/t "% / " 100*sum/it "%)"}' "$1"


