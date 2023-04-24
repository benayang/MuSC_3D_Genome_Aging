sample="aged"

fancplot --width 6 \
-o ./04_FANC/tadtool/$sample.merged.withboundariesAndTAD.png \
chr1:10mb-30mb \
-p triangular ./02_HIC/$sample.merged.40kb.hic@40kb@KR -m 4000000 -vmax 50 \
-p layer ./09_OnTAD/$sample.merged.chr1.bed \
-p scores ./04_FANC/insulation/$sample.merged.insulation -vmin -1 -vmax 1 \
-p bar ./04_FANC/insulation/$sample.merged.boundaries_400kb.bed