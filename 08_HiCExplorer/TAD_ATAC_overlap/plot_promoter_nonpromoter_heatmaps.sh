mainDir=/nas/homes/benyang/HiC

young_promoter_boundary="$mainDir/08_HiCExplorer/TAD expression/young.merged.TAD.promoter.boundaries.knownGenes.bed"
aged_promoter_boundary="$mainDir/08_HiCExplorer/TAD expression/aged.merged.TAD.promoter.boundaries.knownGenes.bed"
young_non_promoter_boundary="$mainDir/08_HiCExplorer/TAD expression/young.merged.TAD.nonpromoter.boundaries.knownGenes.bed"
aged_non_promoter_boundary="$mainDir/08_HiCExplorer/TAD expression/aged.merged.TAD.nonpromoter.boundaries.knownGenes.bed"

young_ATAC=$mainDir/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw
aged_ATAC=$mainDir/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw

#plot ATAC over boundaries
# computeMatrix scale-regions --samplesLabel "Young" "Aged" \
# -p 50 \
# -S $young_ATAC $aged_ATAC \
# -R "$young_promoter_boundary" "$aged_promoter_boundary" "$young_non_promoter_boundary" "$aged_non_promoter_boundary" \
# -o $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.mat.gz \
# --upstream 20000 \
# --regionBodyLength 40000 \
# --downstream 20000 \
# --skipZeros \
# --binSize 100

# plotHeatmap -m $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.mat.gz \
# -out $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.png \
# --regionsLabel "Young Promoter" "Aged Promoter" "Young Non-Promoter" "Aged Non-Promoter" \
# --startLabel "Start" --endLabel "End" \
# --colorMap RdBu \
# --dpi 300 \
# --heatmapHeight 15 \
# --heatmapWidth 7

computeMatrix reference-point --samplesLabel "Young" "Aged" \
-p 50 \
-S $young_ATAC $aged_ATAC \
-R "$young_promoter_boundary" "$aged_promoter_boundary" "$young_non_promoter_boundary" "$aged_non_promoter_boundary" \
-o $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.refPoint.mat.gz \
--referencePoint center \
--upstream 40000 \
--downstream 40000 \
--skipZeros \
--binSize 100

plotHeatmap -m $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.refPoint.mat.gz \
-out $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.refPoint.png \
--regionsLabel "Young Promoter" "Aged Promoter" "Young Non-Promoter" "Aged Non-Promoter" \
--colorMap RdBu \
--dpi 300 \
--heatmapHeight 15 \
--heatmapWidth 7

# computeMatrix reference-point --samplesLabel "Young" "Aged" \
# -p 50 \
# -S $young_ATAC $aged_ATAC \
# -R "$young_promoter_boundary" "$aged_promoter_boundary" "$young_non_promoter_boundary" "$aged_non_promoter_boundary" \
# -o $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.TSS.mat.gz \
# --referencePoint TSS \
# --beforeRegionStartLength 3000 \
# --afterRegionStartLength 3000 \
# --skipZeros \
# --binSize 10

# plotHeatmap -m $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.TSS.mat.gz \
# -out $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.ATAC.knownGenes.TSS.perGroup.png \
# --regionsLabel "Young Promoter" "Aged Promoter" "Young Non-Promoter" "Aged Non-Promoter" \
# --perGroup \
# --plotType se \
# --colorMap RdBu \
# --dpi 300 \
# --heatmapHeight 15 \
# --heatmapWidth 4

# Expressed genes
# young_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"

# computeMatrix reference-point --samplesLabel "Young" "Aged" \
# -p 50 \
# -S $young_ATAC $aged_ATAC \
# -R "$young_A_promoter_boundary" "$young_B_promoter_boundary" \
# "$aged_A_promoter_boundary" "$aged_B_promoter_boundary" \
# "$young_A_promoter_nonBoundary" "$young_B_promoter_nonBoundary" \
# "$aged_A_promoter_nonBoundary" "$aged_B_promoter_nonBoundary" \
# -o $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.expressed.mat.gz \
# --referencePoint TSS \
# --beforeRegionStartLength 1000 \
# --afterRegionStartLength 3000 \
# --skipZeros \
# --binSize 10

# plotHeatmap -m $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.expressed.mat.gz \
# -out $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.expressed.perGroup.png \
# --regionsLabel "Young Boundary (A)" "Young Boundary(B)" \
# "Aged Boundary(A)" "Aged Boundary(B)" \
# "Young Non-Boundary (A)" "Young Non-Boundary(B)" \
# "Aged Non-Boundary(A)" "Aged Non-Boundary(B)" \
# --perGroup \
# --plotType se \
# --colorMap RdBu \
# --dpi 300 \
# --heatmapHeight 15 \
# --heatmapWidth 5

# young_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# young_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"
# aged_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.non_TADBoundary.promoter.knownGenes.expressedGenes.bed"

# all genes
# young_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.TADBoundary.promoter.knownGenes.bed"
# young_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.TADBoundary.promoter.knownGenes.bed"
# aged_A_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.TADBoundary.promoter.knownGenes.bed"
# aged_B_promoter_boundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.TADBoundary.promoter.knownGenes.bed"
# young_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.A.non_TADBoundary.promoter.knownGenes.bed"
# young_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/young.B.non_TADBoundary.promoter.knownGenes.bed"
# aged_A_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.A.non_TADBoundary.promoter.knownGenes.bed"
# aged_B_promoter_nonBoundary="$mainDir/08_HiCExplorer/compartmentOverlap/aged.B.non_TADBoundary.promoter.knownGenes.bed"

# computeMatrix reference-point --samplesLabel "Young" "Aged" \
# -p 50 \
# -S $young_ATAC $aged_ATAC \
# -R "$young_A_promoter_boundary" "$young_B_promoter_boundary" \
# "$aged_A_promoter_boundary" "$aged_B_promoter_boundary" \
# "$young_A_promoter_nonBoundary" "$young_B_promoter_nonBoundary" \
# "$aged_A_promoter_nonBoundary" "$aged_B_promoter_nonBoundary" \
# -o $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.mat.gz \
# --referencePoint TSS \
# --beforeRegionStartLength 1000 \
# --afterRegionStartLength 3000 \
# --skipZeros \
# --binSize 10

plotHeatmap -m $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.mat.gz \
-out $mainDir/08_HiCExplorer/TAD_ATAC_overlap/plot_heatmaps/ATAC.boundary.promoter.compartment.ATAC.knownGenes.TSS.perGroup.png \
--regionsLabel "Young Boundary (A)" "Young Boundary(B)" \
"Aged Boundary(A)" "Aged Boundary(B)" \
"Young Non-Boundary (A)" "Young Non-Boundary(B)" \
"Aged Non-Boundary(A)" "Aged Non-Boundary(B)" \
--perGroup \
--plotType se \
--colorMap RdBu \
--dpi 300 \
--heatmapHeight 15 \
--heatmapWidth 5