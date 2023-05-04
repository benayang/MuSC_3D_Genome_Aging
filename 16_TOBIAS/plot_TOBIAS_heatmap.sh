prefix='/nas/homes/benyang/HiC/16_TOBIAS'

# CEBPB_prefix=$prefix/merged_peakset_annotated/CEBPB_MA0466.3/beds/CEBPB_MA0466.3
# TOBIAS PlotHeatmap \
# --TFBS ${CEBPB_prefix}_young_merged_ATAC_bound.bed ${CEBPB_prefix}_young_merged_ATAC_unbound.bed \
# --TFBS ${CEBPB_prefix}_aged_merged_ATAC_bound.bed ${CEBPB_prefix}_aged_merged_ATAC_unbound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw $prefix/aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/CEBPB_MA0466.3_heatmap.png \
# --signal_labels Young Aged \
# --share_colorbar --sort_by -1

CTCF_prefix=$prefix/merged_peakset_annotated/CTCF_MA0139.1/beds/CTCF_MA0139.1
TOBIAS PlotHeatmap \
--TFBS ${CTCF_prefix}_young_merged_ATAC_bound.bed ${CTCF_prefix}_young_merged_ATAC_unbound.bed \
--TFBS ${CTCF_prefix}_aged_merged_ATAC_bound.bed ${CTCF_prefix}_aged_merged_ATAC_unbound.bed \
--signals $prefix/young_merged_ATAC_corrected.bw $prefix/aged_merged_ATAC_corrected.bw \
--output $prefix/merged_peakset_annotated/CTCF_MA0139.1_heatmap.png \
--signal_labels Young Aged \
--share_colorbar --sort_by -1

# TCFL5_prefix=$prefix/merged_peakset_annotated/TCFL5_MA0632.2/beds/TCFL5_MA0632.2
# TOBIAS PlotHeatmap \
# --TFBS ${TCFL5_prefix}_young_merged_ATAC_bound.bed ${TCFL5_prefix}_young_merged_ATAC_unbound.bed \
# --TFBS ${TCFL5_prefix}_aged_merged_ATAC_bound.bed ${TCFL5_prefix}_aged_merged_ATAC_unbound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw $prefix/aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/TCFL5_MA0632.2_heatmap.png \
# --signal_labels Young Aged \
# --share_colorbar --sort_by -1
