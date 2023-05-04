prefix='/nas/homes/benyang/HiC/16_TOBIAS'

# TCFL5_prefix=$prefix/merged_peakset_annotated/TCFL5_MA0632.2/beds/TCFL5_MA0632.2
# TOBIAS PlotAggregate \
# --TFBS ${TCFL5_prefix}_young_merged_ATAC_bound.bed ${TCFL5_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/TCFL5_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 75 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${TCFL5_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/TCFL5_aggregate.png \
# --share_y both --plot_boundaries --flank 75 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# CTCF_prefix=$prefix/merged_peakset_annotated/CTCF_MA0139.1/beds/CTCF_MA0139.1
# TOBIAS PlotAggregate \
# --TFBS ${CTCF_prefix}_young_merged_ATAC_bound.bed ${CTCF_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/CTCF_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${CTCF_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/CTCF_aggregate.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# Zic2_prefix=$prefix/merged_peakset_annotated/Zic2_MA1629.1/beds/Zic2_MA1629.1
# TOBIAS PlotAggregate \
# --TFBS ${Zic2_prefix}_young_merged_ATAC_bound.bed ${Zic2_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/Zic2_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# CEBPB_prefix=$prefix/merged_peakset_annotated/CEBPB_MA0466.3/beds/CEBPB_MA0466.3
# TOBIAS PlotAggregate \
# --TFBS ${CEBPB_prefix}_young_merged_ATAC_bound.bed ${CEBPB_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/CEBPB_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${CEBPB_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/CEBPB_aggregate.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# KLF15_prefix=$prefix/merged_peakset_annotated/KLF15_MA1513.1/beds/KLF15_MA1513.1
# TOBIAS PlotAggregate \
# --TFBS ${KLF15_prefix}_young_merged_ATAC_bound.bed ${KLF15_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/KLF15_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${KLF15_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/KLF15_aggregate.png \
# --share_y both --plot_boundaries --flank 60 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# KLF7_prefix=$prefix/merged_peakset_annotated/KLF7_MA1959.1/beds/KLF7_MA1959.1
# TOBIAS PlotAggregate \
# --TFBS ${KLF7_prefix}_young_merged_ATAC_bound.bed ${KLF7_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated_expressed_TF/KLF7_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 50 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${KLF7_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated_expressed_TF/KLF7_aggregate.png \
# --share_y both --plot_boundaries --flank 50 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# KLF10_prefix=$prefix/merged_peakset_annotated/KLF10_MA1511.2/beds/KLF10_MA1511.2
# TOBIAS PlotAggregate \
# --TFBS ${KLF10_prefix}_young_merged_ATAC_bound.bed ${KLF10_prefix}_aged_merged_ATAC_bound.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated_expressed_TF/KLF10_aggregate_bound.png \
# --share_y both --plot_boundaries --flank 50 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --TFBS-labels Young Aged \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# TOBIAS PlotAggregate \
# --TFBS ${KLF10_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated_expressed_TF/KLF10_aggregate.png \
# --share_y both --plot_boundaries --flank 50 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# ZNF610_prefix=$prefix/merged_peakset_annotated/ZNF610_MA1713.1/beds/ZNF610_MA1713.1
# TOBIAS PlotAggregate \
# --TFBS ${ZNF610_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/ZNF610_aggregate.png \
# --share_y both --plot_boundaries --flank 100 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

# HES1_prefix=$prefix/merged_peakset_annotated/HES1_MA1099.2/beds/HES1_MA1099.2
# TOBIAS PlotAggregate \
# --TFBS ${HES1_prefix}_all.bed \
# --signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
# --output $prefix/merged_peakset_annotated/HES1_aggregate.png \
# --share_y both --plot_boundaries --flank 50 --log-transform \
# --signal-labels "Young ATAC" "Aged ATAC" \
# --blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 

NR1D1_prefix=$prefix/merged_peakset_annotated_expressed_TF/NR1D1_MA1531.1/beds/NR1D1_MA1531.1
TOBIAS PlotAggregate \
--TFBS ${NR1D1_prefix}_all.bed \
--signals $prefix/young_merged_ATAC_corrected.bw aged_merged_ATAC_corrected.bw \
--output $prefix/merged_peakset_annotated_expressed_TF/NR1D1_aggregate.png \
--share_y both --plot_boundaries --flank 50 --log-transform \
--signal-labels "Young ATAC" "Aged ATAC" \
--blacklist "/nas/homes/benyang/Genome_References/mm10-blacklist.v2.bed" 