prefix=/nas/homes/benyang/HiC

# Get 1kb promoter regions around TSS
bedtools slop -b 1000 -g /nas/homes/benyang/Genome_References/sizes.mm10 -i $prefix/get_tss/tss.gencode.vM25.basic.annotation.filtered.uniq.bed | \
grep -v "Gm" | grep -v "Rik" | uniq > tss.1kb.bed

######### Aged
# H3K4me3+/ATAC+
bedtools intersect -wa -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148118_oQSC_K4m3_MACSwInputControl.bed_peaks.bed | \
bedtools intersect -wa -a stdin -b $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak.gz | \
sort -k1,1 -k2,2n | uniq > aged.tss.1kb.posH3K4me3.posATAC.bed
# H3K4me3+/ATAC-
bedtools intersect -wa -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148118_oQSC_K4m3_MACSwInputControl.bed_peaks.bed | \
bedtools intersect -wa -v -a stdin -b $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak.gz | \
sort -k1,1 -k2,2n | uniq > aged.tss.1kb.posH3K4me3.negATAC.bed
# H3K4me3-
bedtools intersect -wa -v -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148118_oQSC_K4m3_MACSwInputControl.bed_peaks.bed | \
sort -k1,1 -k2,2n | uniq > aged.tss.1kb.negH3K4me3.bed

######### Young
# H3K4me3+/ATAC+
bedtools intersect -wa -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148110_yQ_K4m3_peaks_w_input_control_peaks.bed | \
bedtools intersect -wa -a stdin -b $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak.gz | \
sort -k1,1 -k2,2n | uniq > young.tss.1kb.posH3K4me3.posATAC.bed
# H3K4me3+/ATAC-
bedtools intersect -wa -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148110_yQ_K4m3_peaks_w_input_control_peaks.bed | \
bedtools intersect -wa -v -a stdin -b $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak.gz | \
sort -k1,1 -k2,2n | uniq > young.tss.1kb.posH3K4me3.negATAC.bed
# H3K4me3-
bedtools intersect -wa -v -a tss.1kb.bed -b $prefix/HistoneTracks/Tom_Rando/mm10LiftOver_GSM1148110_yQ_K4m3_peaks_w_input_control_peaks.bed | \
sort -k1,1 -k2,2n | uniq > young.tss.1kb.negH3K4me3.bed