prefix=/nas/homes/benyang/HiC

bedtools annotate -i $prefix/05_Arrowhead_HICCUPS/young_merged_loop_anchors.bed -files $prefix/HistoneTracks/atac.d0.Young.optimal.narrowPeak.gz $prefix/HistoneTracks/mm10LiftOver_GSM1148110_yQ_K4m3_peaks_w_input_control_peaks.bed -names ATAC H3K4me3 -counts | sort -k1,1 -k2,2n > $prefix/05_Arrowhead_HICCUPS/track_signal_at_anchors/young_peaks_per_young_loop_anchors.txt

bedtools annotate -i $prefix/05_Arrowhead_HICCUPS/aged_merged_loop_anchors.bed -files $prefix/HistoneTracks/atac.d0.Aged.optimal.narrowPeak.gz $prefix/HistoneTracks/mm10LiftOver_GSM1148118_oQSC_K4m3_MACSwInputControl.bed_peaks.bed -names ATAC H3K4me3 -counts | sort -k1,1 -k2,2n > $prefix/05_Arrowhead_HICCUPS/track_signal_at_anchors/aged_peaks_per_aged_loop_anchors.txt