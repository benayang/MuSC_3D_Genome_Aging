prefix=/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR

# Rscript plotgardener_generator.R \
# --gene "Myod1" \
# --plot_range "chr7_46023858_46711000" \
# --viewpoint "chr7_46375474_46377474" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 5031 \
# --ATAC_max 659 \
# --H3K4me3_max 669 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/young_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Myod1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Myod1" \
# --plot_range "chr7_46023858_46711000" \
# --viewpoint "chr7_46375474_46377474" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 5031 \
# --ATAC_max 659 \
# --H3K4me3_max 669 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/aged_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Myod1_coverage.png

# --viewpoint "chr10_127289793_127291793" \
# 127290774 127296288

# Rscript plotgardener_generator.R \
# --gene "Ddit3" \
# --viewpoint "chr10_127290774_127296288" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Ddit3_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Ddit3" \
# --viewpoint "chr10_127290774_127296288" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Ddit3_coverage.png

# Kdm5a
# chr6:120364099-120444574 (+)
# id = NM_145997.2

# Rscript plotgardener_generator.R \
# --gene "Kdm5a" \
# --viewpoint "chr6_120364099_120444574" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Kdm5a_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Kdm5a" \
# --viewpoint "chr6_120364099_120444574" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Kdm5a_coverage.png

# Ndufa4l2
# chr10:127514939-127517154 (+)
# id = NM_001098789.1

# Rscript plotgardener_generator.R \
# --gene "Ndufa4l2" \
# --viewpoint "chr10_127514939_127517154" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Ndufa4l2_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Ndufa4l2" \
# --viewpoint "chr10_127514939_127517154" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Ndufa4l2_coverage.png

# Trp53bp2
# chr1:182409167-182462436 (+)
# id = NM_173378.2

# Rscript plotgardener_generator.R \
# --gene "Trp53bp2" \
# --viewpoint "chr1_182409167_182462436" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Trp53bp2_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Trp53bp2" \
# --viewpoint "chr1_182409167_182462436" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Trp53bp2_coverage.png

# Jund
# chr8:70697739-70700616 (+)
# id = NM_001286944.1

Rscript plotgardener_generator.R \
--gene "Jund" \
--viewpoint "chr8_70697739_70700616" \
--res 10000 \
--hic_max 40 \
--scATAC_max 10 \
--ATAC_max 3500 \
--H3K4me3_max 1250 \
--HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
--scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
--bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
--bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
--Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
--conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
--TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
--Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
--Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
--GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
--outfile $prefix/plotgardener/young_Jund_coverage.png

Rscript plotgardener_generator.R \
--gene "Jund" \
--viewpoint "chr8_70697739_70700616" \
--res 10000 \
--hic_max 40 \
--scATAC_max 10 \
--ATAC_max 3500 \
--H3K4me3_max 1250 \
--HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
--scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
--bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
--bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
--Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
--conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
--TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
--Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
--Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
--GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
--outfile $prefix/plotgardener/aged_Jund_coverage.png


# Ctr9
# chr7:111028951-111056377 (+)
# id = NM_009431.2

# Rscript plotgardener_generator.R \
# --gene "Ctr9" \
# --viewpoint "chr7_111028951_111056377" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Ctr9_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Ctr9" \
# --viewpoint "chr7_111028951_111056377" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Ctr9_coverage.png

# Max
# chr12:76937269-76962248 (-)
# id = NM_008558.2

# Rscript plotgardener_generator.R \
# --gene "Max" \
# --viewpoint "chr12_76937269_76962248" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Max_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Max" \
# --viewpoint "chr12_76937269_76962248" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Max_coverage.png

# Smarcc1
# chr9:110132024-110240178 (+)
# id = NM_009211.2

# Rscript plotgardener_generator.R \
# --gene "Smarcc1" \
# --viewpoint "chr9_110132024_110240178" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Smarcc1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Smarcc1" \
# --viewpoint "chr9_110132024_110240178" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Smarcc1_coverage.png

# Kdm5b
# chr1:134560178-134632878 (+)
# id = NM_152895.2

# Rscript plotgardener_generator.R \
# --gene "Kdm5b" \
# --viewpoint "chr1_134560178_134632878" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Kdm5b_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Kdm5b" \
# --viewpoint "chr1_134560178_134632878" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Kdm5b_coverage.png

# Gadd45gip1
# chr8:84832282-84835482 (+)
# id = NM_183358.4

# Rscript plotgardener_generator.R \
# --gene "Gadd45gip1" \
# --viewpoint "chr8_84832282_84835482" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Gadd45gip1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Gadd45gip1" \
# --viewpoint "chr8_84832282_84835482" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Gadd45gip1_coverage.png

# Atp1a2
# chr1:172271709-172298064 (-)
# id = NM_178405.3

# Rscript plotgardener_generator.R \
# --gene "Atp1a2" \
# --viewpoint "chr1_172271709_172298064" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Atp1a2_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Atp1a2" \
# --viewpoint "chr1_172271709_172298064" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Atp1a2_coverage.png

# Myc
# chr15:61985341-61990361 (+)
# id = NM_001177354.1
# --------------
# Exon number: 2
# Amino acid coding number: 204
# chr15:61987510-61988278

# Rscript plotgardener_generator.R \
# --gene "Myc" \
# --viewpoint "chr15_61987510_61988278" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Myc_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Myc" \
# --viewpoint "chr15_61987510_61988278" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Myc_coverage.png

# Ndufc1
# chr3:51405479-51408955 (-)
# id = NM_025523.1

# Rscript plotgardener_generator.R \
# --gene "Ndufc1" \
# --viewpoint "chr3_51405479_51408955" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Ndufc1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Ndufc1" \
# --viewpoint "chr3_51405479_51408955" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Ndufc1_coverage.png

# Kdm2a
# chr19:4316144-4397077 (-)
# id = NM_001001984.2
# --------------
# Exon number: 7
# Amino acid coding number: 189
# chr19:4352268-4352374

# Rscript plotgardener_generator.R \
# --gene "Kdm2a" \
# --viewpoint "chr19_4316144_4397077" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Kdm2a_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Kdm2a" \
# --viewpoint "chr19_4316144_4397077" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Kdm2a_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Trp53" \
# --plot_range "chr11_69000000_71000000" \
# --viewpoint "chr11_69580359_69591873" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/young_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Trp53_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Trp53" \
# --plot_range "chr11_69000000_71000000" \
# --viewpoint "chr11_69580359_69591873" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/aged_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Trp53_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Jund" \
# --plot_range "chr8_70000000_71250000" \
# --viewpoint "chr8_70697739_70700616" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/young_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Jund_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Jund" \
# --plot_range "chr8_70000000_71250000" \
# --viewpoint "chr8_70697739_70700616" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/aged_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Jund_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Samd1" \
# --plot_range "chr8_83787095_84210609" \
# --viewpoint "chr8_83997672_84000386 " \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/young_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Samd1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Samd1" \
# --plot_range "chr8_83787095_84210609" \
# --viewpoint "chr8_83997672_84000386 " \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/aged_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Samd1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Sp1" \
# --plot_range "chr15_101000000_103000000" \
# --viewpoint "chr15_102406316_102436404" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/young_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Sp1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Sp1" \
# --plot_range "chr15_101000000_103000000" \
# --viewpoint "chr15_102406316_102436404" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --Peak2Gene $prefix/plotgardener/aged_links_GI.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Sp1_coverage.png

# chr11:100,253,022-100,800,076
# Acly
# chr11:100476352-100528000 (-)
# id = NM_001199296.1

# Rscript plotgardener_generator.R \
# --gene "Acly" \
# --plot_range "chr11_100253022_100800076" \
# --viewpoint "chr11_100476352_100528000" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Acly_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Acly" \
# --plot_range "chr11_100253022_100800076" \
# --viewpoint "chr11_100476352_100528000" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Acly_coverage.png

# chr16:29,986,154-30,146,864
# Hes1
# chr16:30065357-30067796 (+)
# id = NM_008235.2

#--plot_range "chr16_29744882_30387725" \
# Rscript plotgardener_generator.R \
# --gene "Hes1" \
# --viewpoint "chr16_30065357_30067796" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Young-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/young_MuSC/young_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Hes1_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Hes1" \
# --viewpoint "chr16_30065357_30067796" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 10 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/all_MuSC/GroupBigWigs/Age/Aged-TileSize-100-normMethod-ReadsInTSS-ArchR.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --Peak2Gene $prefix/plotgardener/MuSC_p2g_bedpe.RDS \
# --conns $prefix/aged_MuSC/aged_conns_gi_peakmatrix.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Hes1_coverage.png

# chr17:33,879,339-34,190,220
# Rxrb
# chr17:34031812-34038403 (+)
# id = NM_011306.4

# Rscript plotgardener_generator.R \
# --gene "Rxrb" \
# --plot_range "chr17_33879339_34190220" \
# --viewpoint "chr17_34031812_34038403" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Rxrb_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Rxrb" \
# --plot_range "chr17_33879339_34190220" \
# --viewpoint "chr17_34031812_34038403" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Rxrb_coverage.png

# chr11:100,418,994-101,218,116
# Stat5b
# chr11:100780731-100850585 (-)
# id = NM_001113563.1

# Rscript plotgardener_generator.R \
# --gene "Stat5b" \
# --plot_range "chr11_100418994_101218116" \
# --viewpoint "chr11_100780731_100850585" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Stat5b_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Stat5b" \
# --plot_range "chr11_100418994_101218116" \
# --viewpoint "chr11_100780731_100850585" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Stat5b_coverage.png

# chr1:78,361,827-79,001,582
# Acsl3
# chr1:78657825-78707743 (+)
# id = NM_001136222.1

# Rscript plotgardener_generator.R \
# --gene "Acsl3" \
# --plot_range "chr1_78361827_79001582" \
# --viewpoint "chr1_78657825_78707743" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic" \
# --scATAC $prefix/bigwigs/combined_young_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw" \
# --conns $prefix/plotgardener/young_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/young_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/young_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/young_Acsl3_coverage.png

# Rscript plotgardener_generator.R \
# --gene "Acsl3" \
# --plot_range "chr1_78361827_79001582" \
# --viewpoint "chr1_78657825_78707743" \
# --res 10000 \
# --hic_max 40 \
# --scATAC_max 32500 \
# --ATAC_max 3500 \
# --H3K4me3_max 1250 \
# --HiC "/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic" \
# --scATAC $prefix/bigwigs/aged_MuSC_atac_possorted.rpkm.bw \
# --bulkATAC "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw" \
# --bulkH3K4me3 "/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw" \
# --conns $prefix/plotgardener/aged_conns_GI.RDS \
# --TAD "/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed" \
# --Loop "/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe" \
# --Peaks $prefix/plotgardener/aged_MuSC_peaks.RDS \
# --GeneMatrix $prefix/plotgardener/aged_MuSC_genematrix.RDS \
# --outfile $prefix/plotgardener/aged_Acsl3_coverage.png