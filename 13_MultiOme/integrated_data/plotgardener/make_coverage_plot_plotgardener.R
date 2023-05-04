library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
library(cicero)
library(org.Mm.eg.db)
library(plotgardener)
library(GenomicInteractions)
library(RColorBrewer)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")

projdir = "/nas/homes/benyang/HiC/13_MultiOme"

## Get Cicero connections
load(file.path(projdir, "integrated_data", "all_reps_Young_MuSC_cicero_conns_ccans.RData"))
young_conns = conns
young_ccans = ccans
load(file.path(projdir, "integrated_data", "all_reps_Aged_MuSC_cicero_conns_ccans.RData"))
aged_conns = conns
aged_ccans = ccans

rm(list=c("conns","ccans"))

# Convert conns to GenomicInteractions object
young_peak1_GR = StringToGRanges(young_conns$Peak1)
young_peak2_GR = StringToGRanges(young_conns$Peak2)
young_conns_GI = GenomicInteractions(young_peak1_GR, young_peak2_GR)
mcols(young_conns_GI)$coaccess = young_conns$coaccess

aged_peak1_GR = StringToGRanges(aged_conns$Peak1)
aged_peak2_GR = StringToGRanges(aged_conns$Peak2)
aged_conns_GI = GenomicInteractions(aged_peak1_GR, aged_peak2_GR)
mcols(aged_conns_GI)$coaccess = aged_conns$coaccess

## Import TADs
aged_TADs = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/aged.merged/40kb/aged.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed", sep="\t")
aged_TAD_bed = aged_TADs[, 1:3]
colnames(aged_TAD_bed) = c("chrom","start","end")
young_TADs = read.table("/nas/homes/benyang/HiC/08_HiCExplorer/young.merged/40kb/young.merged_min120kb_max1200kb_step40kb_thresh0.1_delta_0.01_fdr_domains.bed", sep="\t")
young_TAD_bed = young_TADs[, 1:3]
colnames(young_TAD_bed) = c("chrom","start","end")

## Get MACS2 MuSC peaks
aged_MuSC_peaks = readRDS("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_aged_MuSC_MACS2_peaks.RDS")
young_MuSC_peaks = readRDS("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_young_MuSC_MACS2_peaks.RDS")

## Import Loops
aged_loops = read.table("/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/aged_merged_hiccups/merged_loops_links.bedpe", sep="\t")
colnames(aged_loops) = c("chrom1","start1","end1","chrom2","start2","end2","score")
young_loops = read.table("/nas/homes/benyang/HiC/05_Arrowhead_HICCUPS/young_merged_hiccups/merged_loops_links.bedpe", sep="\t")
colnames(young_loops) = c("chrom1","start1","end1","chrom2","start2","end2","score")

## Get MuSC objects with LinkPeaks run
young_MuSC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/all_reps_young_MuSC_subset_linkPeaks.RDS")
aged_MuSC = readRDS("/nas/homes/benyang/HiC/13_MultiOme/all_reps_aged_MuSC_subset_linkPeaks.RDS")

young_MuSC_links = as.data.frame(Links(young_MuSC))
aged_MuSC_links = as.data.frame(Links(aged_MuSC))

## Convert LinkPeaks to bedpe files 
aged_links_bedpe = aged_MuSC_links %>% dplyr::select(seqnames, start, end, score)
aged_links_bedpe = aged_links_bedpe %>% dplyr::mutate(start1 = start, end1 = start+1, start2=end, end2=end+1)
aged_links_bedpe = cbind(aged_links_bedpe[,c("seqnames","start1","end1")], aged_links_bedpe[,c("seqnames","start2","end2","score")])
colnames(aged_links_bedpe)[1] = "chrom1"
colnames(aged_links_bedpe)[4] = "chrom2"
aged_links_bedpe_GI = GenomicInteractions(GRanges(setNames(aged_links_bedpe[,1:3], c("chrom","start","end"))),
                                            GRanges(setNames(aged_links_bedpe[,4:6], c("chrom","start","end"))),
                                            score = aged_links_bedpe[,7])

young_links_bedpe = young_MuSC_links %>% dplyr::select(seqnames, start, end, score)
young_links_bedpe = young_links_bedpe %>% dplyr::mutate(start1 = start, end1 = start+1, start2=end, end2=end+1)
young_links_bedpe = cbind(young_links_bedpe[,c("seqnames","start1","end1")], young_links_bedpe[,c("seqnames","start2","end2","score")])
colnames(young_links_bedpe)[1] = "chrom1"
colnames(young_links_bedpe)[4] = "chrom2"
young_links_bedpe_GI = GenomicInteractions(GRanges(setNames(young_links_bedpe[,1:3], c("chrom","start","end"))),
                                            GRanges(setNames(young_links_bedpe[,4:6], c("chrom","start","end"))),
                                            score = young_links_bedpe[,7])

## Get MACS2 peaks as data.frame
young_peaks = data.frame(strings = rownames(young_MuSC))
young_peaks = young_peaks %>% separate(strings, into=c("chrom","start","end"), remove=T)
young_peaks = young_peaks %>% mutate(start=as.integer(start), end=as.integer(end))

aged_peaks = data.frame(strings = rownames(aged_MuSC))
aged_peaks = aged_peaks %>% separate(strings, into=c("chrom","start","end"), remove=T)
aged_peaks = aged_peaks %>% mutate(start=as.integer(start), end=as.integer(end))

make_plotgardener_plot = function(gene_name, chrom, chromstart, chromend, subset_GR, res, hic_max, scATAC_max, ATAC_max, H3K4me3_max) {
    young_MuSC$Age = "Young"
    young_exp_plt = ExpressionPlot(young_MuSC, features=gene_name, group.by="Age", assay="RNA")
    aged_exp_plt = ExpressionPlot(aged_MuSC, features=gene_name, assay="RNA")

    ## Start assembing plotgardener plot
    params = pgParams(chrom=chrom, chromstart=chromstart, chromend=chromend, x=0.5, width=5, just=c("left","top"), default.inches="inches", assembly="mm10")

    png(file.path(projdir,"integrated_data", "plotgardener",sprintf("all_reps_young_%s_coveragePlt_plotgardener.png",gene_name)), res=300, units="in", width=7, height=10)
    pageCreate(width=7, height=10, default.units="inches")

    young_hic = plotHicTriangle("/nas/homes/benyang/HiC/02_HIC/young.merged/young.merged.hic", resolution=res, height=1.5, params=params, colorTrans="log", palette = colorRampPalette(rev(brewer.pal(n = 11, "RdYlBu"))), bg = brewer.pal(n=11,"RdYlBu")[11], zrange=c(1,hic_max), y=0.25)

    annoHeatmapLegend(
        plot = young_hic, x = 5.5, y = 0.25,
        ticks = T,
        width = 0.13, height = 1.2,
        just = c("right", "top")
    )

    young_sc_ATAC = plotSignal("/nas/homes/benyang/HiC/13_MultiOme/RawData/all_reps_combined_young_MuSC_atac_possorted.rpkm.bw", params=params, y='b0.1', height=0.75, linecolor="#37a7db", fill="#37a7db", range=c(0,scATAC_max), scale=T)

    young_bulk_ATAC = plotSignal("/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Young.bw", params=params, y='b0.1', height=0.75, linecolor="#7ecdbb", fill="#7ecdbb", range=c(0,ATAC_max), scale=T)

    young_bulk_H3K4me3 = plotSignal("/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Young.bw", params=params, y='b0.1', height=0.75, linecolor="black", fill="black", range=c(0,H3K4me3_max), scale=T)

    loops = plotPairsArches(young_loops, params=params, y='b0.1', height=0.75, fill="black", linecolor="black", archHeight="score", alpha=1)

    tmp_links = subsetByOverlaps(young_links_bedpe_GI, subset_GR, use.region="both")
    if(length(tmp_links)==0) { 
        young_links = plotPairsArches(young_links_bedpe_GI, params=params, y='b0.1', height=0.75, fill=colorby("score", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), linecolor="fill", archHeight="score", alpha=0)
    } else {
        young_links = plotPairsArches(tmp_links, params=params, y='b0.1', height=0.75, fill=colorby("score", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), linecolor="fill", archHeight="score", alpha=0.6)

        annoHeatmapLegend(
        plot = young_links, fontcolor = "black",
        ticks = T, digits=2,
        params=params, x=params$width+0.75, y=5.5,
        width = 0.10, height = 0.5, fontsize = 10)
    }

    tmp_conns = subsetByOverlaps(young_conns_GI, subset_GR, use.region="first")
    tmp_conns = tmp_conns[which(mcols(tmp_conns)$coaccess>0),]
    tmp_conns = as.data.frame(tmp_conns)[,c(1:3,6:8,12)]
    young_conns_plt = plotPairsArches(tmp_conns, params=params, y='b0.1', height=1, fill=colorby("coaccess", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), archHeight="coaccess", linecolor="fill", alpha=0.6)

    annoHeatmapLegend(
        plot = young_conns_plt, fontcolor = "black",
        ticks = T, digits=2,
        params=params, x=params$width+0.75, y=6.5,
        width = 0.10, height = 0.5, fontsize = 10
    )

    peaks_plt = plotRanges(young_MuSC_peaks, params=params, fill="blue", y="b0.1", height=0.5)

    TAD_plt = plotRanges(young_TAD_bed, params=params, fill="blue", y="b0.1", height=0.5)

    genes_plt = genesPlot <- plotGenes(
        params = params,
        geneHighlights = data.frame(
            "gene" = gene_name,
            "color" = c("#225EA8")
        ),
        geneBackground = "grey",
        y = 'b0.1',
        fontsize=10
    )

    annoGenomeLabel(plot = genes_plt, params=params, y='b0.1', scale = "Kb")

    plotGG(young_exp_plt, x=params$width+0.5, y=0.25, params=params, width=1, height=1.5)

    pageGuideHide()
    dev.off()

    ## Aged plot
    png(file.path(projdir,"integrated_data", "plotgardener",sprintf("all_reps_aged_%s_coveragePlt_plotgardener.png",gene_name)), res=300, units="in", width=7, height=10)
    pageCreate(width=7, height=10, default.units="inches")

    aged_hic = plotHicTriangle("/nas/homes/benyang/HiC/02_HIC/aged.merged/aged.merged.hic", resolution=res, height=1.5, params=params, colorTrans="log", palette = colorRampPalette(rev(brewer.pal(n = 11, "RdYlBu"))), bg = brewer.pal(n=11,"RdYlBu")[11], zrange=c(1,hic_max), y=0.25)

    annoHeatmapLegend(
        plot = aged_hic, x = 5.5, y = 0.25,
        ticks = T,
        width = 0.13, height = 1.2,
        just = c("right", "top")
    )

    aged_sc_ATAC = plotSignal("/nas/homes/benyang/HiC/13_MultiOme/RawData/all_reps_aged_MuSC_atac_possorted.rpkm.bw", params=params, y='b0.1', height=0.75, linecolor="#37a7db", fill="#37a7db", range=c(0,scATAC_max), scale=T)

    aged_bulk_ATAC = plotSignal("/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/atac.count.rpkm.track.nodup.Aged.bw", params=params, y='b0.1', height=0.75, linecolor="#7ecdbb", fill="#7ecdbb", range=c(0,ATAC_max), scale=T)

    aged_bulk_H3K4me3 = plotSignal("/nas/homes/benyang/HiC/HistoneTracks/get_signal_bigwigs/H3K4me3.count.rpkm.track.nodup.Aged.bw", params=params, y='b0.1', height=0.75, linecolor="black", fill="black", range=c(0,H3K4me3_max), scale=T)

    loops = plotPairsArches(aged_loops, params=params, y='b0.1', height=0.75, fill="black", linecolor="black", archHeight="score", alpha=1)

    tmp_links = subsetByOverlaps(aged_links_bedpe_GI, subset_GR, use.region="both")
    if(length(tmp_links)==0) { 
        aged_links = plotPairsArches(aged_links_bedpe_GI, params=params, y='b0.1', height=0.75, fill=colorby("score", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), linecolor="fill", archHeight="score", alpha=0)
    } else {
        aged_links = plotPairsArches(tmp_links, params=params, y='b0.1', height=0.75, fill=colorby("score", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), linecolor="fill", archHeight="score", alpha=0.6)

        annoHeatmapLegend(
        plot = aged_links, fontcolor = "black",
        ticks = T, digits=2,
        params=params, x=params$width+0.75, y=5.5,
        width = 0.10, height = 0.5, fontsize = 10
    )
    }

    tmp_conns = subsetByOverlaps(aged_conns_GI, subset_GR, use.region="first")
    tmp_conns = tmp_conns[which(mcols(tmp_conns)$coaccess>0),]
    tmp_conns = as.data.frame(tmp_conns)[,c(1:3,6:8,12)]
    aged_conns_plt = plotPairsArches(tmp_conns, params=params, y='b0.1', height=1, fill=colorby("coaccess", palette=colorRampPalette(c("dodgerblue2","firebrick2"))), archHeight="coaccess", linecolor="fill", alpha=0.6)

    annoHeatmapLegend(
        plot = aged_conns_plt, fontcolor = "black",
        ticks = T, digits=2,
        params=params, x=params$width+0.75, y=6.5,
        width = 0.10, height = 0.5, fontsize = 10
    )

    peaks_plt = plotRanges(aged_MuSC_peaks, params=params, fill="blue", y="b0.1", height=0.5)

    TAD_plt = plotRanges(aged_TAD_bed, params=params, fill="blue", y="b0.1", height=0.5)

    genes_plt = genesPlot <- plotGenes(
        params = params,
        geneHighlights = data.frame(
            "gene" = gene_name,
            "color" = c("#225EA8")
        ),
        geneBackground = "grey",
        y = 'b0.1',
        fontsize=10
    )

    annoGenomeLabel(plot = genes_plt, params=params, y='b0.1', scale = "Kb")

    plotGG(aged_exp_plt, x=params$width+0.5, y=0.25, params=params, width=1, height=1.5)

    pageGuideHide()
    dev.off()
}

### Myod1 -----------------------
# make_plotgardener_plot("Myod1", "chr7", 46023858, 46711000, GRanges(seqnames="chr7", IRanges(46375474,46377474)), 10000, 50, 5031, 659, 669)
# make_plotgardener_plot("Ddit3", "chr10", 126900000, 127800000, GRanges(seqnames="chr10", IRanges(127289793,127291793)), 10000, 40, 13965, 2525, 625)
make_plotgardener_plot("Ddit3", "chr10", 126900000, 127800000, GRanges(seqnames="chr10", IRanges(127289793,127291793)), 10000, 40, 32500, 3500, 1250)
