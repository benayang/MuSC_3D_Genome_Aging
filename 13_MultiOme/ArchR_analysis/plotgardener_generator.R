library(optparse)

option_list = list(
    make_option("--gene", type="character", default=NULL,
                    help="Gene name", metavar="Ddit3"),
    make_option("--plot_range", type="character", default=NULL,
                    help="Range for plotgardener plot (chr_start_end)"),
    make_option("--viewpoint", type="character", default=NULL,
                    help="Viewpoint for subsetting Cicero connections object", metavar="character"),
    make_option("--res", type="integer", default=10000,
                    help="HiC matrix resolution", metavar="10000"),
    make_option("--hic_max", type="double", default=40,
                    help="Max HiC value", metavar="40"),
    make_option("--scATAC_max", type="double", default=NULL,
                help="Max scATAC value", metavar="double"),
    make_option("--ATAC_max", type="double", default=NULL,
                help="Max bulk ATAC value", metavar="double"),
    make_option("--H3K4me3_max", type="double", default=NULL,
                help="Max bulk H3K4me3 value", metavar="double"),
    make_option("--HiC", type="character", default=NULL,
                help="HiC file", metavar="character"),
    make_option("--scATAC", type="character", default=NULL,
                help="scATAC bigwig", metavar="character"),
    make_option("--bulkATAC", type="character", default=NULL,
                help="bulk ATAC bigwig", metavar="character"),
    make_option("--bulkH3K4me3", type="character", default=NULL,
                help="bulk H3K4me3 bigwig", metavar="character"),
    make_option("--conns", type="character", default=NULL,
                help="Connections GI file", metavar="character"),
    make_option("--TAD", type="character", default=NULL,
                help="TAD file", metavar="character"),
    make_option("--Loop", type="character", default=NULL,
                help="Loops file", metavar="character"),
    make_option("--Peaks", type="character", default=NULL,
                help="Peaks file", metavar="character"),
    make_option("--Peak2Gene", type="character", default=NULL,
                help="Peak2Gene file", metavar="character"),
    make_option("--GeneMatrix", type="character", default=NULL,
                help="GeneMatrix file", metavar="character"),
    make_option("--outfile", type="character", default=NULL,
                help="Output image file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(S4Vectors))
#suppressPackageStartupMessages(library(cicero))
#suppressPackageStartupMessages(library(ArchR))
#suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(plotgardener))
suppressPackageStartupMessages(library(GenomicInteractions))
suppressPackageStartupMessages(library(RColorBrewer))
#suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))

## Import TADs and loops ----------------------------------------
TADs = read.table(opt$TAD, sep='\t')
TAD_bed = TADs[, 1:3]
colnames(TAD_bed) = c("chrom","start","end")

loops = read.table(opt$Loop, sep='\t')
colnames(loops) = c("chrom1","start1","end1","chrom2","start2","end2","score")

## Get MACS2 peaks as data.frame ---------------------------------
peaks = readRDS(opt$Peaks)
peaks = cbind(chrom = seqnames(peaks), as.data.frame(ranges(peaks))) %>% dplyr::select(chrom,start,end,names)

## Read in required objects --------------------------------------
if (!is.null(opt$Peak2Gene)) {
    Peak2Gene = readRDS(opt$Peak2Gene)
} else {
    Peak2Gene = NULL
}

viewpoint = unlist(strsplit(opt$viewpoint, split="_"))
viewpoint = GRanges(seqnames=viewpoint[1], IRanges(start=as.numeric(viewpoint[2]), end=as.numeric(viewpoint[3])))

conns_GI = readRDS(opt$conns)

## Define plot region based on Peak-to-gene and Cicero connections --------------------
if(!is.null(opt$plot_range)) {
    plot_region = unlist(strsplit(opt$plot_range, split="_"))
} else {
    p2g_bound = c(min(start(anchorOne(Peak2Gene[which(Peak2Gene$gene==opt$gene), ]))),
                    max(end(anchorTwo(Peak2Gene[which(Peak2Gene$gene==opt$gene), ]))))
    p2g_region = GRanges(data.frame(seqnames=as.character(seqnames(viewpoint)),
                                        start = p2g_bound[1], 
                                        end = p2g_bound[2]))
    subset_conns = subsetByOverlaps(conns_GI, p2g_region, ignore.strand=T)
    subset_conns = subset_conns[which(subset_conns$coaccess>0.1), ]
    conns_bound = c(min(start(anchorOne(subset_conns))), max(end(anchorTwo(subset_conns))))

    plot_region = c(as.character(seqnames(viewpoint)),
                    min(p2g_bound[1], conns_bound[1]) - 5e4, 
                    max(p2g_bound[2], conns_bound[2]) + 5e4)
}

## read in bigwigs to define local ymin/ymax limits --------------------

params = pgParams(chrom=plot_region[1], chromstart=as.numeric(plot_region[2]), chromend=as.numeric(plot_region[3]), 
                    x=0.5, width=5, just=c("left","top"), default.inches="inches", assembly="mm10")

scATAC_bw = readBigwig(opt$scATAC, params=params)
bulkATAC_bw = readBigwig(opt$bulkATAC, params=params)
bulkH3K4me3_bw = readBigwig(opt$bulkH3K4me3, params=params)

assign_ylim = function(bw) {
    ylim <- c(0,quantile(bw$score, probs=c(0.999)))
    bw$score[bw$score < ylim[1]] <- ylim[1]
    bw$score[bw$score > ylim[2]] <- ylim[2]
    return(bw)
}

scATAC_bw = assign_ylim(scATAC_bw)
bulkATAC_bw = assign_ylim(bulkATAC_bw)
bulkH3K4me3_bw = assign_ylim(bulkH3K4me3_bw)

## Start assembing plotgardener plot ---------------------------------
print("Finished loading data, making plot now")

png(opt$outfile, res=300, units="in", width=6, height=10)
pageCreate(width=6, height=10, default.units="inches")

hic = plotHicTriangle(
    opt$HiC, 
    resolution=opt$res, height=1.75, 
    params=params, colorTrans="log", 
    palette = colorRampPalette(rev(brewer.pal(n = 11, "RdYlBu"))), 
    bg = brewer.pal(n=11,"RdYlBu")[11], 
    zrange=c(1,opt$hic_max), y=0.25
)

annoHeatmapLegend(
    plot = hic, x = 5.5, y = 0.25,
    ticks = T, fontcolor="black",
    width = 0.13, height = 1.2,
    just = c("right", "top")
)

TAD_plt = plotRanges(TAD_bed, params=params, fill="blue", y="b0", height=0.5)

loops_plt = plotPairsArches(loops, params=params, y='b0.1', height=0.75, fill="black", linecolor="black", archHeight="score", alpha=1)

bulk_H3K4me3 = plotSignal(bulkH3K4me3_bw, params=params, y='b0.1', height=0.75, linecolor="#661100", fill="#661100", 
                            range=c(0,ceiling(max(bulkH3K4me3_bw$score))), scale=T)

bulk_ATAC = plotSignal(bulkATAC_bw, params=params, y='b0.1', height=0.75, linecolor="#117733", fill="#117733", 
                        range=c(0,ceiling(max(bulkATAC_bw$score))), scale=T)

scATAC = plotSignal(scATAC_bw, params=params, y='b0.1', height=0.75, linecolor="#37a7db", fill="#37a7db", 
                    range=c(0,ceiling(max(scATAC_bw$score))), scale=T)

peaks_plt = plotRanges(peaks, params=params, fill="red1", y="b0.1", collapse=T, height=0.125)

if(!is.null(Peak2Gene)) {
    subset_p2g = Peak2Gene[which(Peak2Gene$gene == opt$gene), ]
    subset_p2g$h = pairdist(subset_p2g)
    subset_p2g$h = subset_p2g$h / max(subset_p2g$h)
} else {
    subset_p2g <- GRanges(data.frame(seqnames=plot_region[1], start=0, end=0, Correlation=0))
    subset_p2g$h = 0
}

p2g_plt = plotPairsArches(
    subset_p2g, 
    flip = T,
    lwd = 2,
    params=params, y='b0.1', height=1, 
    fill=colorby("Correlation", palette=colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black")),
                range = c(0.45, max(subset_p2g$Correlation))), 
    linecolor="fill", archHeight="h", alpha=0.6)

annoHeatmapLegend(
    plot = p2g_plt, fontcolor = "black",
    ticks = F, digits=2,
    params=params,
    x=as.numeric(gsub("inches","",p2g_plt$width,fixed=T)) + 0.5,
    y=as.numeric(gsub("inches","",p2g_plt$y,fixed=T)) + 0.125,
    width = 0.10, height = 0.75, fontsize = 10
)

conns_colnames = c("seqnames1","start1","end1","seqnames2","start2","end2","coaccess","distance")
subset_conns = as.data.frame(subset_conns)[,conns_colnames]

if(nrow(subset_conns)==0) {
    subset_conns = as.data.frame(conns_GI[which(mcols(conns_GI)$coaccess>0.1), ])[,conns_colnames]
    subset_conns$h = subset_conns$distance / max(subset_conns$distance)
    conns_plt = plotPairsArches(
        as.data.frame(conns_GI[which(mcols(conns_GI)$coaccess>0.1), ])[,conns_colnames], 
        params=params, 
        flip = T, lwd = 2,
        y='b0.1', height=1, 
        fill=colorby("coaccess", palette=colorRampPalette(c("dodgerblue2","firebrick2")), 
                    range=c(0.1, max(subset_conns$coaccess))), 
        archHeight="h", 
        linecolor="fill", 
        alpha=0.6)   
} else {
    print("subset_conns")
    subset_conns$h = subset_conns$distance / max(subset_conns$distance)
    conns_plt = plotPairsArches(
        subset_conns, 
        flip = T, lwd = 2,
        params=params, y='b0.1', height=1, 
        fill=colorby("coaccess", palette=colorRampPalette(c("dodgerblue2","firebrick2")), 
                    range=c(0.1, max(subset_conns$coaccess))), 
        archHeight="h", 
        linecolor="fill", 
        alpha=0.6)

    annoHeatmapLegend(
    plot = conns_plt, fontcolor = "black",
    ticks = F, digits=2,
    params=params,
    x=as.numeric(gsub("inches","",conns_plt$width,fixed=T)) + 0.5,
    y=as.numeric(gsub("inches","",conns_plt$y,fixed=T)) + 0.125,
    width = 0.10, height = 0.75, fontsize = 10)
}

genes_plt = plotGenes(
    params = params,
    geneHighlights = data.frame(
        "gene" = opt$gene,
        "color" = c("#225EA8")
    ),
    geneBackground = "grey",
    y = 'b0.1',
    fontsize=10
)

annoGenomeLabel(plot = genes_plt, fontsize=14, params=params, y='b0.1', scale = "Kb")

# plotGG(young_exp_plt, x=params$width+0.5, y=0.25, params=params, width=1, height=1.5)

pageGuideHide()
dev.off()
