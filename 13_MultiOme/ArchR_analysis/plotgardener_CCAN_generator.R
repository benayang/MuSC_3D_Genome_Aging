library(optparse)

option_list = list(
    make_option("--gene", type="character", default=NULL,
                    help="Gene name", metavar="Ddit3"),
    make_option("--plot_range", type="character", default=NULL,
                    help="Range for plotgardener plot (chr_start_end)"),
    make_option("--viewpoint", type="character", default=NULL,
                    help="Viewpoint for subsetting Cicero connections object", metavar="character"),
    make_option("--aged_conns", type="character", default=NULL,
                help="Connections GI file", metavar="character"),
    make_option("--young_conns", type="character", default=NULL,
                help="Connections GI file", metavar="character"),
    make_option("--all_conns", type="character", default=NULL,
                help="Connections GI file", metavar="character"),
    make_option("--Peaks", type="character", default=NULL,
                help="Peaks file", metavar="character"),
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

## Get MACS2 peaks as data.frame ---------------------------------
peaks = readRDS(opt$Peaks)
peaks = cbind(chrom = seqnames(peaks), as.data.frame(ranges(peaks))) %>% dplyr::select(chrom,start,end,names)

## Read in required objects --------------------------------------
viewpoint = unlist(strsplit(opt$viewpoint, split="_"))
viewpoint = GRanges(seqnames=viewpoint[1], IRanges(start=as.numeric(viewpoint[2]), end=as.numeric(viewpoint[3])))

young_conns_gi = readRDS(opt$young_conns)
aged_conns_gi = readRDS(opt$aged_conns)
all_conns_gi = readRDS(opt$all_conns)


## Start assembing plotgardener plot ---------------------------------

plot_region = unlist(strsplit(opt$plot_range, sep="_"))

params = pgParams(chrom=plot_region[1], chromstart=as.numeric(plot_region[2]), chromend=as.numeric(plot_region[3]), 
                    x=0.5, width=5, just=c("left","top"), default.inches="inches", assembly="mm10")

print("Finished loading data, making plot now")

png(opt$outfile, res=300, units="in", width=6, height=10)
pageCreate(width=6, height=10, default.units="inches")

conns_colnames = c("seqnames1","start1","end1","seqnames2","start2","end2","coaccess","distance")
young_subset_conns = subsetByOverlaps(young_conns_gi, )
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
