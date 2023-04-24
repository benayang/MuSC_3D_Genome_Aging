library(plotgardener)
library(ggplot2)
library(tidyr)
library(dplyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

DLR = read.table("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\15_HOMER_HIC\\AgedVsYoung.DLR.merged.bedgraph",
                 sep="\t", header=F, skip=1)
colnames(DLR) = c("chrom","start","end","score")
DLR_ranges = GRanges(DLR)

ICF = read.table("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\15_HOMER_HIC\\AgedVsYoung.ICF.merged.bedgraph",
                 sep="\t", header=F, skip=1)
colnames(ICF) = c("chrom","start","end","score")
ICF_ranges = GRanges(ICF)

params = pgParams(assembly="mm10", chrom="chr1", chromstart=78144036, chromend=81380134, x=0.5,
                  just = c("left", "top"), default.units = "inches", resolution=10000)

pageCreate(width = 5, height = 6, default.units = "inches")

youngHicPlot <- plotHicSquare(
  data = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC\\young.merged.hic", 
  params = params, colorTrans = "log2",
  half="top",
  zrange = c(0.0001,75),
  matrix = "observed
  ",
  y = 0.5, width = 3, height = 3
)
agedHicPlot <- plotHicSquare(
  data = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC\\aged.merged.hic", 
  params = params, colorTrans = "log2",
  half="bottom",
  zrange = c(0.0001,50),
  matrix = "log2oe", palette = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu"))),
  y = 0.5, width = 3, height = 3
)

annoHeatmapLegend(
  plot = hicPlot, x = 3.75, y = 0.5,
  width = 0.12, height = 1.2,
  just = c("left", "top"), default.units = "inches"
)

#DLR_ranges[unique(findOverlaps(DLR_ranges, type = "any", select = "first"))]
DLRPlot = plotSignal(
  data = DLR_ranges,
  params = params, negData = T,
  y = 'b0.1', fill="#37a7db",
  width=3, height=0.5 
)

ICFPlot = plotSignal(
  data = ICF_ranges,
  params = params, negData = T,
  y = 'b0.1', linecolor="red", fill="red",
  width=3, height=0.5 
)
