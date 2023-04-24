library(ggplot2)
library(dplyr)
library(tidyr)
library(plotgardener)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\Virtual4C"

young.domain.bed = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                        "young.TAD.domain.classified.bed"))
aged.domain.bed = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                       "aged.TAD.domain.classified.bed"))


# Define gene parameters --------------------------------------------------

gene = "Mpc2"
chrom="chr1"

aged = as.data.frame(t(read.table(file.path(projdir,paste0(gene,"_aged_data.txt")), sep="\t", header=F, row.names=c("start","score"), colClasses = c("numeric","numeric"))))
young = as.data.frame(t(read.table(file.path(projdir,paste0(gene,"_young_data.txt")), sep="\t", header=F, row.names=c("start","score"), colClasses = c("numeric","numeric"))))

young_bed = data.frame(chr=chrom, 
                       start = seq(min(as.integer(young$start)),max(as.integer(young$start)),5000))
young_bed$end = young_bed$start + 4999
young_bed = left_join(young_bed, young, by=c("start"))
#young_bed$smooth_scores = smooth.spline(young_bed$score)

aged_bed = data.frame(chr=chrom, 
                      start = seq(min(as.integer(aged$start)),max(as.integer(aged$start)),5000))
aged_bed$end = aged_bed$start + 4999
aged_bed = left_join(aged_bed, aged, by=c("start"))
#aged_bed$smooth_scores = sgolayfilt(aged_bed$score, p=3)

# Restrict plot to TADs surrounding genes ---------------------------------
# Pax7
# gene_start = 139737059; gene_end = 139833026
# young_bed = young_bed %>% dplyr::filter(start>138868292 & end<140531190)
# aged_bed = aged_bed %>% dplyr::filter(start>138868292 & end<140531190)
# Mef2a
# gene_start = 67231163; gene_end = 67372858;
# young_bed = young_bed %>% dplyr::filter(start>66668231 & end<68374112)
# aged_bed = aged_bed %>% dplyr::filter(start>66668231 & end<68374112)
# Kdm5b
# gene_start = 134560178; gene_end = 134632878;
# young_bed = young_bed %>% dplyr::filter(start>134113889 & end<135443577)
# aged_bed = aged_bed %>% dplyr::filter(start>134113889 & end<135443577)
# Mtor
# gene_start = 148448582; gene_end = 148557685;
# young_bed = young_bed %>% dplyr::filter(start>147824771 & end<149729874)
# aged_bed = aged_bed %>% dplyr::filter(start>147824771 & end<149729874)
# Foxo3
# gene_start = 42185786; gene_end = 42276742;
# young_bed = young_bed %>% dplyr::filter(start>40991533 & end<42769626)
# aged_bed = aged_bed %>% dplyr::filter(start>40991533 & end<42769626)
# Sesn2
# gene_start = 132492804; gene_end = 132510456;
# young_bed = young_bed %>% dplyr::filter(start>131831552 & end<133564320)
# aged_bed = aged_bed %>% dplyr::filter(start>131831552 & end<133564320)
# Mxi1
# gene_start = 53310506; gene_end = 53375810;
# young_bed = young_bed %>% dplyr::filter(start>50127905 & end<54067414)
# aged_bed = aged_bed %>% dplyr::filter(start>50127905 & end<54067414)
# Ddit3
# gene_start = 127288793; gene_end = 127298288;
# young_bed = young_bed %>% dplyr::filter(start>126873091 & end<127501575)
# aged_bed = aged_bed %>% dplyr::filter(start>126873091 & end<127501575)
# Spry2
# gene_start = 105891947; gene_end = 105896819;
# young_bed = young_bed %>% dplyr::filter(start>105201592 & end<106403934)
# aged_bed = aged_bed %>% dplyr::filter(start>105201592 & end<106403934)
# Stat1
# gene_start = 52119438; gene_end = 52161865;
# young_bed = young_bed %>% dplyr::filter(start>51699207 & end<52716926)
# aged_bed = aged_bed %>% dplyr::filter(start>51699207 & end<52716926)
# Kcnn2
# gene_start = 45268860; gene_end = 45685883;
# young_bed = young_bed %>% dplyr::filter(start>44560504 & end<46449466)
# aged_bed = aged_bed %>% dplyr::filter(start>44560504 & end<46449466)
# Mpc1
# gene_start = 8283813; gene_end = 8297661;
# young_bed = young_bed %>% dplyr::filter(start>8166577 & end<8413850)
# aged_bed = aged_bed %>% dplyr::filter(start>8166577 & end<8413850)
# Mpc2
gene_start = 165461208; gene_end = 165481214;
young_bed = young_bed %>% dplyr::filter(start>165216530 & end<165725389)
aged_bed = aged_bed %>% dplyr::filter(start>165216530 & end<165725389)


# Make plot ---------------------------------------------------------------

young_zero_idx = which(young_bed$score==0)
aged_zero_idx = which(aged_bed$score==0)

young_bed[young_zero_idx, "score"] = (young_bed[young_zero_idx-1, "score"] + young_bed[young_zero_idx+1, "score"])/2
aged_bed[aged_zero_idx, "score"] = (aged_bed[aged_zero_idx-1, "score"] + aged_bed[aged_zero_idx+1, "score"])/2
max_signal = max(c(young_bed$score, aged_bed$score), na.rm=T) + 1

params <- pgParams(chrom = chrom, chromstart = min(aged_bed$start), chromend = max(aged_bed$end),
                   x = 0.5, width = 3, assembly = "mm10")

png(file.path(projdir,sprintf("plotgardener_virtual4C_%s.png",gene)),res=300,width=3.5,height=2,units="in")
## Create a plotgardener page
pageCreate(width = 3.5, height = 2, default.units = "inches")

young_signal <- plotSignal(
  data = young_bed, params = params, binSize=5000, ymax = 1.1, 
  range = c(0, max_signal),
  y = 0.25, width = 3, height = 0.75, linecolor = "#0072B5FF",
  just = c("left", "top"), default.units = "inches"
)
annoYaxis(
  plot = young_signal, axisLine = T,
  at = pretty(seq(0, max_signal, length.out=4)),
  fontsize = 10
)
aged_signal <- plotSignal(
  data = aged_bed, params = params, binSize=5000, ymax = 1.1,
  range = c(0, max_signal), 
  y = 0.25, width = 3, height = 0.75, linecolor = "#BC3C29FF",
  just = c("left", "top"), default.units = "inches"
)

annoHighlight(
  plot = young_signal, 
  chrom=chrom, chromstart=gene_start, chromend=gene_end,
  fill = "#7ecdbb", alpha=0.4,
  y = 0.25, height = 0.75, just = c("left", "top"),
  default.units = "inches"
)

plotText(
  label = "Young", fonsize = 10, fontcolor = "#0072B5FF",
  params=params,
  x=3.5, y = 0.25, just = c("right", "top"),
  default.units = "inches"
)
plotText(
  label = "Aged", fonsize = 10, fontcolor = "#BC3C29FF",
  params=params,
  x=3.5, y = 0.45, just = c("right", "top"),
  default.units = "inches"
)

young_domainsPlot = plotRanges(young.domain.bed, linecolor="black",
                               params = params, collapse=T, 
                               fill = sapply(strsplit(young.domain.bed$V9, ","), 
                                             function(x) rgb(x[1], x[2], x[3], maxColorValue=255)),
                               y="0.025b", width=3, height=0.075, just = c("left", "top"),
                               default.units = "inches")
aged_domainsPlot = plotRanges(aged.domain.bed, linecolor="black",
                              params = params, collapse=T, 
                              fill = sapply(strsplit(aged.domain.bed$V9, ","), 
                                            function(x) rgb(x[1], x[2], x[3], maxColorValue=255)),
                              y="0b", width=3, height=0.075, just = c("left", "top"),
                              default.units = "inches")

genesPlot <- plotGenes(
  params = params,
  geneHighlights = data.frame(
    "gene" = c(gene),
    "color" = c("#225EA8")
  ),
  fontsize=8, stroke=0.01,
  geneBackground = "grey",
  y = "0.025b", height = 0.4,
  just = c("left", "top"), default.units = "inches"
)

plotGenomeLabel(
  params=params, scale = "Mb",
  assembly = "mm10",
  y = "0.025b", length = 3,
  default.units = "inches"
)

pageGuideHide()
dev.off()