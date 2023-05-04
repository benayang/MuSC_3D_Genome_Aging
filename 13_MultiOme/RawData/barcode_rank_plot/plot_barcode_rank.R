#!/usr/bin/env Rscript
library(optparse)
library(BiocParallel)
library(EnsDb.Mmusculus.v79)

option_list = list(
    make_option("--indir", type="character", default=NULL, 
                help="input directory", metavar="character"),
    make_option("--ATAC", type="character", default=NULL, 
                help="input solo ATAC file with MACS2 peaks", metavar="character"),            
    make_option("--outdir", type="character", default=NULL, 
                help="output directory", metavar="character"),
    make_option("--prefix", type="character", default=NULL, 
                help="sample prefix for file outputs", metavar="character"),
    make_option("--maxlines", type="integer", default=NULL, 
                help="max.lines for Fragment files", metavar="integer"),
    make_option("--threads", type="integer", default=NULL, 
                help="number of threads", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);

param <- SnowParam(workers = args$threads, type = "SOCK")

suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))

#######################################
#### Get RNA Data
#######################################

raw = Read10X(file.path(args$indir,"raw_feature_bc_matrix"))
raw_RNA = CreateSeuratObject(raw[['Gene Expression']])
raw_RNA_sce = as.SingleCellExperiment(raw_RNA)

#######################################
#### Get ATAC Data
#######################################

# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# seqlevelsStyle(annotations) <- 'UCSC'
# genome(annotations) <- "mm10"

# grange.counts <- StringToGRanges(rownames(raw[['Peaks']]), sep = c(":", "-"))
# grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
# atac_counts <- raw[['Peaks']][as.vector(grange.use), ]

# frag.file <- file.path(args$indir,"atac_fragments.tsv.gz")
# chrom_assay <- CreateChromatinAssay(
#     counts = atac_counts,
#     sep = c(":", "-"),
#     genome = 'mm10',
#     fragments = frag.file,
#     annotation = annotations,
#     max.lines = args$max.lines
# )

# raw_ATAC = CreateSeuratObject(chrom_assay)
# raw_ATAC_sce = as.SingleCellExperiment(raw_ATAC)

raw_ATAC = readRDS(args$ATAC)
DefaultAssay(raw_ATAC) = "ATAC"
raw_ATAC = DietSeurat(raw_ATAC, assay="ATAC")
raw_ATAC_sce = as.SingleCellExperiment(raw_ATAC)

#######################################
#### barcode ranks
#######################################

rna.br.out <- barcodeRanks(raw_RNA_sce, BPPARAM = param)
atac.br.out <- barcodeRanks(raw_ATAC_sce, BPPARAM = param)

saveRDS(rna.br.out, file.path(args$outdir, sprintf("%s_RNA_barcode_rank.RDS", args$prefix)))
saveRDS(atac.br.out, file.path(args$outdir, sprintf("%s_ATAC_barcode_rank.RDS", args$prefix)))

png(file.path(args$outdir, sprintf("%s_RNA_barcode_rank.png", args$prefix)), res=300, units='in', width=5, height=5)
plot(rna.br.out$rank, rna.br.out$total, log="xy", xlab="Barcode Rank", ylab="Total UMI Per Barcode")
o <- order(rna.br.out$rank)
lines(rna.br.out$rank[o], rna.br.out$fitted[o], col="red")
abline(h=metadata(rna.br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(rna.br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
dev.off()

png(file.path(args$outdir, sprintf("%s_ATAC_barcode_rank.png", args$prefix)), res=300, units='in', width=5, height=5)
plot(atac.br.out$rank, atac.br.out$total, log="xy", xlab="Barcode Rank", ylab="Total ATAC Transposition Events Per Barcode")
o <- order(atac.br.out$rank)
lines(atac.br.out$rank[o], atac.br.out$fitted[o], col="red")
abline(h=metadata(atac.br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(atac.br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
    legend=c("knee", "inflection"))
dev.off()

#######################################
#### Find empty droplets
#######################################

set.seed(100)
rna.e.out <- emptyDrops(raw_RNA_sce, BPPARAM = param)
rna.is.cell <- rna.e.out$FDR <= 0.01
print(table(Limited=rna.e.out$Limited, Significant=rna.is.cell))
saveRDS(rna.e.out, file.path(args$outdir, sprintf("%s_RNA_empty_droplets.RDS", args$prefix)))

atac.e.out <- emptyDrops(raw_ATAC_sce, BPPARAM = param)
atac.is.cell <- atac.e.out$FDR <= 0.01
print(table(Limited=atac.e.out$Limited, Significant=atac.is.cell))
saveRDS(atac.e.out, file.path(args$outdir, sprintf("%s_ATAC_empty_droplets.RDS", args$prefix)))

png(file.path(args$outdir, sprintf("%s_RNA_empty_droplets.png", args$prefix)), res=300, units='in', width=5, height=5)
plot(rna.e.out$Total, -rna.e.out$LogProb, col=ifelse(rna.is.cell, "red", "black"),
    xlab="Total UMI count", ylab="-Log Probability")
dev.off()

png(file.path(args$outdir, sprintf("%s_ATAC_empty_droplets.png", args$prefix)), res=300, units='in', width=5, height=5)
plot(atac.e.out$Total, -atac.e.out$LogProb, col=ifelse(atac.is.cell, "red", "black"),
    xlab="Total ATAC Transposition Event Count", ylab="-Log Probability")
dev.off()