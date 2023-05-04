#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("Need to pass at least one argument.", call.=F)
} else {
    df = read.table(args[1], sep='\t', header=F)
    df$V7 = unlist(sapply(df$V7, function(x) ifelse(is.na(x), NA, unlist(strsplit(x, split=".", fixed=T))[1])))
    write.table(df, args[2], sep='\t', row.names=F, col.names=F, quote=F)
}