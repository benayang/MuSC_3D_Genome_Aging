library(cicero) # needs to be version 1.3.5
library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(Signac)

projdir = '/nas/homes/benyang/HiC/13_MultiOme'

load("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_young_MuSC_cicero.RData")
load("/nas/homes/benyang/HiC/13_MultiOme/integrated_data/all_reps_aged_MuSC_cicero.RData")

# get the chromosome sizes from the Seurat object
genome <- seqlengths(aged_MuSC_data[["ATAC"]])

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
run_cicero = function(cicero_obj, genome.df, age) {
    conns <- cicero::run_cicero(cicero_obj, genomic_coords = genome.df, window = 5e5, sample_num = 100)
    ccans <- cicero::generate_ccans(conns)
    links <- ConnectionsToLinks(conns = conns, ccans = ccans, sep=c("-","-"))

    save(conns, ccans, links, file=paste0(age,"_MuSC_cicero_conns_ccans.RData"))
}

run_cicero(young_MuSC_cds.cicero, genome.df, "all_reps_Young")
run_cicero(aged_MuSC_cds.cicero, genome.df, "all_reps_Aged")