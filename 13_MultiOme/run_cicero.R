library(cicero) # needs to be version 1.3.5
library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(Signac)

projdir = '/nas/homes/benyang/HiC/13_MultiOme'

load(file.path(projdir,"MuSC_aged_subset_singlets.cicero.RData"))

# get the chromosome sizes from the Seurat object
genome <- seqlengths(MuSC_aged_subset_singlets)

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(MuSC_cds.cicero, genomic_coords = genome.df, window = 5e5, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans, sep=c("_","_"))
Links(MuSC_aged_subset_singlets) <- links

save(conns, ccans, MuSC_aged_subset_singlets, file=file.path(projdir,"MuSC_aged_subset_singlets.cicero.connections.RData"))