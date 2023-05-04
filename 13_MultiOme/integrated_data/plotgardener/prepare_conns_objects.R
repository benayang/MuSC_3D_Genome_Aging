library(dplyr)
library(tidyr)
library(Signac)
library(Seurat)
library(cicero)

projdir = "/nas/homes/benyang/HiC/13_MultiOme/integrated_data"

load(file.path(projdir,"all_reps_Aged_MuSC_cicero_conns_ccans.RData"))
aged_conns = conns
aged_ccans = ccans
load(file.path(projdir,"all_reps_Young_MuSC_cicero_conns_ccans.RData"))
young_conns = conns
young_ccans = ccans

rm(list=c("conns","ccans"))
gc()

get_overlapping_coords = function(conns, region) {
    target_conns = find_overlapping_coordinates(conns$Peak1, region)
    target_conns_idx = which(conns$Peak1 %in% target_conns)
    target_conns = conns[target_conns_idx, ]
    target_conns = conns %>% 
                    separate(Peak1, into=c("chrom1","start1","end1"), remove=T) %>% 
                    separate(Peak2, into=c("chrom2","start2","end2"), remove=T)
    stopifnot(identical(target_conns$chrom1, target_conns$chrom2))
    target_conns = target_conns %>% mutate(start1=as.integer(start1), end1=as.integer(end1), start2=as.integer(start2), end2=as.integer(end2))
    target_conns = target_conns[target_conns$coaccess>0, ]

    return(target_conns)
}

aged_myod_conns = get_overlapping_coords(aged_conns, "chr7-46350000-46479092")
young_myod_conns = get_overlapping_coords(young_conns, "chr7-46350000-46479092")

write.table(aged_myod_conns, file.path(projdir,"plotgardener","aged_myod_conns.txt"), sep='\t', quote=F, row.names=F, col.names=F)
write.table(young_myod_conns, file.path(projdir,"plotgardener","young_myod_conns.txt"), sep='\t', quote=F, row.names=F, col.names=F)