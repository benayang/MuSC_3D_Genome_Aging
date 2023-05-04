library(monocle3)
library(cicero)
library(tidyr)
library(dplyr)

projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR'

projAging6 = readRDS("/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/Save-projAging6-01/Save-ArchR-Project.rds")

young_atac_cds.cicero = readRDS(file.path(projdir, "young_MuSC", "young_MuSC_cicero.RDS"))
aged_atac_cds.cicero = readRDS(file.path(projdir, "aged_MuSC", "aged_MuSC_cicero.RDS"))
MuSC_atac_cds.cicero = readRDS(file.path(projdir, "all_MuSC", "all_MuSC_cicero.RDS"))

genome <- as.data.frame(ArchR::getChromSizes(projAging6))

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = genome$seqnames, "length" = genome$width)

# run cicero
run_cicero = function(cicero_obj, genome.df, age, outdir) {
    conns <- cicero::run_cicero(cicero_obj, genomic_coords = genome.df, window = 5e5, sample_num = 100)
    ccans <- cicero::generate_ccans(conns)

    save(conns, ccans, file=file.path(outdir, paste0(age,"_MuSC_cicero_conns_ccans.RData")))
}

run_cicero(young_atac_cds.cicero, genome.df, "Young", file.path(projdir, "young_MuSC"))
run_cicero(aged_atac_cds.cicero, genome.df, "Aged", file.path(projdir, "aged_MuSC"))
run_cicero(MuSC_atac_cds.cicero, genome.df, "all_MuSC", file.path(projdir, "all_MuSC"))