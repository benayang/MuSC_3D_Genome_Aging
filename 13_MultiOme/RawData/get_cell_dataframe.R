projdir='/nas/homes/benyang/HiC/13_MultiOme'

aged_MuSC = readRDS(file.path(projdir, "all_reps_aged_MuSC_subset_linkPeaks.RDS"))
aged_cell_df = data.frame(cells=Cells(aged_MuSC), group="all_reps_aged_MuSC")
write.table(aged_cell_df, file.path(projdir, "RawData","all_reps_aged_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)

young_MuSC = readRDS(file.path(projdir, "all_reps_young_MuSC_subset_linkPeaks.RDS"))
young_cell_df = data.frame(cells=colnames(young_MuSC)[grep("_1",colnames(young_MuSC),fixed=T)], group="all_reps_young_MuSC")
young_cell_df$cells = sapply(young_cell_df$cells, function(x) gsub("_1","",x,fixed=T))
young_v2_cell_df = data.frame(cells=colnames(young_MuSC)[grep("_2",colnames(young_MuSC),fixed=T)], group="all_reps_young_v2_MuSC")
young_v2_cell_df$cells = sapply(young_v2_cell_df$cells, function(x) gsub("_2","",x,fixed=T))
write.table(young_cell_df, file.path(projdir, "RawData","all_reps_young_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(young_v2_cell_df, file.path(projdir, "RawData","all_reps_young_v2_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)