projdir='/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis/MuSC_ArchR'

MuSC_projAging1 = readRDS(file.path(projdir,"all_MuSC","Save-ArchR-Project.rds"))
cellcoldata = getCellColData(MuSC_projAging1)

aged_cell_df = data.frame(cells=rownames(cellcoldata)[cellcoldata$Sample=="Aged"], group="Aged_MuSC")
aged_cell_df[,'cells'] = gsub("Aged#","",aged_cell_df[,'cells'],fixed=T)
young_cell_df = data.frame(cells=rownames(cellcoldata)[cellcoldata$Sample=="Young"], group="Young_MuSC")
young_cell_df[,'cells'] = gsub("Young#","",young_cell_df[,'cells'],fixed=T)
young_v2_cell_df = data.frame(cells=rownames(cellcoldata)[cellcoldata$Sample=="Young_v2"], group="Young_v2_MuSC")
young_v2_cell_df[,'cells'] = gsub("Young_v2#","",young_v2_cell_df[,'cells'],fixed=T)

write.table(aged_cell_df, file.path(projdir, "bigwigs","aged_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(young_cell_df, file.path(projdir, "bigwigs","young_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(young_v2_cell_df, file.path(projdir, "bigwigs","young_v2_MuSC_cell_barcodes.tsv"), sep="\t", row.names=F, col.names=F, quote=F)