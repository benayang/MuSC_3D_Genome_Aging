library(ArchR)
library(dplyr)
library(tidyr)
library(parallel)
library(RColorBrewer)
library(scales)

#ArchR::installExtraPackages()
projdir = '/nas/homes/benyang/HiC/13_MultiOme/ArchR_analysis'

addArchRGenome("mm10")
addArchRThreads(threads = 45) 

load(file.path(projdir,"MuSC_ArchR","all_MuSC","all_MuSC_MuSC_cicero_conns_ccans.RData"))
MuSC_cicero = readRDS(file.path(projdir,"MuSC_ArchR","all_MuSC","all_MuSC_cicero.RDS"))
MuSC_projAging1 = readRDS(file.path(projdir, "Save-MuSC_projAging1-01", "Save-ArchR-Project.rds"))

# use a manual threshold
ccans2 = cicero::generate_ccans(conns, 0.1)
ccans2 = ccans2 %>% separate(Peak, into=c("seqnames","start","end"), sep="_", remove=F)
ccans2_gr = GRanges(ccans2[,c("seqnames","start","end")])
mcols(ccans2_gr)$CCAN = ccans2$CCAN

MuSC_peakset = getPeakSet(MuSC_projAging1)
MuSC_peakmatrix = getMatrixFromProject(MuSC_projAging1, "PeakMatrix")

MuSC_peakset_string = data.frame(MuSC_peakset) %>% mutate(peak = paste(seqnames, start, end, sep="_")) %>% pull(peak)

MuSC_peakmatrix_data = assay(MuSC_peakmatrix)
rownames(MuSC_peakmatrix_data) = MuSC_peakset_string

young_MuSC_cellnames = rownames(colData(MuSC_peakmatrix))[colData(MuSC_peakmatrix)$Age != "Aged"]
aged_MuSC_cellnames = rownames(colData(MuSC_peakmatrix))[colData(MuSC_peakmatrix)$Age == "Aged"]

# MuSC_peakmatrix_comparison = data.frame(peak = MuSC_peakset_string,
#                                         young = Matrix::rowMeans(MuSC_peakmatrix_data[,young_MuSC_cellnames]),
#                                         aged = Matrix::rowMeans(MuSC_peakmatrix_data[,aged_MuSC_cellnames]))

# conns_filt = conns[conns$coaccess>0.1, ]
# comparison_filt = MuSC_peakmatrix_comparison[MuSC_peakmatrix_comparison$peak %in% conns_filt$Peak1, ]

# all(conns_filt$Peak1 %in% MuSC_peakmatrix_comparison$peak)
# all(conns_filt$Peak2 %in% MuSC_peakmatrix_comparison$peak)

# conns_peaks = data.frame(peakstring = unique(conns_filt$Peak1))
# conns_peaks = conns_peaks %>% drop_na() %>% tidyr::separate(peakstring, into=c("chrom","start","end"), sep="_", remove=F) %>% GRanges()

MuSC_projAging1 = addFeatureMatrix(
  input = MuSC_projAging1,
  features = ccans2_gr,
  matrixName = "CCANFeatureMatrix",
  ceiling = 10^9,
  binarize = TRUE,
  verbose = TRUE,
  threads = getArchRThreads(),
  parallelParam = NULL,
  force = TRUE,
  logFile = createLogFile("addFeatureMatrix")
)

ccan_matrix = getMatrixFromProject(MuSC_projAging1, "CCANFeatureMatrix", binarize=T)
ccan_matrix_data = assay(ccan_matrix)
rownames(ccan_matrix_data) = rownames(ccans2_gr)

ccan_matrix_data_summary = summary(ccan_matrix_data)
ccan_matrix_data_summary$CCAN = mcols(ccans2_gr)$CCAN[ccan_matrix_data_summary$i]
ccan_matrix_data_summary$Age = colData(ccan_matrix)$Age[ccan_matrix_data_summary$j]

ccan_matrix_pheatmap = ccan_matrix_data_summary %>% group_by(Age,CCAN) %>% summarise(num_peaks=n()) %>% pivot_wider(names_from="Age", values_from="num_peaks")
pheatmap::pheatmap(log10(ccan_matrix_pheatmap[,2:3]), cluster_cols=F, show_rownames=F, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(200), filename=file.path(projdir,"MuSC_ArchR","all_cells","Plots","num_peaks_per_CCAN.png"))

MuSC_peakmatrix_comparison = data.frame(young = Matrix::rowMeans(ccan_matrix_data[,young_MuSC_cellnames]),
                                        aged = Matrix::rowMeans(ccan_matrix_data[,aged_MuSC_cellnames]))
rownames(MuSC_peakmatrix_comparison) = mcols(conns_peaks)$peakstring