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

MuSC_projAging1 = readRDS(file.path(projdir, "Save-MuSC_projAging1-01", "Save-ArchR-Project.rds"))

exp_mat = getMatrixFromProject(MuSC_projAging1, 'GeneIntegrationMatrix')
score_mat = getMatrixFromProject(MuSC_projAging1, 'GeneScoreMatrix')

common_genes = intersect(rowData(exp_mat)$name, rowData(score_mat)$name)
exp_common_idx = na.omit(match(rowData(exp_mat)$name, common_genes))
score_common_idx = na.omit(match(rowData(score_mat)$name, common_genes))

exp_mat_avg = Matrix::rowMeans(assay(exp_mat)[exp_common_idx, ])
score_mat_avg = Matrix::rowMeans(assay(score_mat)[score_common_idx, ])

plt_df = data.frame(avg_exp = exp_mat_avg, avg_score = score_mat_avg)
plt_df$score_percentile_bin = cut(rank(-plt_df$avg_score)/nrow(plt_df), breaks=seq(0,1,by=0.25))

ggplot(plt_df, aes(x=score_percentile_bin, y=log2(exp_mat_avg+1))) +
stat_summary(fun.data="mean_se", geom="pointrange") +
theme_bw()
ggsave(file.path(projdir, "MuSC_ArchR", "all_MuSC", "Plots", "expression_vs_binned_activity.png"), dpi=300, width=4, height=4)