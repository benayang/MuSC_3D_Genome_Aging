library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\17_TOBIAS"

results = read.table(file.path(projdir, "merged_peakset_annotated_expressed_TF", "bindetect_merged_peakset_expressed_results.txt"), sep='\t', header=T)
results = results %>% mutate(change_percentile = rank(young_merged_ATAC_aged_merged_ATAC_change)/n(),
                             pval_percentile = rank(-log10(young_merged_ATAC_aged_merged_ATAC_pvalue))/n())

sig_id = results %>% dplyr::filter(change_percentile > 0.95 | change_percentile < 0.05 | pval_percentile > 0.95) %>% pull(output_prefix)
# max_id = results %>% dplyr::slice_max(order_by=young_merged_ATAC_aged_merged_ATAC_change, prop=0.05) %>% pull(output_prefix)
# min_id = results %>% dplyr::slice_min(order_by=young_merged_ATAC_aged_merged_ATAC_change, prop=0.05) %>% pull(output_prefix)

results$sig = "NS"
results$sig[results$output_prefix %in% sig_id & results$young_merged_ATAC_aged_merged_ATAC_change > 0] = "Young"
results$sig[results$output_prefix %in% sig_id & results$young_merged_ATAC_aged_merged_ATAC_change < 0] = "Aged"

summary(results$total_tfbs[results$sig!="NS"])

results %>% dplyr::filter(sig=="Young") %>% pull(name) %>% writeClipboard()

ggplot(results, aes(x = young_merged_ATAC_aged_merged_ATAC_change, 
                    y = -log10(young_merged_ATAC_aged_merged_ATAC_pvalue))) +
  geom_point(aes(color=sig), alpha=0.8) +
  geom_label_repel(data=results[results$sig != "NS", ], aes(label=name, fill=sig), 
                  segment.size = 0.5, color="white", segment.color="black",
                  max.overlaps=35, force = 5, size=3) +
  scale_color_manual(values=c(pal_nejm()(2)[1],"black",pal_nejm()(2)[2])) +
  scale_fill_manual(values=pal_nejm()(2)) +
  scale_x_continuous(breaks=seq(-0.4,0.4,0.2)) +
  #scale_y_continuous(breaks=seq(0,1500,250)) +
  theme_pubr() + guides(fill="none", color=guide_legend(title="")) +
  labs(x="Differential TF Binding Score (Young/Aged)", y=expression(bold(paste("-lo","g"["10"],"(p-value)")))) +
  theme(axis.title = element_text(size=12, face="bold"),
        axis.text = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = c(0.1,0.15))
ggsave(file.path(projdir, "merged_peakset_annotated_expressed_volcano.png"), dpi=300, width=6, height=6)

filtered_TF = c(grep("Zic",results$name,ignore.case=T),
                grep("CTCF",results$name,ignore.case=T),
                grep("Nr1d",results$name,ignore.case=T),
                grep("Ror",results$name,ignore.case=T),
                grep("Dbp",results$name,ignore.case=T),
                grep("KLF",results$name,ignore.case=T),
                grep("HES",results$name,ignore.case=T),
                grep("^SP",results$name))

ggplot(results, aes(x = young_merged_ATAC_aged_merged_ATAC_change, 
                    y = -log10(young_merged_ATAC_aged_merged_ATAC_pvalue))) +
  geom_point(aes(color=sig), alpha=0.8) +
  geom_label_repel(data=results[filtered_TF, ] %>% dplyr::filter(sig!="NS"), 
                   aes(label=name, fill=sig), 
                   segment.size = 0.5, color="white", segment.color="black",
                   max.overlaps=35, force = 5, size=3) +
  scale_color_manual(values=c(pal_nejm()(2)[1],"black",pal_nejm()(2)[2])) +
  scale_fill_manual(values=pal_nejm()(2)) +
  scale_x_continuous(breaks=seq(-0.4,0.4,0.2)) +
  #scale_y_continuous(breaks=seq(0,1500,250)) +
  theme_pubr() + guides(fill="none", color=guide_legend(title="")) +
  labs(x="Differential TF Binding Score (Young/Aged)", y=expression(bold(paste("-lo","g"["10"],"(p-value)")))) +
  theme(axis.title = element_text(size=12, face="bold"),
        axis.text = element_text(size=12, face="bold"),
        legend.text = element_text(size=10, face="bold"),
        legend.position = c(0.15,0.175))
ggsave(file.path(projdir, "merged_peakset_annotated_expressed_filtered_volcano.png"), dpi=300, width=4, height=4)

# Gene expression ---------------------------------------------------------

library(Hmisc)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC"

tpm.data = readRDS(file.path(projdir, "filtered_tpm_data.RDS"))

# tpm.data = read.table(file.path(projdir, "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
# gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
# gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
# duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]
# 
# tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
# rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
# tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
# tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

DE.genes = unlist(read.table(file.path(projdir, "compartmentExpression", "DE_0.05padj_log2(1.5)LFC_genenames.tsv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)
DE.data = read.delim2("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\diff_d0_Y_vs_A.tsv")

sig_TFs = results %>% dplyr::filter(sig!="NS") %>% pull(name) %>% unique() %>% tolower() %>% capitalize()
sig_TFs[sig_TFs %in% DE.genes.name]
# sig_TFs = c(sig_TFs, "Nr1h3", "Rxra", "Fosl2", "Jund", "Zic1", "Zic2")
# sig_TFs = unique(sig_TFs)

pheatmap(tpm.data['Myod1',c(1:3,5:7)], cluster_rows = F, cluster_cols = F, scale="row")

exp_df <- tpm.data[sig_TFs, ] %>% 
  dplyr::select(1:3,5:7) %>% 
  drop_na()
exp_df <- exp_df[!rownames(exp_df)=='Hes1',]
exp_df$Aged <- rowMeans(exp_df[,1:3])
exp_df$Young <- rowMeans(exp_df[,5:7])
pheatmap(log2(exp_df[,c("Aged","Young")]), cluster_cols=F, color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200), 
         cellwidth = 30, cellheight = 10, scale="row",
         filename = file.path(projdir, "merged_peakset_annotated_expressed_sigTFs_pooled_hmp.png"))

tpm.data[sig_TFs, ] %>% drop_na() %>% dplyr::select(1:3,5:7) %>% 
  pheatmap(cluster_cols=F, scale="row", color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200), 
           cellwidth = 30, cellheight = 10, 
           filename = file.path(projdir, "merged_peakset_annotated_expressed_sigTFs_hmp.png"))

pheatmap(tpm.data[capitalize(tolower(results[filtered_TF, ]$name)), c(1:3,5:7)] %>% drop_na(),
         cluster_cols=F, scale="row", color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200))
pheatmap(tpm.data[c("Tcfl5","Tcf7l1","Mybl1","Srebf1","Srebf2","Thap1","Pou2f2","Foxc1","Klf16","Foxo3","Max","Maf","Bcl6"), c(1:3,5:7)] %>% drop_na(),
         cluster_cols=F, scale="row", color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200))

# sig_TF_df = data.frame(young = tpm.data[sig_TFs,] %>% drop_na() %>% dplyr::select(1:3) %>% rowMeans(),
#                        aged = tpm.data[sig_TFs,] %>% drop_na() %>% dplyr::select(5:7) %>% rowMeans()) %>%
#   mutate(log2fc = log2(young/aged))
# 
# sig_TF_plt_df = data.frame(x=1, y=rownames(sig_TF_df), value=sig_TF_df$log2fc)

sig_TF_plt_df = results[results$name %in% toupper(rownames(sig_TF_df)), ] %>% dplyr::filter(sig!="NS") %>%
  dplyr::select(name, output_prefix, young_merged_ATAC_aged_merged_ATAC_change) %>%
  dplyr::filter(!(output_prefix %in% c("CTCF_MA1929.1","CTCF_MA1930.1"))) %>%
  distinct() %>%
  dplyr::select(name, young_merged_ATAC_aged_merged_ATAC_change) %>%
  rename(young_merged_ATAC_aged_merged_ATAC_change = "TF") %>%
  tibble::column_to_rownames("name") 
rownames(sig_TF_plt_df) = capitalize(tolower(rownames(sig_TF_plt_df)))
  
tf_plt = pheatmap(sig_TF_plt_df, cluster_cols = F, color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200))
rownames(sig_TF_plt_df)[tf_plt$tree_row[["order"]]]

pheatmap(data.frame(FC=sig_TF_df[rownames(sig_TF_plt_df)[tf_plt$tree_row[["order"]]],"log2fc"], row.names=rownames(sig_TF_plt_df)[tf_plt$tree_row[["order"]]]), 
         cluster_cols = F, cluster_rows=F,
         color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200))

tpm.data[sig_TFs, c(1:3,5:7)] %>% drop_na() %>% 
  pheatmap(cluster_cols = F, scale="row",
           color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200),
           width=4, height=14, border_color = NULL, fontsize = 12,
           filename = file.path(projdir,"TF_expression.png"))


# Get TF families ---------------------------------------------------------

motif_family = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\HOMER\\call_homer",
                                    "tf_to_family_map.txt"), sep="\t", as.is = T, fill = T, header = T, col.names=c("TF","Family","Domain","Domain2"))
motif_family = motif_family %>% distinct()

family_TFs = motif_family %>%
  dplyr::filter(Family %in% unique(motif_family$Family[motif_family$TF %in% sig_TFs])) %>%
  pull(TF) %>% unique()

tpm.data[family_TFs, c(1:3,5:7)] %>% drop_na() %>% 
  pheatmap(scale="row", color = colorRampPalette(rev(brewer.pal("RdBu",n=11)))(200),
           width=4, height=14, border_color = NULL, fontsize_row = 8,
           filename = file.path(projdir,"TF_family_expression.png"))


# TF clustering -----------------------------------------------------------

library(dendextend)
library(pheatmap)

clustering = read.delim2(file.path(projdir,"merged_peakset_annotated_expressed_TF","bindetect_merged_peakset_expressed_distances.txt"))

idx = which(colnames(clustering) %in% results$output_prefix[results$sig!="NS"])

clustering_mat = matrix(as.numeric(unlist(clustering[idx, idx])), nrow=length(idx), ncol=length(idx))
rownames(clustering_mat) = colnames(clustering)[idx]
colnames(clustering_mat) = colnames(clustering)[idx]

clustering_mat_plt = pheatmap(clustering_mat, clustering_method = "ward.D2")

plot(clustering_mat_plt$tree_col)

dend = clustering_mat_plt$tree_col %>% as.dendrogram()

dend_colors = data.frame(TF = clustering_mat_plt$tree_row$labels[clustering_mat_plt$tree_row$order])
dend_colors$Direction = results[results$sig!="NS", ]$sig[match(dend_colors$TF,results[results$sig!="NS", ]$output_prefix)]
dend_colors$color = sapply(dend_colors$Direction, function(x) ifelse(x=="Young",pal_nejm()(2)[2],pal_nejm()(2)[1]))

png(file.path(projdir,"merged_peakset_annotated_expressed_TF","sig_TF_dendrogram.png"), res=300, units='in', width=5, height=9)
par(cex=0.7, cex.axis=1.5, cex.lab=1.5, mar=c(5,3,0,10))
dend %>% set("labels_col", value = dend_colors$color) %>% plot(horiz=T, xlab="Complete Distance")
dev.off()
