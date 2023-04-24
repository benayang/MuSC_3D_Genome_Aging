library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggsci)
library(pheatmap)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\HOMER\\call_homer"


# Get expression data -----------------------------------------------------

tpm.data = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC", 
                                "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

avg.tpm.data = data.frame(A=rowMeans(tpm.data[,1:3]),
                          Y=rowMeans(tpm.data[,4:8]))

DE.genes = unlist(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC", 
                                       "compartmentExpression", "DE_0.05padj_log2(1.5)LFC_genenames.csv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)


# Get motifs --------------------------------------------------------------

motif_family = read.table(file.path(projdir,"tf_to_family_map.txt"), sep="\t", as.is = T, fill = T, header = T, col.names=c("TF","Family","Domain","Domain2"))
motif_family = motif_family %>% distinct()

# For loop anchors in TADs
aged_motifs = c("Zeb2","Zbtb7a","Neurog2","Pou2f2")
young_motifs = c("Ctcf","Tcf7l2","Mxi1")

motif_df = tpm.data[c(young_motifs,aged_motifs), c(1:3,5:7)] %>% drop_na()

pheatmap(as.matrix(motif_df), cluster_cols=F, scale="row", color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)))

avg_motif_df = data.frame(A=rowMeans(motif_df[,1:3]),
                          Y=rowMeans(motif_df[,4:6]))         
pheatmap(as.matrix(avg_motif_df), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)))

motif_fam_df = motif_family[motif_family$TF %in% c(aged_motifs,young_motifs),]

motif_family[motif_family$TF=="Pou2f1",]

curated_fam = c(unique(motif_fam_df$Family),"Pou")
curated_fam_df = motif_family[motif_family$Family %in% curated_fam,]

curated_fam_tpm = tpm.data[curated_fam_df$TF,] %>% drop_na()
curated_fam_tpm[rownames(curated_fam_tpm) %in% DE.genes.name,]
pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% DE.genes.name,c(1:3,5:7)]), 
         cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), scale="row", 
         cellwidth = 15, cellheight=10)

pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% DE.genes.name,c(5:7,1:3)]), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), scale="row", 
         filename = file.path(projdir,"TF_fam_DE_expression.png"), cellwidth = 15, cellheight=17, treeheight_row = 25)
pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% rownames(tpm.data),c(5:7,1:3)]), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), scale="row", 
         filename = file.path(projdir,"TF_fam_expression.png"), cellwidth = 15, cellheight=10)

pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% rownames(tpm.data),c(5,6,8,1:3)]), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), scale="row", 
         cellwidth = 15, cellheight=10)

# Get motifs with average TPM ---------------------------------------------

curated_fam_tpm = avg.tpm.data[curated_fam_df$TF,] %>% drop_na() %>% mutate(A=log2(A+1), Y=log2(Y+1))
curated_fam_tpm[rownames(curated_fam_tpm) %in% DE.genes.name,]
pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% rownames(tpm.data),]), 
         cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), 
         cellwidth = 15, cellheight=10)

pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% DE.genes.name,]), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), 
         filename = file.path(projdir,"TF_fam_DE_expression_average.png"), cellwidth = 15, cellheight=15)
pheatmap(as.matrix(curated_fam_tpm[rownames(curated_fam_tpm) %in% rownames(tpm.data),]), cluster_cols=F, color=rev(colorRampPalette(brewer.pal(9,"RdBu"))(200)), 
         filename = file.path(projdir,"TF_fam_expression_average.png"), cellwidth = 15, cellheight=10)
