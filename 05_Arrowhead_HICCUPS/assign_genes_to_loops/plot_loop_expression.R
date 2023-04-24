library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(rstatix)
library(VennDiagram)
library(tidyr)
library(scales)
library(dplyr)
library(ggpattern)

# Import data -------------------------------------------------------------

expression_projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer"
projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS"

tpm.data = read.table(file.path(expression_projdir, "TAD expression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with non-negative TPM  

DE.genes = unlist(read.table(file.path(expression_projdir, "TAD expression", "DE_0.05padj_log2(1.5)LFC_genenames.csv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)

# Plot TPM within loop domains  ------------------------------------------------

genenames.dir = file.path(projdir, "assign_genes_to_loops")

gene_df = data.frame()
for(suffix in c(".loop_anchor.knownGenes.bed",".loop_domain.knownGenes.bed",
                ".diff_loop_anchor.knownGenes.bed",".diff_loop_domain.knownGenes.bed")) {
  tasks=grep(pattern=suffix,list.files(genenames.dir),fixed=T,value=T)
  for(taskname in tasks) {
    genes = unique(read.table(file.path(genenames.dir, taskname))$V4)
    Age = ifelse(grepl(pattern="aged",taskname,fixed=T), "Aged", "Young")
    LoopType = ifelse(grepl(pattern="diff",taskname,fixed=T), "Diff", "Reg")
    LoopRegion = ifelse(grepl(pattern="domain",taskname,fixed=T), "Domain", "Anchor")
    tmp = data.frame(filename=taskname, genes=genes, Age=Age, LoopType=LoopType, LoopRegion=LoopRegion)
    gene_df = rbind(gene_df, tmp)
  }
}
gene_df = gene_df %>% filter(genes %in% rownames(tpm.data))
avg.tpm = data.frame(Aged=rowMeans(tpm.data[,1:3]),
                     Young=rowMeans(tpm.data[,5:7]))
#gene_df = cbind(gene_df,tpm.data[gene_df$genes,1:8])
gene_df$Young = avg.tpm[gene_df$genes,"Young"]
gene_df$Aged = avg.tpm[gene_df$genes,"Aged"]

gene_df.plt.df = gene_df %>% 
  pivot_longer(cols=c(Young,Aged), names_to="TPM.Age", values_to="TPM") %>%
  dplyr::select(genes,Age,LoopType,LoopRegion,TPM.Age,TPM)

uncovered.aged.genes = avg.tpm[!(rownames(avg.tpm) %in% pull(gene_df.plt.df[gene_df.plt.df$Age=="Aged" & gene_df.plt.df$LoopType=="Reg","genes"])), ] %>% 
  tibble::rownames_to_column("genes") %>% mutate(Age="Aged")
uncovered.young.genes = avg.tpm[!(rownames(avg.tpm) %in% pull(gene_df.plt.df[gene_df.plt.df$Age=="Young" & gene_df.plt.df$LoopType=="Reg","genes"])), ] %>% 
  tibble::rownames_to_column("genes") %>% mutate(Age="Young")
uncovered.genes = rbind(uncovered.aged.genes, uncovered.young.genes) %>% distinct()
uncovered.genes = uncovered.genes %>% pivot_longer(cols=c(Young,Aged),names_to="TPM.Age", values_to="TPM")
uncovered.genes = uncovered.genes %>% mutate(LoopType="Uncovered",LoopRegion="Uncovered") %>% dplyr::select(genes,Age,LoopType,LoopRegion,TPM.Age,TPM)

#all.genes.df = rbind(gene_df.plt.df, uncovered.genes) %>% filter(Age==TPM.Age) %>% mutate(Age=factor(Age,levels=c("Young","Aged"))) %>% distinct()
all.genes.df = rbind(gene_df.plt.df, uncovered.genes) %>% mutate(Age=factor(Age,levels=c("Young","Aged")), Group=paste(Age,LoopType,LoopRegion,sep="_")) %>% distinct()

all.genes.df %>% distinct() %>% ggplot(aes(x=TPM.Age,y=log2(TPM))) + facet_wrap(~LoopRegion+LoopType) + geom_boxplot(aes(fill=TPM.Age))

all.genes.df %>% filter(LoopType=="Diff" & Age=="Young") %>% pivot_wider(names_from=TPM.Age,values_from=TPM,id_cols = genes) %>% 
  tibble::column_to_rownames("genes") %>% mutate(Young=log2(Young),Aged=log2(Aged)) %>% 
  pheatmap(cluster_cols = F,angle_col = 45, file=file.path(projdir,"Young_diff_loops_hmp.png"),dpi=300,width=4,height=13)
all.genes.df %>% filter(LoopType=="Diff" & Age=="Aged") %>% pivot_wider(names_from=TPM.Age,values_from=TPM,id_cols = genes) %>% 
  tibble::column_to_rownames("genes") %>% mutate(Young=log2(Young),Aged=log2(Aged)) %>% 
  pheatmap(cluster_cols = F,angle_col = 45, file=file.path(projdir,"Aged_diff_loops_hmp.png"),dpi=300,width=4,height=8)

gene.age.test = all.genes.df %>%
  filter(LoopType != "Diff" & LoopRegion!="Anchor") %>%
  mutate(log2tpm=log2(TPM)) %>%
  group_by(LoopType) %>%
  wilcox_test(log2tpm~Age) %>%
  add_xy_position(x="Age",group="LoopType",dodge=0.9)

gene.loop.test = all.genes.df %>%
  filter(LoopType != "Diff" & LoopRegion!="Anchor") %>%
  mutate(log2tpm=log2(TPM)) %>%
  group_by(Age) %>%
  wilcox_test(log2tpm~LoopType) %>%
  add_xy_position(x="Age",group="LoopType",dodge=0.9)
gene.loop.test = gene.loop.test %>% mutate(new.label=rep("<2e-16",2))

all.genes.df %>%
  filter(LoopType != "Diff" & LoopRegion!="Anchor") %>%
  group_by(Age,LoopType,TPM.Age) %>%
  summarise(count=n())

all.genes.df %>%
  #filter(LoopType != "Diff" & LoopRegion!="Anchor") %>%
  filter(LoopRegion=="Anchor") %>%
  ggplot(aes(x=Age,y=log2(TPM))) +
  geom_violin_pattern(aes(fill=Age,pattern=LoopType), color="black",
                      pattern_density = 0.1,
                      pattern_fill    = 'white',
                      pattern_colour  = 'black',
                      position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,LoopType)), color="black", width=0.2, outlier.shape=NA, position=position_dodge(0.9)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  #stat_pvalue_manual(gene.age.test, tip.length = 0.01, bracket.nudge.y = 3.5, step.increase = 0.05) +
  #stat_pvalue_manual(gene.loop.test, label = "new.label", tip.length = 0.01, bracket.nudge.y = 0.5) +
  #scale_fill_manual(values=rev(pal_nejm()(2))) +
  #scale_pattern_manual(values=c("none","stripe"), labels=c("Loop","Non-loop")) +
  theme_bw() + guides(fill="none", pattern=guide_legend(title=NULL)) +
  labs(x=NULL,y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(legend.position = "top",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"loop_expression.png"),dpi=300,width=3,height=4)


# Write expression of genes -----------------------------------------------

all.genes.df %>% filter(LoopRegion=="Domain" & LoopType=="Diff" & Age=="Aged") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","aged_diff_loops_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Domain" & LoopType=="Diff" & Age=="Young") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","young_diff_loops_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Anchor" & LoopType=="Diff" & Age=="Aged") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","aged_diff_loops_anchor_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Anchor" & LoopType=="Diff" & Age=="Young") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","young_diff_loops_anchor_TPM.tsv"),sep="\t",quote=F,row.names=F)

all.genes.df %>% filter(LoopRegion=="Domain" & LoopType=="Reg" & Age=="Aged") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","aged_loops_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Domain" & LoopType=="Reg" & Age=="Young") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","young_loops_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Anchor" & LoopType=="Reg" & Age=="Aged") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","aged_loops_anchor_TPM.tsv"),sep="\t",quote=F,row.names=F)
all.genes.df %>% filter(LoopRegion=="Anchor" & LoopType=="Reg" & Age=="Young") %>% pivot_wider(names_from=TPM.Age, values_from=TPM) %>% select(genes,Young,Aged) %>%
  mutate(zscore=scale(Aged/Young)) %>% write.table(file.path(projdir,"assign_genes_to_loops","young_loops_anchor_TPM.tsv"),sep="\t",quote=F,row.names=F)

# Look at expression of quiescence genes ----------------------------------

library(pheatmap)
quiescence.genes = unlist(read.table("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\quiescence_genes.csv",sep=","), use.names = F)
qui.df = all.genes.df[all.genes.df$genes %in% quiescence.genes,] %>% filter(Age==TPM.Age)

qui.mat = qui.df %>% distinct()
qui.mat$log2TPM = log2(qui.mat$TPM)
qui.mat %>% pivot_wider(names_from="Age",values_from="log2TPM",id_cols=c(genes,LoopType)) %>% select(genes,Aged,Young,LoopType) %>% 
  filter(LoopType=="Reg") %>% tibble::column_to_rownames("genes") %>% select(Aged,Young) %>% as.matrix() %>% 
  pheatmap(na_col = "black", cluster_rows = F, cluster_cols = F, angle_col = 45, 
           file=file.path(projdir,"loop_quiescence_genes.png"), width=4,height=10)

# ChIPSeeker --------------------------------------------------------------

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

young.bed = readPeakFile(file.path(projdir,"young_merged_loop_anchors.bed"))
aged.bed = readPeakFile(file.path(projdir,"aged_merged_loop_anchors.bed"))

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
young.tagMatrix = getTagMatrix(young.bed, windows=promoter)
aged.tagMatrix = getTagMatrix(young.bed, windows=promoter)

tagHeatmap(aged.tagMatrix, xlim=c(-5000, 5000), color="red")

# plotAvgProf(young.tagMatrix, xlim=c(-5000, 5000),
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",
#             conf = 0.95, resample = 1000)

genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")
young.anchor.matrix_extension <- getTagMatrix(young.bed,windows = genebody, nbin = 800,
                                       upstream = 3000,downstream = 3000)
aged.anchor.matrix_extension <- getTagMatrix(aged.bed,windows = genebody, nbin = 800,
                                      upstream = 3000,downstream = 3000)
profile_plot = plotPeakProf(list(Y=young.anchor.matrix_extension,
                                 A=aged.anchor.matrix_extension),
                            conf = 0.95, ylab="Loop anchor count frequency") 
profile_plot +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)], labels=c("Y"="Young", "A"="Aged")) +
  scale_color_manual(values=pal_nejm()(2)[c(2,1)], labels=c("Y"="Young", "A"="Aged")) +
  ylab("Loop Anchor Count Frequency") +
  theme(legend.position = "top",
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=14))
ggsave(file.path(projdir,"loop_anchor_frequency_genebodies.png"), width=6, height=4, dpi=300)

young.peakAnno <- annotatePeak(young.bed, tssRegion=c(-1000, 1000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
aged.peakAnno <- annotatePeak(aged.bed, tssRegion=c(-1000, 1000),
                              TxDb=txdb, annoDb="org.Mm.eg.db")

saveRDS(young.peakAnno, file.path(projdir,"young_peakAnno.RDS"))
saveRDS(aged.peakAnno, file.path(projdir,"aged_peakAnno.RDS"))

young.peakAnno.df = data.frame(young.peakAnno@annoStat)
aged.peakAnno.df = data.frame(aged.peakAnno@annoStat)
young.peakAnno.df$Count = 0.01*young.peakAnno.df$Frequency * young.peakAnno@peakNum
aged.peakAnno.df$Count = 0.01*aged.peakAnno.df$Frequency * aged.peakAnno@peakNum
young.peakAnno.df$Age = "Y"
aged.peakAnno.df$Age = "A"
peakAnno.df = rbind(young.peakAnno.df, aged.peakAnno.df)
peakAnno.df$Feature = factor(peakAnno.df$Feature)
peakAnno.df$Feature = factor(peakAnno.df$Feature, levels=rev(levels(peakAnno.df$Feature)))

ggplot(peakAnno.df, aes(y=Age, x=Count)) +
  geom_col(aes(fill=Feature)) +
  scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,6000,1000)) +
  scale_y_discrete(expand=expansion(mult=c(0,0.1))) +
  scale_fill_manual(values=rev(brewer.pal(length(unique(peakAnno.df$Feature)),"Paired"))) +
  theme_bw() + 
  #xlab("Frequency (%)") +
  xlab("# Loop Anchors") +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.text=element_text(size=14, face="bold", color="black"),
        axis.title=element_text(size=14, face="bold", color="black"),
        legend.title=element_text(face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir, "loop_anchor_anno_counts.png"), dpi=300, width=6, height=3)

png(file.path(projdir,"young_loop_anchor_anno_pie.png"), width=6, height=4, res=300, units="in")
plotAnnoPie(young.peakAnno)
dev.off()
png(file.path(projdir,"aged_loop_anchor_anno_pie.png"), width=6, height=4, res=300, units="in")
plotAnnoPie(aged.peakAnno)
dev.off()
ChIPseeker::upsetplot(young.peakAnno, vennpie=T)

plotDistToTSS(young.peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
plotDistToTSS(aged.peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

# WebGestalt --------------------------------------------------------------

library(WebGestaltR)
