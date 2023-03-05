library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(rstatix)
library(VennDiagram)
library(tidyr)
library(dplyr)
library(ggpattern)

# Import data -------------------------------------------------------------

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC"

tpm.data <- readRDS(file.path(projdir, "filtered_tpm_data.RDS"))
# tpm.data = read.table(file.path(projdir, "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
# gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
# gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
# duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]
# 
# tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
# rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
# tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
# tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

DE.genes = read.table(file.path(projdir, "compartmentExpression", "DE_0.05padj_log2(1.5)LFC_genenames.tsv"), sep="\t", as.is=TRUE)
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)
DE.data = read.table(file.path(projdir, "compartmentExpression", "DE_0.05padj_log2(1.5)LFC.tsv"), sep="\t")
DE.data$symbol = sapply(rownames(DE.data), function(x) unlist(strsplit(x, split="_", fixed=TRUE))[2])

# Plot TPM per compartment ------------------------------------------------

genenames.dir = file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "tss")
#young.genes = read.table(file.path(genenames.dir, "young.tss.near.ab.bed"), sep="\t")
#aged.genes = read.table(file.path(genenames.dir, "aged.tss.near.ab.bed"), sep="\t")

#young.a.genes = young.genes %>% filter(V5=="A") %>% pull(V4) %>% unique()
#young.b.genes = young.genes %>% filter(V5=="B") %>% pull(V4) %>% unique()
#aged.a.genes = aged.genes %>% filter(V5=="A") %>% pull(V4) %>% unique()
#aged.b.genes = aged.genes %>% filter(V5=="B") %>% pull(V4) %>% unique()

young.a.genes = unlist(read.table(file.path(genenames.dir, "young.A.1kb.tss.bed"))$V4, use.names=F)
young.b.genes = unlist(read.table(file.path(genenames.dir, "young.B.1kb.tss.bed"))$V4, use.names=F)
aged.a.genes = unlist(read.table(file.path(genenames.dir, "aged.A.1kb.tss.bed"))$V4, use.names=F)
aged.b.genes = unlist(read.table(file.path(genenames.dir, "aged.B.1kb.tss.bed"))$V4, use.names=F)

young.A.tpm = tpm.data[rownames(tpm.data) %in% young.a.genes, ]
young.B.tpm = tpm.data[rownames(tpm.data) %in% young.b.genes, ]
aged.A.tpm = tpm.data[rownames(tpm.data) %in% aged.a.genes, ]
aged.B.tpm = tpm.data[rownames(tpm.data) %in% aged.b.genes, ]

young.A.tpm.avg = rowMeans(young.A.tpm[,5:7])
young.B.tpm.avg = rowMeans(young.B.tpm[,5:7])
aged.A.tpm.avg = rowMeans(aged.A.tpm[,1:3])
aged.B.tpm.avg = rowMeans(aged.B.tpm[,1:3])

tpm.avg.list = list(young.A.tpm.avg, young.B.tpm.avg, aged.A.tpm.avg, aged.B.tpm.avg)
names(tpm.avg.list) = c("young.A","young.B","aged.A","aged.B")
tpm.avg.list = melt(tpm.avg.list)
colnames(tpm.avg.list) = c("tpm","group")

#tpm.avg.list = tpm.avg.list[tpm.avg.list$tpm>0,]
#tpm.avg.list[tpm.avg.list$tpm<=0,"tpm"] = 0.001 # keep only genes with non-negative TPM  

tpm.avg.list = tpm.avg.list %>% 
  separate(group, into=c("Age","compartment"), remove=F) %>%
  mutate(log2tpm = log2(tpm)) %>%
  mutate_if(is.character, as.factor)
levels(tpm.avg.list$Age) = list(Y="young", A="aged")

tpm.avg.list.test = tpm.avg.list %>%
  group_by(compartment) %>%
  wilcox_test(log2tpm ~ Age) %>%
  add_xy_position(x="compartment", dodge=1)
tpm.avg.list.age.test = tpm.avg.list %>%
  group_by(Age) %>%
  wilcox_test(log2tpm ~ compartment) %>%
  add_xy_position(x="compartment", group="Age", dodge=1)
tpm.avg.list.age.test$y.position = c(17,20)

tpm.avg.list %>%
  group_by(compartment, Age) %>%
  summarise(avg=mean(tpm),
            count=n())

#png(file.path(projdir, "compartmentExpression", "100kb", "100kb_compartment_expression_minTPM0.png"), res=300, units="in", width=3, height=4.5)
ggplot(tpm.avg.list, aes(x=compartment, y=log2tpm)) +
  geom_violin(aes(fill=Age), color="black", position=position_dodge(1)) +
  geom_boxplot(aes(group=interaction(Age,compartment)), color="black", width=0.2, outlier.shape=NA, fill="white", position=position_dodge(1)) +
  stat_pvalue_manual(tpm.avg.list.test, tip.length = 0.01, label.size = 4, bracket.nudge.y = 1) +
  stat_pvalue_manual(tpm.avg.list.age.test, tip.length = 0.01, label.size = 4, bracket.nudge.y = 1) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-15,25,5)) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw() +
  ylab(expression(bold(paste("lo","g"["2"],"(TPM)")))) + xlab(NULL) +
  theme(legend.position = "top",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
#dev.off()
ggsave(file.path(projdir, "compartmentExpression", "100kb", "100kb_compartment_expression_minTPM0.png"), width=3, height=4, dpi=300)

# Plot TPM per compartment switch -----------------------------------------

#a.to.b.genes = unlist(read.table(file.path(genenames.dir, "A_to_B.500bp.promoter.genenames.txt")), use.names=F)
#b.to.a.genes = unlist(read.table(file.path(genenames.dir, "B_to_A.500bp.promoter.genenames.txt")), use.names=F)
#static.genes = unlist(read.table(file.path(genenames.dir, "static.500bp.promoter.genenames.txt")), use.names=F)

a.to.b.genes = unlist(read.table(file.path(genenames.dir, "A_to_B.1kb.tss.bed"))$V4, use.names=F)
b.to.a.genes = unlist(read.table(file.path(genenames.dir, "B_to_A.1kb.tss.bed"))$V4, use.names=F)
static.A.genes = unlist(read.table(file.path(genenames.dir, "static.A.1kb.tss.bed"))$V4, use.names=F)
static.B.genes = unlist(read.table(file.path(genenames.dir, "static.B.1kb.tss.bed"))$V4, use.names=F)

#static.genes = static.genes[!((static.genes %in% a.to.b.genes) | (static.genes %in% b.to.a.genes))]

venn.diagram(x=list(a.to.b.genes, b.to.a.genes, static.genes),
             category.names = c("A to B","B to A","Static"),
             filename=file.path(projdir,"compartmentExpression","100kb","gene.venndiagram.png"),
             output=F,
             
             imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = pal_jama()(3),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 180),
             cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "sans",
             rotation = 1)

a.to.b.tpm = tpm.data[rownames(tpm.data) %in% a.to.b.genes, ]
b.to.a.tpm = tpm.data[rownames(tpm.data) %in% b.to.a.genes, ]
static.A.tpm = tpm.data[rownames(tpm.data) %in% static.A.genes, ]
static.B.tpm = tpm.data[rownames(tpm.data) %in% static.B.genes, ]

venn.diagram(x=list(rownames(a.to.b.tpm), rownames(b.to.a.tpm), rownames(static.tpm)),
             category.names = c("A to B","B to A","Static"),
             filename=file.path(projdir,"compartmentExpression","100kb","expressed.genes.venndiagram.png"),
             output=F,
             
             imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             
             # Circles
             lwd = 1,
             lty = 'blank',
             fill = pal_jama()(3),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 180),
             cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "sans",
             rotation = 1)

a.to.b.tpm.avg = data.frame(A = rowMeans(a.to.b.tpm[,1:3]), Y = rowMeans(a.to.b.tpm[,5:7]))
b.to.a.tpm.avg = data.frame(A = rowMeans(b.to.a.tpm[,1:3]), Y = rowMeans(b.to.a.tpm[,5:7]))
static.A.tpm.avg = data.frame(A = rowMeans(static.A.tpm[,1:3]), Y = rowMeans(static.A.tpm[,5:7]))
static.B.tpm.avg = data.frame(A = rowMeans(static.B.tpm[,1:3]), Y = rowMeans(static.B.tpm[,5:7]))

#a.to.b.tpm.avg = subset(a.to.b.tpm.avg, A>0 & Y>0)
#b.to.a.tpm.avg = subset(b.to.a.tpm.avg, A>0 & Y>0)
#static.tpm.avg = subset(static.tpm.avg, A>0 & Y>0)

# save each TPM list
write.table(a.to.b.tpm[,c(1:3,5:7)] %>% tibble::rownames_to_column("Gene_Symbol"),
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "a_to_b_tpm.tsv"),
            quote = F, sep="\t", col.names=T, row.names=F)
write.table(static.A.tpm[,c(1:3,5:7)] %>% tibble::rownames_to_column("Gene_Symbol"),
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "static_A_tpm.tsv"),
            quote = F, sep="\t", col.names=T, row.names=F)

# save each avg list
write.table(a.to.b.tpm.avg,
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "a_to_b_tpm_avg.tsv"),
            quote = F, sep="\t", col.names=NA)
write.table(b.to.a.tpm.avg,
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "b_to_a_tpm_avg.tsv"),
            quote = F, sep="\t", col.names=NA)
write.table(static.A.tpm.avg,
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "staticA_tpm_avg.tsv"),
            quote = F, sep="\t", col.names=NA)
write.table(static.B.tpm.avg,
            file.path(projdir, "compartmentExpression", "compartmentSwitchTPM", "staticB_tpm_avg.tsv"),
            quote = F, sep="\t", col.names=NA)

writeClipboard(rownames(static.B.tpm)[apply(static.B.tpm, 1, function(x) any(x)>=1)])

tpm.avg.df = rbind(log2(static.A.tpm.avg),
                   log2(static.B.tpm.avg),
                   log2(b.to.a.tpm.avg),
                   log2(a.to.b.tpm.avg))
tpm.avg.df = rbind(static.A.tpm.avg,
                   static.B.tpm.avg,
                   b.to.a.tpm.avg,
                   a.to.b.tpm.avg)
tpm.avg.df$group = factor(c(rep("Static.A", nrow(static.A.tpm.avg)),
                            rep("Static.B", nrow(static.B.tpm.avg)),
                            rep("B.to.A", nrow(b.to.a.tpm.avg)),
                            rep("A.to.B", nrow(a.to.b.tpm.avg))),
                          levels=c("Static.A","Static.B","B.to.A","A.to.B"))
  
tpm.avg.df[apply(tpm.avg.df, 1, function(x) any(x>=1)),] %>% dplyr::filter(group=="B.to.A") %>% mutate(log2fc = log2(Y/A),
                                                                                                       zscore = scale(Y/A, center=T, scale=T)) %>% 
  dplyr::select(zscore) %>% write.table(file="clipboard", sep='\t', row.names=T, col.names=F, quote=F)

tpm.avg.df[apply(tpm.avg.df, 1, function(x) any(x>=1)),] %>% dplyr::filter(group=="Static.B") %>% mutate(log2fc = log2((Y/A)+1),
                                                                                                         zscore = scale(Y/A, center=T, scale=T)) %>%
  dplyr::filter(is.na(log2fc)) %>% head()

tpm.avg.melt = melt(tpm.avg.df, variable.name = "Age", value.name = "tpm")
tpm.avg.melt$Age = factor(tpm.avg.melt$Age, levels=c("Y","A"))
tpm.avg.melt = tpm.avg.melt[complete.cases(tpm.avg.melt),]

tpm.facet.age.test = tpm.avg.melt %>%
  group_by(Age) %>%
  wilcox_test(tpm ~ group) %>%
  add_xy_position(x="group")
  #add_xy_position(x="group",group="Age",dodge=0.9) 
#tpm.avg.df.age.test[tpm.avg.df.age.test$Age=="A","y.position"] = tpm.avg.df.age.test[tpm.avg.df.age.test$Age=="A","y.position"] + 6

tpm.avg.melt %>%
  group_by(group, Age) %>%
  summarise(avg = mean(tpm),
            count = n())

tpm.avg.melt %>%
  ggplot(aes(x=group, y=log2(tpm))) +
  facet_grid(cols=vars(Age)) +
  geom_violin(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,group)), color="black", width=0.2, outlier.shape=NA, fill="white", position=position_dodge(0.9)) +
  #stat_pvalue_manual(tpm.facet.age.test, tip.length = 0.01, label.size = 4, label="p.adj", bracket.nudge.y = 0.5, step.increase = 0.075, step.group.by = "Age") +
  scale_x_discrete(labels=c(TeX("Static A",bold=T),TeX("Static B",bold=T),TeX("$B \\rightarrow A$", bold=T),TeX("$A \\rightarrow B$", bold=T))) +
  #scale_x_discrete(labels=c(expression(bold("Static")),expression(bold(B %->% A), cex),expression(bold(A %->% B)))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-10,30,5)) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw(base_size=14) +
  ylab(expression(bold(paste("lo","g"["2"],"(TPM)")))) + xlab(NULL) +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        strip.text.x = element_blank() , 
        strip.background = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ),
        axis.text = element_text(face="bold", color="black"),
        axis.title = element_text(face="bold"))
ggsave(file.path(projdir, "compartmentExpression", "100kb", "100kb_compartment_switch_expression_facet_minTPM0.png"), dpi=300, width=6, height=4)

tpm.age.test = tpm.avg.melt %>%
  group_by(group) %>%
  mutate(log2tpm = log2(tpm)) %>%
  wilcox_test(log2tpm ~ Age) %>%
  add_xy_position(x="group",group="Age",dodge=0.9) %>%
  mutate(p.label = ifelse(p<2.2e-16, "<2.2e-16", p))

tpm.avg.melt %>%
  ggplot(aes(x=group, y=log2(tpm))) +
  geom_violin(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,group)), color="black", width=0.2, outlier.shape=NA, fill="white", position=position_dodge(0.9)) +
  stat_pvalue_manual(tpm.age.test, label = "p.label", tip.length = 0.01, label.size = 4, bracket.nudge.y = 0.05) +
  scale_x_discrete(labels=c(TeX("Static A",bold=T),TeX("Static B",bold=T),TeX("$B \\rightarrow A$", bold=T),TeX("$A \\rightarrow B$", bold=T))) +
  #scale_x_discrete(labels=c(expression(bold("Static")),expression(bold(B %->% A), cex),expression(bold(A %->% B)))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-10,30,5)) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)], labels=c("Y"="Young", "A"="Aged")) +
  theme_bw() +
  ylab(expression(bold(paste("lo","g"["2"],"(TPM)")))) + xlab(NULL) +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        strip.text.x = element_blank() , 
        strip.background = element_blank(),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title.y = element_text(size=12, face="bold", margin = margin(t = 0, r = 0, b = 0, l = 0)))
ggsave(file.path(projdir, "compartmentExpression", "100kb", "100kb_compartment_switch_expression_minTPM0.png"), dpi=300, width=3.5, height=4)


# Export genes for GSEA and ORA -------------------------------------------

# no results for log2fc or zscore
a.to.b.tpm.avg %>% mutate(fc = A/Y, log2fc = log2(A/Y), zscore = scale(A/Y)) %>% 
  dplyr::select(fc) %>%
  write.table(file="clipboard", col.names=F, quote=F, sep='\t')
rownames(a.to.b.tpm.avg)[rownames(a.to.b.tpm.avg) %in% DE.genes.name]
writeClipboard(rownames(a.to.b.tpm.avg))

b.to.a.tpm.avg %>% mutate(fc = A/Y, log2fc = log2(A/Y), zscore = scale(A/Y)) %>% 
  dplyr::select(fc) %>%
  write.table(file="clipboard", col.names=F, quote=F, sep='\t')

static.A.tpm.avg %>% 
  dplyr::filter(A>=1 | Y>=1) %>%
  mutate(fc = A/Y, log2fc = log2(A/Y), zscore = scale(A/Y)) %>% 
  dplyr::select(log2fc) %>%
  write.table(file="clipboard-1000", col.names=F, quote=F, sep='\t')

static.A.tpm.avg %>% 
  dplyr::filter(A>=1 | Y>=1) %>%
  mutate(fc = A/Y, log2fc = log2(A/Y), zscore = scale(A/Y)) %>% 
  dplyr::select(zscore) %>%
  write.table(file="clipboard-1000", col.names=F, quote=F, sep='\t')
writeClipboard(rownames(static.A.tpm)[rownames(static.A.tpm) %in% DE.genes.name])

static.A.hmp = pheatmap(static.A.tpm[rownames(static.A.tpm) %in% DE.genes.name, c(1:3,5:7)], cluster_cols=F, scale="row", kmeans_k = 6)
writeClipboard(names(static.A.hmp$kmeans$cluster[static.A.hmp$kmeans$cluster!=6]))
sort(cutree(static.A.hmp$tree_row, k=2))
plot(static.A.hmp$tree_row)

static.B.tpm.avg %>% mutate(fc = A/Y, log2fc = log2(A/Y), zscore = scale(A/Y)) %>% 
  dplyr::select(fc) %>%
  write.table(file="clipboard-1000", col.names=F, quote=F, sep='\t')

## signal to noise metric for GSEA ------------------------------------------

static_A_s2n = data.frame(aged_avg = rowMeans(static.A.tpm[,1:3]),
                          young_avg = rowMeans(static.A.tpm[,5:7]),
                          aged_sd = matrixStats::rowSds(as.matrix(static.A.tpm[,1:3])),
                          young_sd = matrixStats::rowSds(as.matrix(static.A.tpm[,5:7])))
static_A_s2n %>% 
  mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
  dplyr::select(s2n) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  #dplyr::filter(gene %in% DE.genes.name) %>%
  write.table(file="clipboard-1000", sep='\t', col.names=F, row.names=F, quote=F)

static_B_s2n = data.frame(aged_avg = rowMeans(static.B.tpm[,1:3]),
                          young_avg = rowMeans(static.B.tpm[,5:7]),
                          aged_sd = matrixStats::rowSds(as.matrix(static.B.tpm[,1:3])),
                          young_sd = matrixStats::rowSds(as.matrix(static.B.tpm[,5:7])))
static_B_s2n %>% 
  mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
  dplyr::select(s2n) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  #dplyr::filter(gene %in% DE.genes.name) %>%
  write.table(file="clipboard-1000", sep='\t', col.names=F, row.names=F, quote=F)

b.to.a_s2n = data.frame(aged_avg = rowMeans(b.to.a.tpm[,1:3]),
                        young_avg = rowMeans(b.to.a.tpm[,5:7]),
                        aged_sd = rowSds(as.matrix(b.to.a.tpm[,1:3])),
                        young_sd = rowSds(as.matrix(b.to.a.tpm[,5:7])))
b.to.a_s2n %>% 
  mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
  dplyr::select(s2n) %>%
  as.data.frame() %>%
  write.table(file="clipboard-1000", sep='\t', col.names=F, row.names=T, quote=F)

static_A_ttest = apply(static.A.tpm[,c(1:3,5:7)], 1, function(x) t.test(x[1:3], x[4:6])$statistic)
static_A_ttest_df = data.frame(ttest=static_A_ttest) %>% tibble::rownames_to_column("gene")
static_A_ttest_df = static_A_ttest_df[static_A_ttest_df$gene %in% DE.genes.name,]
write.table(static_A_ttest_df, file="clipboard-1000", sep='\t', col.names=F, row.names=F, quote=F)

static_B_ttest = apply(static.B.tpm[,c(1:3,5:7)], 1, function(x) t.test(x[1:3], x[4:6])$statistic)
static_B_ttest_df = as.data.frame(static_B_ttest)
write.table(static_B_ttest_df, file="clipboard-1000", sep='\t', col.names=F, row.names=T, quote=F)

# Intersection between compartment switching genes and gene lists ---------

library(readxl)
genesets = read_xlsx("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Sestrins\\Mouse RNA-Seq\\TargetedGeneLists.xlsx")

quiescence_genes = genesets %>% pull(Quiescence) %>% unique()
cellcycle_genes = genesets %>% pull(Cell_Cycle) %>% unique()

length(intersect(quiescence_genes, static.A.genes)) / length(quiescence_genes)
length(intersect(quiescence_genes, static.B.genes)) / length(quiescence_genes)

length(intersect(cellcycle_genes, static.A.genes)) / length(cellcycle_genes)
length(intersect(cellcycle_genes, static.B.genes)) / length(cellcycle_genes)


# of DE genes per compartment group -------------------------------------

table(DE=rownames(static.A.tpm.avg) %in% DE.genes.name, up_in_aged=with(static.A.tpm.avg, log2(A/Y)>0))
table(DE=rownames(static.B.tpm.avg) %in% DE.genes.name, up_in_aged=with(static.B.tpm.avg, log2(A/Y)>0))
table(rownames(a.to.b.tpm.avg) %in% DE.genes.name)
table(rownames(b.to.a.tpm.avg) %in% DE.genes.name)

# Try a heatmap -----------------------------------------------------------

library(pheatmap)
tpm.df = rbind(log2(static.A.tpm),
                   log2(static.B.tpm),
                   log2(b.to.a.tpm),
                   log2(a.to.b.tpm))
tpm.avg.df$group = factor(c(rep("Static.A", nrow(static.A.tpm.avg)),
                            rep("Static.B", nrow(static.B.tpm.avg)),
                            rep("B.to.A", nrow(b.to.a.tpm.avg)),
                            rep("A.to.B", nrow(a.to.b.tpm.avg))),
                          levels=c("Static.A","Static.B","B.to.A","A.to.B"))
tpm.avg.mat = as.matrix(tpm.avg.df[,1:2])
pheatmap(tpm.avg.mat, cluster_cols = F, show_rownames = F)
library(ComplexHeatmap)
hmp.plt = Heatmap(tpm.avg.mat, cluster_columns = F, row_split = tpm.avg.df$group, show_row_names = F)
draw(hmp.plt)

# DE genes ----------------------------------------------------------------

a.to.b.DE = a.to.b.genes[a.to.b.genes %in% DE.genes.name]
b.to.a.DE = b.to.a.genes[b.to.a.genes %in% DE.genes.name]
static.A.DE = static.A.genes[static.A.genes %in% DE.genes.name]
static.B.DE = static.B.genes[static.B.genes %in% DE.genes.name]

a.to.b.DE.tpm = tpm.data[rownames(tpm.data) %in% a.to.b.DE, ]
b.to.a.DE.tpm = tpm.data[rownames(tpm.data) %in% b.to.a.DE, ]
static.A.DE.tpm = tpm.data[rownames(tpm.data) %in% static.A.DE, ]
static.B.DE.tpm = tpm.data[rownames(tpm.data) %in% static.B.DE, ]

a.to.b.DE.tpm.avg = data.frame(A = rowMeans(a.to.b.DE.tpm[,1:3]), Y = rowMeans(a.to.b.DE.tpm[,4:8]))
b.to.a.DE.tpm.avg = data.frame(A = rowMeans(b.to.a.DE.tpm[,1:3]), Y = rowMeans(b.to.a.DE.tpm[,4:8]))
static.A.DE.tpm.avg = data.frame(A = rowMeans(static.A.DE.tpm[,1:3]), Y = rowMeans(static.A.DE.tpm[,4:8]))
static.B.DE.tpm.avg = data.frame(A = rowMeans(static.B.DE.tpm[,1:3]), Y = rowMeans(static.B.DE.tpm[,4:8]))

tpm.avg.df = rbind(static.A.DE.tpm.avg,
                   static.B.DE.tpm.avg,
                   b.to.a.DE.tpm.avg,
                   a.to.b.DE.tpm.avg)
tpm.avg.df$group = factor(c(rep("Static.A", nrow(static.A.DE.tpm.avg)),
                            rep("Static.B", nrow(static.B.DE.tpm.avg)),
                            rep("B.to.A", nrow(b.to.a.DE.tpm.avg)),
                            rep("A.to.B", nrow(a.to.b.DE.tpm.avg))),
                          levels=c("Static.A","Static.B","B.to.A","A.to.B"))
#tpm.avg.df = subset(tpm.avg.df, A>0 & Y>0)

tpm.avg.df[apply(tpm.avg.df, 1, function(x) any(x>=1)),] %>% dplyr::filter(group=="Static.A") %>% mutate(log2fc = log2(Y/A),
                                                                                                       zscore = scale(Y/A, center=T, scale=T)) %>% 
  dplyr::select(zscore) %>% write.table(file="clipboard-1000", sep='\t', row.names=T, col.names=F, quote=F)


tpm.avg.df[,c("A","Y")] = data.frame(apply(tpm.avg.df[,c("A","Y")], 2, log2))


tpm.avg.melt = melt(tpm.avg.df, variable.name = "Age", value.name = "log2tpm")
tpm.avg.melt$Age = factor(tpm.avg.melt$Age, levels=c("Y","A"))
tpm.avg.melt = tpm.avg.melt[complete.cases(tpm.avg.melt),]

DE.tpm.avg.df.age.test = tpm.avg.melt %>%
  group_by(Age) %>%
  wilcox_test(log2tpm ~ group) %>%
  add_xy_position(x="group")

tpm.avg.melt %>%
  group_by(group, Age) %>%
  summarise(avg = mean(log2tpm),
            count = n(),
            summary = list(quantile(log2tpm)))

venn.diagram(x=list(rownames(a.to.b.DE.tpm.avg), rownames(b.to.a.DE.tpm.avg), rownames(static.DE.tpm.avg)),
             category.names = c("A to B","B to A","Static"),
             filename=file.path(projdir,"compartmentExpression","100kb","DE.genes.venndiagram.png"),
             output=F,
             
             imagetype="png" ,
             height = 600 , 
             width = 600 , 
             resolution = 300,
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = pal_jama()(3),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 180),
             cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "sans",
             rotation = 1)

#png(file.path(projdir, "compartmentExpression", "compartment_switch_expression_0.05padj_1LFC_DE.png"), res=200, units="in", width=2.5, height=4)
tpm.avg.melt %>%
  ggplot(aes(x=group, y=log2tpm)) +
  facet_grid(cols=vars(Age)) +
  geom_violin(aes(fill=Age), color="black") +
  geom_boxplot(aes(group=interaction(group,Age)), color="black", outlier.shape=NA, position=position_dodge(0.9), width=0.2) +
  #geom_hline(yintercept=median(with(tpm.avg.melt,log2tpm[group=="Static" & Age=="Y"])), linetype="dashed") +
  stat_pvalue_manual(DE.tpm.avg.df.age.test, tip.length = 0.01, label.size = 4, bracket.nudge.y = 0.25, step.group.by = "Age", step.increase = 0.05) +
  scale_x_discrete(labels=c(TeX("Static",bold=T),TeX("$B \\rightarrow A$", bold=T),TeX("$A \\rightarrow B$", bold=T))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-5,20,2.5)) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw(base_size=14) +
  ylab(expression(bold(paste("lo","g"["2"],"(TPM)")))) + xlab(NULL) +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        plot.margin = unit( c(0,0,0,0) , units = "lines" ),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(face="bold"))
#dev.off()
ggsave(file.path(projdir, "compartmentExpression", "100kb", "compartment_switch_expression_0.05padj_log2(1.5)LFC_DE.png"), dpi=300, width=4.5, height=4)

# Fold change -------------------------------------------------------------

a.to.b.tpm.fc = with(a.to.b.tpm.avg, log2(A/Y))
b.to.a.tpm.fc = with(b.to.a.tpm.avg, log2(A/Y))
static.A.tpm.fc = with(static.A.tpm.avg, log2(A/Y))

tpm.fc.df = data.frame(group = c(rep("Static", length(static.tpm.fc)),
                                 rep("B.to.A", length(b.to.a.tpm.fc)),
                                 rep("A.to.B", length(a.to.b.tpm.fc))),
                       fc = c(static.tpm.fc, b.to.a.tpm.fc, a.to.b.tpm.fc))
tpm.fc.df$group = factor(tpm.fc.df$group, levels=c("Static", "B.to.A", "A.to.B"))

tpm.fc.df.test = tpm.fc.df %>%
  t_test(fc ~ group, ref.group="Static") %>%
  add_xy_position(x="group")

tpm.fc.df %>%
  group_by(group) %>%
  summarise(mean(fc, na.rm=T))

png(file.path(projdir, "compartmentExpression", "compartment_switch_fc.png"), res=200, units="in", width=2.5, height=4)
ggplot(tpm.fc.df, aes(x=group, y=fc)) +
  geom_boxplot(aes(fill=group), color="black", outlier.size = 0.1) +
  #stat_pvalue_manual(tpm.fc.df.test, tip.length = 0.01, label.size = 5) +
  scale_fill_jama() +
  theme_bw() + 
  xlab("") + ylab("log2(Aged TPM / Young TPM)") +
  theme(legend.position = "none",
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
dev.off()