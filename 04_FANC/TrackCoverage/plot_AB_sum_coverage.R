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
projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC"

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


# Get bigwig sums ---------------------------------------------------------

files = list.files(file.path(projdir,"track_coverage","sum_coverage"), pattern="sumCoverage")

bw_df = list()
for(f in files) {
  fname = gsub(pattern=".txt",replacement="",x=basename(f))
  tmp = read.table(file.path(projdir,"track_coverage","sum_coverage",f),header=T)
  tmp$fname = fname
  #Age = ifelse(grepl(pattern="aged",taskname,fixed=T), "Aged", "Young")
  #LoopType = ifelse(grepl(pattern="diff",taskname,fixed=T), "Diff", "Reg")
  #tmp = data.frame(filename=taskname, genes=genes, Age=Age, LoopType=LoopType)
  bw_df[[fname]] = tmp
}
bw_df = bind_rows(bw_df)
bw_df = bw_df %>% separate(fname, into=c("Age","Type",NA,NA), sep="\\.", remove=F)

bw_plt_df = bw_df %>%
  pivot_longer(cols=c(H4K20me1,H3K27me3,H3K4me3,ATAC), names_to="Mark", values_to="RPKM") %>%
  mutate(meanRPKM=RPKM/(end-start),
         Age=factor(Age,levels=c("young","aged"))) %>% 
  distinct()


# Only plot A/B compartments ----------------------------------------------


ab.test = bw_plt_df %>%
  mutate(log2mean=log2(meanRPKM+1)) %>%
  filter((Type %in% c('A','B'))) %>%
  group_by(Mark,Age) %>%
  wilcox_test(log2mean~Type) %>% 
  add_xy_position(x="Type",group="Age",dodge=0.9,fun="mean_se")

ab.age.test = bw_plt_df %>%
  mutate(log2mean=log2(meanRPKM+1)) %>%
  filter((Type %in% c('A','B'))) %>%
  group_by(Type,Mark) %>%
  wilcox_test(log2mean~Age) %>% 
  add_xy_position(x="Type",group="Age",dodge=0.9,fun="mean_se")

bw_plt_df %>%
  filter((Type %in% c('A','B'))) %>%
  ggplot(aes(x=Type,y=log2(meanRPKM+1))) +
  facet_wrap(~Mark, scales="free_y") +
  stat_summary(aes(group=Age),fun.data="mean_se",geom="errorbar",width=0.5,position=position_dodge(0.9)) +
  stat_summary(aes(fill=Age),fun.data="mean_se",geom="bar",position=position_dodge(0.9)) +
  #geom_boxplot(aes(fill=Age), color="black", outlier.size=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  stat_pvalue_manual(ab.age.test, tip.length=0.01) +
  stat_pvalue_manual(ab.test, tip.length=0.01, bracket.nudge.y = 1, step.group.by = c("Age","Mark"), step.increase = 0.1) +
  labs(x=NULL) +
  theme(axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=12,face="bold",color="black"),
        strip.text=element_text(size=12,face="bold",color="black"))
ggsave(file.path(projdir,"track_coverage","sum_coverage","sum_coverage_AB.png"),dpi=300,width=7,height=5)

# Plot compartment switches -----------------------------------------------

ab.switch.test = bw_plt_df %>%
  mutate(log2mean=log2(meanRPKM+1)) %>%
  filter(!(Type %in% c('A','B'))) %>%
  group_by(Mark,Age) %>%
  wilcox_test(log2mean~Type, ref.group = "static") %>% 
  add_xy_position(x="Type",group="Age",dodge=0.9,fun = "mean_se")

ab.age.switch.test = bw_plt_df %>%
  mutate(log2mean=log2(meanRPKM+1)) %>%
  filter(!(Type %in% c('A','B'))) %>%
  group_by(Type,Mark) %>%
  wilcox_test(log2mean~Age) %>% 
  add_xy_position(x="Type",group="Age",dodge=0.9,fun="mean_se")

bw_plt_df %>%
  filter(!(Type %in% c('A','B'))) %>%
  ggplot(aes(x=Type,y=log2(meanRPKM+1))) +
  facet_wrap(~Mark, scales="free_y") +
  stat_summary(aes(group=Age),fun.data="mean_se",geom="errorbar",width=0.5,position=position_dodge(0.9)) +
  stat_summary(aes(fill=Age),fun.data="mean_se",geom="bar",position=position_dodge(0.9)) +
  #geom_boxplot(aes(fill=Age), color="black", outlier.size=0.1, position=position_dodge(0.9)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  stat_pvalue_manual(ab.age.switch.test, tip.length=0.01, label.size=3) +
  stat_pvalue_manual(ab.switch.test, tip.length=0.01, bracket.nudge.y = 0.75, label="p.adj", label.size = 3, step.group.by = c("Mark"), step.increase = 0.075) +
  labs(x=NULL) +
  theme(axis.text=element_text(size=12,face="bold",color="black"),
        axis.title=element_text(size=12,face="bold",color="black"),
        strip.text=element_text(size=12,face="bold",color="black"))
ggsave(file.path(projdir,"track_coverage","sum_coverage","sum_coverage_ABswitch.png"),dpi=300,width=7,height=5)

# Fold change -------------------------------------------------------------



bw_plt_df %>% pivot_wider(names_from=Age,values_from=RPKM,id_cols=c(chr,start,end,Mark,Type)) %>% mutate(FC=log2(aged/young+1)) %>% drop_na() %>% distinct() %>%
  ggplot(aes(x=Type,y=FC)) +
  facet_wrap(~Mark,scales = "free_y") +
  geom_boxplot()
