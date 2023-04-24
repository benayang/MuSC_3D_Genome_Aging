library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(rstatix)
library(tidyr)
library(scales)
library(dplyr)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C"


# Get expression data -----------------------------------------------------

tpm.data = read.table(file.path(projdir, "04_FANC", "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with non-negative TPM  


# Plot coverage -----------------------------------------------------------

files = list.files(file.path(projdir,"08_HiCExplorer","TAD histone overlap","coverage"), pattern = ".tsv")
fnames = sapply(files, function(x) unlist(strsplit(x,split=".",fixed=T)))
fnames = sapply(fnames, function(x) paste(head(x[-c(1:4)], -1), collapse = '.'))

coverage.list = list()
for(f in files) {
  tmp = read.table(file.path(projdir,"08_HiCExplorer","TAD histone overlap","coverage",f), header=T, sep='\t')
  tmp$Group = fnames[f]
  tmp = tmp[complete.cases(tmp),]
  coverage.list[[fnames[f]]] = tmp
}

all.coverage = do.call(rbind, coverage.list)
all.coverage.log1p = all.coverage
all.coverage.log1p[,4:12] = apply(all.coverage.log1p[,4:12], 2, function(x) log10(x+1))

all.coverage.log1p.long = all.coverage.log1p %>%
  pivot_longer(cols=4:12, names_to="histone", values_to="coverage")

all.coverage.log1p.long %>%
  group_by(Group, histone) %>%
  summarise(count = n()) %>%
  as.data.frame()

coverage.test = all.coverage.log1p.long %>%
  #filter(Group %in% c("A_to_B","B_to_A","static")) %>%
  group_by(histone) %>%
  t_test(coverage ~ Group) %>%
  add_xy_position(x="histone", group = "Group")

all.coverage.log1p.long %>%
  filter(Group %in% c("A_to_B","B_to_A","static") & coverage>0) %>%
  ggplot(aes(x=histone, y=coverage)) +
  geom_boxplot(aes(fill=Group), color="black", outlier.size=0.1) + 
  #stat_pvalue_manual(coverage.test, tip.length = 0.01, bracket.nudge.y = 0.5, step.increase = 0.1, step.group.by = "histone") +
  scale_fill_npg() +
  theme_bw() +
  ylab(expression(bold(paste("lo","g"["10"],"(Signal + 1)")))) + xlab("") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"08_HiCExplorer","TAD histone overlap","coverage","compartments_h3k27me3_coverage_allgenes.png"), width=2.5, height=4.5, dpi=300)

# Coverage in expressed genes ---------------------------------------------

expressed.genes = rownames(tpm.data)
sum(all.coverage.log1p.long$Gene %in% expressed.genes)/nrow(all.coverage.log1p.long)

exp.all.coverage.log1p = all.coverage.log1p.long %>%
  filter(Gene %in% expressed.genes & Group %in% c("A_to_B","B_to_A","static")) 

expressed.test = exp.all.coverage.log1p %>%
  filter(Group %in% c("A_to_B","B_to_A","static") & coverage>0 & histone %in% c("youngGSEH3K27me3", "youngRPKMH3K27me3", "agedRPKMH3K27me3", "youngH4K20me1", "agedH4K20me1")) %>%
  mutate(histone = factor(histone, levels=c("youngGSEH3K27me3", "youngH3K27me3", "agedH3K27me3", "youngH4K20me1", "agedH4K20me1"))) %>%
  group_by(histone) %>%
  t_test(coverage ~ Group) %>%
  add_xy_position(x="histone", group="Group")

exp.all.coverage.log1p %>%
  filter(Group %in% c("A_to_B","B_to_A","static") & coverage>0 & histone %in% c("youngGSEH3K27me3", "youngRPKMH3K27me3", "agedRPKMH3K27me3", "youngH4K20me1", "agedH4K20me1")) %>%
  mutate(histone = factor(histone, levels=c("youngGSEH3K27me3", "youngH3K27me3", "agedH3K27me3", "youngH4K20me1", "agedH4K20me1"))) %>%
  ggplot(aes(x=histone, y=coverage, color=Group)) +
  geom_boxplot(aes(fill=Group), color="black", outlier.size=0.1) + 
  #geom_dotplot(binaxis="y", stackdir="center", dotsize=0.1, binwidth = 0.1, position=position_dodge(0.9)) + 
  stat_pvalue_manual(expressed.test, tip.length=0.01) +
  scale_fill_npg() +
  theme_bw() +
  ylab(expression(bold(paste("lo","g"["10"],"(Signal + 1)")))) + xlab("") +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        axis.text.x = element_text(angle=15, hjust=1))
ggsave(file.path(projdir,"08_HiCExplorer","TAD histone overlap","coverage","compartments_histone_coverage_minTPM0.png"), width=7, height=5, dpi=300)
