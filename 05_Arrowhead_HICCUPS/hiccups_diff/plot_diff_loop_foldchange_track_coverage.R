library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
library(rstatix)
library(ggsci)
library(ggpubr)
library(latex2exp)
library(HelloRanges)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\hiccups_diff\\loop_domain_track_coverage"

tasks=list.files(projdir, pattern=".tab")

data=list()
for(taskname in tasks){
  print(taskname)
  tmp = read.table(file.path(projdir,taskname),header=F,col.names=c("chr","start","end","ATAC","H3K4me3"))
  fname_list = unlist(strsplit(taskname, split="_", fixed=T))
  fname = gsub(".tab", "", paste(fname_list[4:length(fname_list)], collapse="_"))
  tmp$Group = taskname
  tmp$Diff = ifelse(grepl("non", fname, fixed=T), "Non-Diff", "Diff")
  tmp$loopAge = ifelse(grepl("young", fname, fixed=T), "young", "aged")
  tmp$markAge = fname_list[[1]]
  data[[taskname]] = tmp
  print(head(tmp))
}

df=bind_rows(data)


# Look at compartment switch ----------------------------------------------

comp.switch.test = df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  mutate(loopAge=factor(loopAge,levels=c("young","aged"))) %>%
  group_by(loopAge,BW) %>%
  wilcox_test(FC~Diff) %>%
  add_xy_position(x="loopAge",group="Diff",dodge=0.9,fun="mean_se")

comp.switch.win.age.test = df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  group_by(Group,BW) %>%
  wilcox_test(FC~loopAge) %>%
  add_xy_position(x="Diff",group="loopAge",fun="mean_se",dodge=0.9)

df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  group_by(Group,BW,Age) %>%
  summarise(count=n())

df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  dplyr::filter(loopAge==markAge) %>%
  mutate(loopAge=factor(loopAge,levels=c("young","aged"))) %>%
  group_by(loopAge,BW,Diff) %>%
  summarise(FC_avg = mean(FC), FC_sem = sd(FC)/sqrt(n())) %>%
  ggplot(aes(x=loopAge,y=FC_avg, group=Diff)) + 
  facet_grid(rows=vars(BW), scales="free_y") +
  #stat_summary(aes(group=Diff), fun.data="mean_se",geom="errorbar", width=0.5, position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=FC_avg-FC_sem, ymax=FC_avg+FC_sem), width=0.5, position=position_dodge(0.9)) +
  geom_col_pattern(aes(fill=loopAge, pattern=Diff), 
                   pattern_fill="black", 
                   pattern_color="black",
                   pattern_density = 0.1,
                   color="black", 
                   position=position_dodge(0.9)) +
  scale_fill_manual(values=rev(pal_nejm()(2)), labels=c("Y","A")) +
  scale_pattern_manual(values=c("stripe","none"), labels=c("Unique", "Shared")) +
  stat_pvalue_manual(comp.switch.test, tip.length=0.1, label.size=5, bracket.nudge.y = 0.1) +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_y_continuous(expand = expansion(mult=c(0,0.2))) +
  theme_bw() + guides(fill="none", pattern=guide_legend(title="Loop", override.aes=list(fill="white"))) +
  labs(x=NULL, y="Fold Change") +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.margin = margin(t=0,r=0,b=0,l=0,unit="pt"),
        legend.position = "top",
        legend.text=element_text(size=13),
        legend.title=element_text(size=13, face="bold"))
ggsave(file.path(projdir,"loop_domain_track_fold_change.png"),dpi=300,width=3,height=4)


# Gene expression in diff vs non-diff loops -------------------------------


tpm.data = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression", 
                                "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

avg.tpm.data = data.frame(A = rowMeans(tpm.data[,1:3]), Y = rowMeans(tpm.data[,5:7])) %>% tibble::rownames_to_column("gene")

DE.genes = unlist(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression", 
                                       "DE_0.05padj_log2(1.5)LFC_genenames.csv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)

# get 1kb TSS promoters
tss_slop_1kb = GRanges(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                            "TAD expression","tss_1kb_slop.bed"), col.names = c("chr","start","end","gene")))

diff_aged_loops = df %>% dplyr::filter(loopAge=="aged" & Diff=="Diff") %>% dplyr::select(chr, start, end) %>% GRanges()
diff_young_loops = df %>% dplyr::filter(loopAge=="young" & Diff=="Diff") %>% dplyr::select(chr, start, end) %>% GRanges()
non_diff_young_loops = df %>% dplyr::filter(loopAge=="young" & Diff=="Non-Diff") %>% dplyr::select(chr, start, end) %>% GRanges()
non_diff_aged_loops = df %>% dplyr::filter(loopAge=="aged" & Diff=="Non-Diff") %>% dplyr::select(chr, start, end) %>% GRanges()

diff_aged_loop_genes = unique(eval(R_bedtools_intersect(tss_slop_1kb, diff_aged_loops, wa=T))) %>% as_tibble() %>% mutate(Diff = "Diff", Age="A")
diff_young_loop_genes = unique(eval(R_bedtools_intersect(tss_slop_1kb, diff_young_loops, wa=T))) %>% as_tibble() %>% mutate(Diff = "Diff", Age="Y")
non_diff_aged_loop_genes = unique(eval(R_bedtools_intersect(tss_slop_1kb, non_diff_aged_loops, wa=T))) %>% as_tibble() %>% mutate(Diff = "Non-Diff", Age="A")
non_diff_young_loop_genes = unique(eval(R_bedtools_intersect(tss_slop_1kb, non_diff_young_loops, wa=T))) %>% as_tibble() %>% mutate(Diff = "Non-Diff", Age="Y")

gene_df = bind_rows(diff_aged_loop_genes, diff_young_loop_genes, non_diff_aged_loop_genes, non_diff_young_loop_genes)
gene_df = left_join(gene_df, avg.tpm.data, by="gene") %>% drop_na()

gene_df_test = gene_df %>%
  pivot_longer(cols=c("A","Y"),names_to="TPM.Age",values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(log2tpm = log2(TPM),
         Age = factor(Age, levels=c("Y","A"))) %>%
  group_by(Age) %>%
  wilcox_test(log2tpm ~ Diff) %>%
  add_xy_position(x="Age", group="Diff", dodge=0.9)

gene_df_age_test = gene_df %>%
  pivot_longer(cols=c("A","Y"),names_to="TPM.Age",values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(log2tpm = log2(TPM),
         Age = factor(Age, levels=c("Y","A"))) %>%
  group_by(Diff) %>%
  wilcox_test(log2tpm ~ Age) %>%
  add_xy_position(x="Age", group="Diff", dodge=0.9)

gene_df %>%
  pivot_longer(cols=c("A","Y"),names_to="TPM.Age",values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(Age = factor(Age, levels=c("Y","A"))) %>%
  ggplot(aes(x=Age, y=log2(TPM))) +
  geom_violin_pattern(aes(group=interaction(Age,Diff), fill=Age, pattern=Diff), 
                      color="black",
                      pattern_fill = "black",
                      pattern_color = "black",
                      pattern_density = 0.1,
                      position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,Diff)), width=0.2, color="black", fill="white", outlier.shape=NA, position=position_dodge(0.9)) +
  stat_pvalue_manual(gene_df_test, tip.length=0.01, position=position_dodge(0.9), bracket.nudge.y = 0.5, size=4) +
  stat_pvalue_manual(gene_df_age_test, tip.length=0.01, position=position_dodge(0.9), bracket.nudge.y = 4, size=4, step.increase = 0.05) +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1)), breaks=seq(-10,20,5)) +
  theme_bw() + labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  guides(fill="none", pattern=guide_legend(title="Loop", override.aes=list(fill="white"))) +
  scale_pattern_manual(values = c("stripe","none"), labels=c("Unique","Shared")) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text = element_text(size=14, color="black", face="bold"),
        axis.title = element_text(size=14, face="bold"),
        legend.position = "top",
        legend.margin = margin(t=0,r=0,b=0,l=-10,unit='pt'),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13))
ggsave(file.path(projdir,"loop_domain_log2TPM.png"),dpi=300,width=3,height=4)


# All genes in loop domains -----------------------------------------------

all_aged_loop_genes = unique(c(diff_aged_loop_genes$gene, non_diff_aged_loop_genes$gene))
all_aged_non_loop_genes = unique(rownames(tpm.data)[!(rownames(tpm.data) %in% unique(c(diff_aged_loop_genes$gene, non_diff_aged_loop_genes$gene)))])
all_young_loop_genes = unique(c(diff_young_loop_genes$gene, non_diff_young_loop_genes$gene))
all_young_non_loop_genes = unique(rownames(tpm.data)[!(rownames(tpm.data) %in% unique(c(diff_young_loop_genes$gene, non_diff_young_loop_genes$gene)))])

all_loop_tpm = rbind(data.frame(avg.tpm.data[avg.tpm.data$gene %in% all_aged_loop_genes, ], Age="A", Loop="Loop"),
                     data.frame(avg.tpm.data[avg.tpm.data$gene %in% all_aged_non_loop_genes, ], Age="A", Loop="Non-Loop"),
                     data.frame(avg.tpm.data[avg.tpm.data$gene %in% all_young_loop_genes, ], Age="Y", Loop="Loop"),
                     data.frame(avg.tpm.data[avg.tpm.data$gene %in% all_young_non_loop_genes, ], Age="Y", Loop="Non-Loop"))

all_loop_tpm_test = all_loop_tpm %>%
  pivot_longer(c("A","Y"), names_to="TPM.Age", values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(Age = factor(Age, levels=c("Y","A")),
         log2tpm = log2(TPM)) %>%
  group_by(Age) %>%
  wilcox_test(log2tpm ~ Loop) %>%
  add_xy_position(x="Age", group="Loop", dodge=0.9)

all_loop_tpm_age_test = all_loop_tpm %>%
  pivot_longer(c("A","Y"), names_to="TPM.Age", values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(Age = factor(Age, levels=c("Y","A")),
         log2tpm = log2(TPM)) %>%
  group_by(Loop) %>%
  wilcox_test(log2tpm ~ Age) %>%
  add_xy_position(x="Age", group="Loop", dodge=0.9)

all_loop_tpm %>%
  pivot_longer(c("A","Y"), names_to="TPM.Age", values_to="TPM") %>%
  dplyr::filter(Age==TPM.Age) %>%
  mutate(Age = factor(Age, levels=c("Y","A"))) %>%
  ggplot(aes(x=Age, y=log2(TPM))) +
  geom_violin_pattern(aes(fill=Age, pattern=Loop), 
                    color="black",
                    pattern_fill = "black",
                    pattern_color = "black",
                    pattern_density = 0.1,
                    position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,Loop)), fill="white", color="black", outlier.shape=NA, width=0.2, position=position_dodge(0.9)) +
  stat_pvalue_manual(all_loop_tpm_test, tip.length=0.01, position=position_dodge(0.9), bracket.nudge.y = 1) +
  stat_pvalue_manual(all_loop_tpm_age_test, tip.length=0.01, position=position_dodge(0.9), bracket.nudge.y = 4, step.increase=0.05) +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1)), breaks=seq(-10,20,5)) +
  theme_bw() + labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  guides(fill="none", pattern=guide_legend(title="Loop", override.aes=list(fill="white"))) +
  scale_pattern_manual(values = c("stripe","none"), labels=c("Loop","Non-Loop")) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text = element_text(size=14, color="black", face="bold"),
        axis.title = element_text(size=14, face="bold"),
        legend.position = "top",
        legend.margin = margin(t=0,r=0,b=0,l=-40,unit='pt'),
        legend.title = element_text(size=13, face="bold"),
        legend.text = element_text(size=13))
ggsave(file.path(projdir,"loop_vs_nonloop_log2TPM.png"),dpi=300,width=3,height=4)
