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

tpm.data = read.table(file.path(projdir, "compartmentExpression", "RNA_transformed_tpm_minus_sva_contribs.txt"))
gene.sym = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
gene.sym.tbl = data.frame(table(gene.sym)) # find duplicate gene symbols
duplicate.genes = gene.sym.tbl[gene.sym.tbl$Freq>1,"gene.sym"]

tpm.data = tpm.data[-which(gene.sym %in% duplicate.genes), 1:8] # isolate day0 data for young and aged and remove gene duplicates
rownames(tpm.data) = unlist(lapply(rownames(tpm.data), function(x) unlist(strsplit(x, "_", fixed=T))[2]))
tpm.data = tpm.data[complete.cases(tpm.data), ] # remove NA rows
tpm.data = tpm.data[which(apply(tpm.data, 1, function(x) all(x>=0))), ] # keep only genes with TPM>=0 in either young or aged

DE.genes = unlist(read.table(file.path(projdir, "compartmentExpression", "DE_0.05padj_log2(1.5)LFC_genenames.csv"), sep="\t", as.is=T))
DE.genes.name = unlist(lapply(DE.genes, function(x) unlist(strsplit(x,"_",fixed=T))[2]), use.names=F)


# Make plots --------------------------------------------------------------

data.frame(tpm = as.matrix(unlist(tpm.data["Mpc1", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Mpc1_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Mpc2", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Mpc2_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Sesn2", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Sesn2_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Foxo3", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Foxo3_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Mef2a", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Mef2a_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Kdm5b", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Kdm5b_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Ddit3", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Ddit3_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Spry2", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Spry2_expression_bar.png"), dpi=300, width=2, height=1.5)

data.frame(tpm = as.matrix(unlist(tpm.data["Mxi1", c(1:3, 5:7)]), nrow=6, ncol=1), 
           Age = factor(rep(c("A","Y"), each=3), levels=c("A","Y"))) %>%
  ggplot(aes(x=tpm, y=Age)) + 
  stat_summary(geom="errorbar", width=0.5, fun.data="mean_se") + 
  stat_summary(aes(fill=Age), geom="bar", fun.data="mean_se") + 
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) + 
  theme_bw() + scale_fill_manual(values=pal_nejm()(2)) +
  labs(x="TPM", y=NULL) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg,
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold", color="black"),
        legend.position = "none")
ggsave(file.path(projdir, "Virtual4C", "Ddit3_expression_bar.png"), dpi=300, width=2, height=1.5)
