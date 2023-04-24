library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggsci)
library(ggpubr)
library(rstatix)
library(tidyr)
library(scales)
library(dplyr)
library(ggpattern)
library(latex2exp)

# Import data -------------------------------------------------------------

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS"

zscore.dirs = grep(".zip",list.files(file.path(projdir,"loop_GSEA"), pattern="loop_anchor_ORA_"), invert=T, value=T, fixed=T)
zscore.files = sapply(zscore.dirs, function(x) list.files(file.path(projdir,"loop_GSEA",x), pattern="enrichment_results_wg_result", full.names = T))
zscore.file.df = data.frame(group=names(zscore.files), fname=unlist(zscore.files,use.names=F))

zscore.list = list()
for(g in zscore.file.df$group) {
  tmp = read.delim2(zscore.file.df[zscore.file.df$group==g,"fname"], header=T, as.is=T)
  tmp = tmp %>% mutate(Group=g)
  zscore.list[[g]] = tmp
}
zscore.df = do.call(rbind,zscore.list)
zscore.df = zscore.df %>% 
  separate(Group, into=c("Age",NA,NA,NA,"geneset"), sep="_", remove=F) %>% 
  mutate(FDR = as.numeric(FDR), enrichmentRatio = as.numeric(enrichmentRatio))

zscore.df %>%
  group_by(Age) %>%
  arrange(enrichmentRatio) %>%
  mutate(description=factor(description,levels=description)) %>%
  filter(Age=="young") %>%
  ggplot(aes(x=enrichmentRatio,y=description)) +
  geom_col(aes(fill=FDR)) +
  facet_grid(rows=vars(Age), space="free", scales="free_y") +
  scale_fill_distiller(palette="RdBu")

zscore.df.filtered = zscore.df %>% 
  filter((Age=="aged" & (description %in% c("Notch signaling pathway","p53 signaling pathway", "TGF-beta signaling pathway",
                                            "FoxO signaling pathway", "Wnt signaling pathway", "HIF-1 signaling pathway", "Hippo signaling pathway", "PI3K-Akt signaling pathway"))) | 
           (Age=="young" & (description %in% c("FoxO signaling pathway", "cGMP-PKG signaling pathway", "stem cell differentiation",
                                               "Wnt signaling pathway", "Ras signaling pathway")))) %>%
  mutate(Age=factor(Age))

zscore.df.filtered %>%
  group_by(Age) %>%
  arrange(enrichmentRatio, .by_group=T) %>%
  mutate(description=factor(description,levels=unique(description))) %>%
  ungroup() %>%
  ggplot(aes(x=enrichmentRatio,y=description)) +
  geom_col(aes(fill=FDR)) +
  scale_x_continuous(expand=expansion(mult=c(0,0.1))) +
  facet_grid(rows=vars(Age), space="free", scales="free_y", labeller=labeller(Age=c(aged="Aged",young="Young"))) +
  scale_fill_distiller(palette="RdBu") +
  theme_bw() + labs(x="Enrichment Ratio", y=NULL) +
  theme(axis.text=element_text(size=12,color="black",face="bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12),
        strip.text=element_text(size=14,face="bold"))
ggsave(file.path(projdir,"loop_GSEA","loop_anchor_ORA.png"),dpi=300,width=5,height=6)
