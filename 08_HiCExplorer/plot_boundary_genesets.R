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

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer"

files = c(file.path(projdir,"WebGestalt","Lost_TAD_Boundary_expressed_promoters_ORA","enrichment_results_wg_result1662711318.txt"),
          file.path(projdir,"WebGestalt","Shared_TAD_Boundary_expressed_promoters_ORA","enrichment_results_wg_result1650836505.txt"))
names(files) = c("Lost","Shared")

file.list = list()
for(f in files) {
  tmp = read.delim2(f, header=T, as.is=T)
  tmp = tmp %>% mutate(Group=f)
  file.list[[f]] = tmp
}
file.df = do.call(rbind,file.list)


# Shared TAD boundaries ---------------------------------------------------

shared.filtered = file.df %>%
  dplyr::filter(Group == files[2]) %>% 
  dplyr::filter(unname(unlist(description)) %in% c("regulation of nuclease activity",
                                   "respiratory electron transport chain",
                                   "DNA damage checkpoint",
                                   "ERAD pathway",
                                   "tRNA metabolic process",
                                   "glycosyl compound metabolic process",
                                   "DNA biosynthetic process",
                                   "regulation of cell cycle phase transition",
                                   "transport along microtubule",
                                   "protein folding",
                                   "protein transport",
                                   "Post-translational protein modification",
                                   "Membrane Trafficking",
                                   "mRNA Splicing - Major Pathway",
                                   "rRNA processing in the nucleus and cytosol",
                                   "ribonucleoprotein complex assembly",
                                   "nucleotide metabolic process")) %>%
  mutate(FDR = as.numeric(FDR), enrichmentRatio = as.numeric(enrichmentRatio)) %>%
  arrange(enrichmentRatio) %>%
  mutate(description = factor(description, levels=description))

ggplot(shared.filtered, aes(y=description, x=enrichmentRatio)) +
  geom_col(aes(fill=FDR)) +
  scale_fill_distiller(palette="RdBu") +
  guides(fill = guide_colourbar(barwidth = 10, barheight = 1)) +
  theme_bw() + labs(x="Normalized Enrichment Score", y=NULL) +
  theme(axis.text=element_text(size=12,color="black",face="bold"),
        axis.title=element_text(size=12,face="bold"),
        plot.margin = margin(0,1,0,0,unit="pt"),
        legend.position="top",
        legend.margin = margin(t=0, r=0, b=0, l=-100, unit="pt"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12))
ggsave(file.path(projdir, "WebGestalt", "shared_TAD_boundary_geneset.png"), dpi=300, width=5, height=5)

# Lost TAD boundaries ---------------------------------------------------

lost.filtered = file.df %>%
  dplyr::filter(Group == files[1]) %>% 
  dplyr::filter(unname(unlist(description)) %in% c("DNA Replication Pre-Initiation",
                                                   "RNA Polymerase II Transcription",
                                                   "Orc1 removal from chromatin",
                                                   "Mitotic G1-G1/S",
                                                   "G1/S Transition",
                                                   "G2/M Checkpoints")) %>%
  mutate(FDR = as.numeric(FDR), enrichmentRatio = as.numeric(enrichmentRatio)) %>%
  arrange(enrichmentRatio) %>%
  mutate(description = factor(description, levels=description))

ggplot(lost.filtered, aes(y=description, x=enrichmentRatio)) +
  geom_col(aes(fill=FDR)) +
  scale_fill_distiller(palette="RdBu") +
  guides(fill = guide_colourbar(barwidth = 8)) +
  theme_bw() + labs(x="Normalized Enrichment Score", y=NULL) +
  theme(axis.text=element_text(size=12,color="black",face="bold"),
        axis.title=element_text(size=12,face="bold"),
        plot.margin = margin(0,1,0,0,unit="pt"),
        legend.position="top",
        legend.margin = margin(t=0, r=0, b=0, l=-100, unit="pt"),
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12))
ggsave(file.path(projdir, "WebGestalt", "lost_TAD_boundary_geneset.png"), dpi=300, width=5, height=3)
