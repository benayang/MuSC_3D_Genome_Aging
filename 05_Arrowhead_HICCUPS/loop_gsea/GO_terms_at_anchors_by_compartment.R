library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(tidyverse)

projdir <- "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C"


# Import gene expression, gene coordinates, and loop positions ---------------------------------------------------

gene_data <- genes(EnsDb.Mmusculus.v79)
seqlevelsStyle(gene_data) <- 'UCSC'
filt_gene_data <- gene_data[grep("pseudogene", gene_data$gene_biotype, invert=TRUE), ]

tpm.data <- readRDS(file.path(projdir, "04_FANC", "filtered_tpm_data.RDS"))

young_loops <- read.table(file.path(projdir,"05_Arrowhead_HICCUPS","young_merged_loops_noHeader.bedpe"),
                          col.names = c("chromosome1","x1","x2",
                                        "chromosome2","y1","y2",
                                        "blank1","blank2","blank3","blank4",
                                        "color","observed",
                                        "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                                      "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                                        "number_collapsed","centroid1","centroid2","radius"))
aged_loops <- read.table(file.path(projdir,"05_Arrowhead_HICCUPS","aged_merged_loops_noHeader.bedpe"),
                         col.names = c("chromosome1","x1","x2",
                                       "chromosome2","y1","y2",
                                       "blank1","blank2","blank3","blank4",
                                       "color","observed",
                                       "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                                       "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                                       "number_collapsed","centroid1","centroid2","radius"))

young_loops$Age <- "Y"
aged_loops$Age <- "A"

# stable = inner_join(aged_loops[,1:6], young_loops[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
# unique_aged = anti_join(aged_loops[,1:6], young_loops[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
# unique_young = anti_join(young_loops[,1:6], aged_loops[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
# 
# gained_loops <- unique_aged[,c("chromosome1","x1","y2")] %>%
#   setNames(c("seqnames","start","end")) %>%
#   GRanges()
# 
# lost_loops <- unique_young[,c("chromosome1","x1","y2")] %>%
#   setNames(c("seqnames","start","end")) %>%
#   GRanges()

common_loops <- inner_join(aged_loops[,1:6], young_loops[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
common_domains <- common_loops[,c("chromosome1","x1","y2")] %>%
  setNames(c("seqnames","start","end")) %>%
  GRanges()

young_domains <- young_loops[,c("chromosome1","x1","y2")] %>%
  setNames(c("seqnames","start","end")) %>%
  GRanges()

aged_domains <- aged_loops[,c("chromosome1","x1","y2")] %>%
  setNames(c("seqnames","start","end")) %>%
  GRanges()

domain_overlaps <- findOverlaps(aged_domains, young_domains)
gained_loops <- aged_domains[-from(domain_overlaps), ]
lost_loops <- young_domains[-to(domain_overlaps), ]

# young_loops <- GenomicInteractions(young_loops[,1:3] %>% setNames(c("seqnames", "start", "end")) %>% GRanges(), 
#                                    young_loops[,4:6] %>% setNames(c("seqnames", "start", "end")) %>% GRanges(), 
#                                    young_loops$observed)
# aged_loops <- GenomicInteractions(aged_loops[,1:3] %>% setNames(c("seqnames", "start", "end")) %>% GRanges(), 
#                                   aged_loops[,4:6] %>% setNames(c("seqnames", "start", "end")) %>% GRanges(), 
#                                   aged_loops$observed)

# Overlap loops with compartments -----------------------------------------

ab.switch = read.table(file.path(projdir, "04_FANC", "compartmentExpression", "compartmentBed", "100kb", "ab.switch.bed"), 
                       sep="\t", col.names=c('chrom','start','end','group'))
ab.switch = GRanges(ab.switch)

table(ab.switch$group)

get_gene_hits <- function(loops, ab_group) {
  hits <- findOverlaps(loops, ab.switch[ab.switch$group==ab_group,], type = "within")
  gene_hits <- findOverlaps(filt_gene_data, unique(loops[from(hits),]), select = "all", type = "any")
  gene_hit_data <- unique(drop_na(tpm.data[filt_gene_data$gene_name[from(gene_hits)], ]))
  return(gene_hit_data)
}

pheatmap(get_gene_hits(gained_loops, "staticA")[,c(1:3, 5:7)], scale = "row", show_rownames = FALSE)


# GSEA analysis -----------------------------------------------------------

write_s2n <- function(gene_hits) {
  s2n = data.frame(aged_avg = rowMeans(gene_hits[,1:3]),
                   young_avg = rowMeans(gene_hits[,5:7]),
                   aged_sd = matrixStats::rowSds(as.matrix(gene_hits[,1:3])),
                   young_sd = matrixStats::rowSds(as.matrix(gene_hits[,5:7])))
  s2n %>% 
    mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
    dplyr::select(s2n) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    #dplyr::filter(gene %in% DE.genes.name) %>%
    write.table(file="clipboard-1000", sep='\t', col.names=F, row.names=F, quote=F)
}

write_s2n(get_gene_hits(gained_loops, "staticA")) # some GO terms, no reactome terms
write_s2n(get_gene_hits(gained_loops, "staticB")) # no terms
write_s2n(get_gene_hits(lost_loops, "staticA")) # some GO terms, no reactome terms
write_s2n(get_gene_hits(lost_loops, "staticB")) # some GO terms, no reactome terms


# GO analysis -------------------------------------------------------------

writeClipboard(rownames(get_gene_hits(gained_loops, "staticA")))
writeClipboard(rownames(get_gene_hits(gained_loops, "staticB")))
writeClipboard(rownames(get_gene_hits(lost_loops, "staticA")))
writeClipboard(rownames(get_gene_hits(lost_loops, "staticB")))
writeClipboard(rownames(get_gene_hits(common_domains, "staticA")))
writeClipboard(rownames(get_gene_hits(common_domains, "staticB")))


# Heatmap of gene expression in terms -------------------------------------

hmp_genes <- c("Vegfa","Akt1","Itgb1","Cdkn1a","Trp53", "Rb1",
               "Foxo1","Foxo3","Foxo4","Foxo6",
               "Notch1","Mef2a","Mef2c",
               "Dab2ip","Arrb2","Ptk2",
               "Camk2b", "Camk2d", "Dlg1", "Dusp6",
               "Pik3ca","Pik3r1","Prkca",
               "Grb2","Fgfr1","Fgfr2","Fgfr3","Fgfr4",
               "Frs2","Gab1","Gab2",
               "Cdc42", "Rac1","Pip4k2c",
               "Mapk1","Mapk3","Mapk7","Mapk14",
               "Stat3")

hmp_genes <- c("Cnn1","F11r","Foxp1","Itgb1","Frmd7","Myh10","Mylk3","Mypn","Rock1")

filt_hmp_genes <- hmp_genes[hmp_genes %in% rownames(get_gene_hits(common_domains, "staticA"))]

pheatmap(tpm.data[hmp_genes, c(1:3, 5:7)],
         scale = "row", color = rev(colorRampPalette(colors = brewer.pal(11, "RdBu"))(200)))


# Plot of terms -----------------------------------------------------------

lost_terms <- read.delim2(file.path(projdir, "05_Arrowhead_HICCUPS", "loop_GSEA", "lost_loops_ORA_staticA", "enrichment_results_wg_result1670728994.txt"), 
                          sep = '\t', header=T, as.is=T)
stable_terms <- read.delim2(file.path(projdir, "05_Arrowhead_HICCUPS", "loop_GSEA", "stable_loops_ORA_staticA", "enrichment_results_wg_result1670731125.txt"), 
                            sep = '\t', header=T, as.is=T)
lost_terms = lost_terms %>% 
  mutate(FDR = as.numeric(FDR), enrichmentRatio = as.numeric(enrichmentRatio)) %>%
  dplyr::filter(description %in% c("VEGFA-VEGFR2 Pathway", "Signaling by VEGF", "actomyosin structure organization", "lamellipodium organization")) %>%
  mutate(group = "Lost")
stable_terms = stable_terms %>% 
  mutate(FDR = as.numeric(FDR), enrichmentRatio = as.numeric(enrichmentRatio)) %>%
  dplyr::filter(description %in% c("MAPK1/MAPK3 signaling", "RAF/MAP kinase cascade", "cell-cell signaling wnt", "Signaling by EGFR",
                                   "ERK1 and ERK2 cascade", "MAPK family signaling cascades", "PIP3 activates AKT signaling")) %>%
  mutate(group = "Stable")

terms <- rbind(lost_terms, stable_terms)

terms %>%
  arrange(enrichmentRatio) %>%
  mutate(description=factor(description,levels=description)) %>%
  ggplot(aes(x=enrichmentRatio, y=description)) +
  facet_grid(rows=vars(group), scales="free_y", space="free") +
  geom_col(aes(fill=FDR), color="black") +
  #scale_x_continuous(breaks=seq(-2,2,0.5)) +
  #guides(fill = guide_colourbar(barwidth = 10, barheight = 1, title.vjust=0.9)) +
  scale_fill_distiller(palette="RdBu") +
  theme_bw() + labs(x="Normalized Enrichment Score", y=NULL) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,face="bold"),
        plot.margin = margin(0,1,0,0,unit="pt"),
        legend.position="right",
        legend.title=element_text(size=12,face="bold"),
        legend.text=element_text(size=12),
        strip.text=element_text(size=14,face="bold"))
ggsave(file.path(projdir, "05_Arrowhead_HICCUPS", "loop_GSEA", "loop_domain_staticA_ORA.png"),dpi=300,width=8,height=4)

# Plot Vegfa expression ---------------------------------------------------

tpm.data["Vegfa", c(1:3,5:7)] %>%
  rownames_to_column("gene") %>%
  pivot_longer(!gene) %>%
  mutate(Age = sapply(name, function(x) unlist(strsplit(x, split="_", fixed=TRUE))[2])) %>%
  mutate(Age = sapply(Age, function(x) ifelse(x=="A", "Aged", "Young"))) %>%
  mutate(Age = factor(Age, levels=c("Young", "Aged"))) %>%
  ggplot(aes(x=Age, y=value)) +
  stat_summary(aes(color=Age), fun.data="mean_se", geom="crossbar", width=0.5) +
  scale_color_manual(values = pal_nejm()(2)[2:1]) +
  geom_point() +
  theme_bw() + 
  ylab("Vegfa (TPM)") +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5, size=5) +
  theme(legend.position = "none",
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        axis.title.x = element_blank())
ggsave(file.path(projdir, "05_Arrowhead_HICCUPS", "Vegfa_expression.png"), dpi=300, width=3, height=4)
