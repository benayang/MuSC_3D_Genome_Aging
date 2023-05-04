library(dplyr)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(tidyr)
library(ggsci)
library(rstatix)
library(plotgardener)

projdir <- "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\21_MDkNN"

# Read in necessary objects -----------------------------------------------

TAD_dir <- "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer"

aged.domain <- read.table(file.path(TAD_dir,"aged.merged","40kb","aged.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))
young.domain <- read.table(file.path(TAD_dir,"young.merged","40kb","young.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))
bookend.bound_df <- readRDS(file.path(TAD_dir,"bookend.bound_df.RDS"))
merged.classified.stats <- readRDS(file.path(TAD_dir,"hicInterIntraTAD","merged.classified.stats.RDS"))
tpm.data <- readRDS(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C",
                              "04_FANC", "tpm.data.RDS"))

# Get expressed genes in each TAD type ------------------------------------

library(GenomicRanges)
library(HelloRanges)

ab.switch <- read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression\\compartmentBed\\100kb",
                                 "ab.switch.bed")) %>% setNames(c("chr","start","end","group")) %>% GRanges()

# Get all gene promoters within aged and young TAD domains
tss_slop_1kb <- GRanges(read.table(file.path(TAD_dir,"TAD expression","tss_1kb_slop.bed"), col.names = c("chr","start","end","gene")))
get_young_promoter_code <- R_bedtools_intersect(a = tss_slop_1kb, 
                                               b = young.domain %>% dplyr::select(1:4) %>% setNames(c("chr","start","end","ID")) %>% GRanges(),
                                               wa = T, wb = T)
young_overlap_promoters <- eval(get_young_promoter_code)
get_aged_promoter_code <- R_bedtools_intersect(a = tss_slop_1kb, 
                                              b = aged.domain %>% dplyr::select(1:4) %>% setNames(c("chr","start","end","ID")) %>% GRanges(),
                                              wa = T, wb = T)
aged_overlap_promoters <- eval(get_aged_promoter_code)

merged.classified.stats <- merged.classified.stats %>% 
  group_by(Age) %>%
  mutate(bins_v2 = cut_interval(rank(domain_score)/n(), length=0.25))

table(merged.classified.stats$bins_v2, merged.classified.stats$Age)

ggplot(merged.classified.stats, aes(x=bins_v2, y=domain_score)) + geom_point(aes(color=Age))

# Add rearranged TAD type, domain score, and binned domain scores 
young_overlap_promoters@second$Type <- bookend.bound_df[match(young_overlap_promoters@second$ID, bookend.bound_df$ID.Young), ]$Type
young_overlap_promoters@second$domain_score <- merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$domain_score
young_overlap_promoters@second$bins <- merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$bins
young_overlap_promoters@second$bins_v2 <- merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$bins_v2

aged_overlap_promoters@second$Type <- bookend.bound_df[match(aged_overlap_promoters@second$ID, bookend.bound_df$ID.Aged), ]$Type
aged_overlap_promoters@second$domain_score <- merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$domain_score
aged_overlap_promoters@second$bins <- merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$bins
aged_overlap_promoters@second$bins_v2 <- merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$bins_v2

young_overlap_promoters_type <- cbind(gene = young_overlap_promoters@first$gene, 
                                     as.data.frame(young_overlap_promoters@second))
aged_overlap_promoters_type <- cbind(gene = aged_overlap_promoters@first$gene, 
                                    as.data.frame(aged_overlap_promoters@second))

# pheatmap(drop_na(tpm.data[unique(overlap_promoters$gene), c(1:3,5:7)]), scale="row", cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(200))
avg.tpm.data <- data.frame(A=rowMeans(tpm.data[,1:3]), Y=rowMeans(tpm.data[,5:7])) %>% tibble::rownames_to_column("gene")

overlap_promoters_df <- rbind(left_join(young_overlap_promoters_type, avg.tpm.data[,c("gene","Y")], by="gene") %>% 
                               dplyr::rename(TPM="Y") %>% mutate(Age="Young") %>% drop_na(),
                             left_join(aged_overlap_promoters_type, avg.tpm.data[,c("gene","A")], by="gene") %>% 
                               dplyr::rename(TPM="A") %>% mutate(Age="Aged") %>% drop_na()) %>%
  mutate(Type = factor(Type, levels=c("Shared","Shifted","Merged","Split","Combination")),
         Type_v2 = if_else(Type=="Shared","Shared","Unique"),
         Age = factor(Age, levels=c("Young","Aged")), 
         log2tpm = log2(TPM)) %>%
  group_by(Age) %>%
  mutate(tpm_percentile = rank(TPM)/length(TPM)) %>%
  ungroup()

# Co-regulation gene expression scores ------------------------------------

coregulation <- rbind(left_join(young_overlap_promoters_type, 
                               avg.tpm.data[,c("gene","Y","A")], by="gene") %>% mutate(Age = "Young"),
                     left_join(aged_overlap_promoters_type, 
                               avg.tpm.data[,c("gene","Y","A")], by="gene") %>% mutate(Age = "Aged")) %>%
  dplyr::filter(!is.na(Y) & !is.na(A)) %>%
  mutate(log2fc = log2(A/Y),
         Age = factor(Age, levels=c("Young", "Aged"))) %>%
  group_by(Age) %>%
  distinct(gene, .keep_all = T)
# coregulation %>% dplyr::filter(Type.Y != Type.A) %>% head()
# stopifnot(identical(coregulation$Type.Y, coregulation$Type.A))
coregulation %>%
  group_by(Age, ID) %>%
  summarise(coreg = abs(mean(sign(log2fc)))) %>%
  left_join(overlap_promoters_df %>% dplyr::select(ID, Type, Age), by=c("ID","Age")) %>%
  ggplot(aes(x=Type, y=coreg)) + 
  geom_violin(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age, Type)), color="black", width=0.2, position=position_dodge(0.9)) + 
  scale_fill_manual(values=pal_nejm()(2)[2:1]) +
  theme_pubr()

# Degree of disorder in TADs (DoD) ----------------------------------------

young_mdknn <- read.table(file.path(projdir, "young_individual", "MDkNN_young_10kb_individual.txt"), sep="\t", header=T)
aged_mdknn <- read.table(file.path(projdir, "aged_individual", "MDkNN_aged_10kb_individual.txt"), sep="\t", header=T)

merged_mdknn <- full_join(young_mdknn, aged_mdknn, by=c("ChromID","Start","End"), suffix=c(".Young",".Aged")) %>%
  left_join(merged.classified.stats %>% dplyr::select(chr,start,end,Type,Age,bins,domain_score), 
            by=c("ChromID"="chr","Start"="start","End"="end")) %>%
  pivot_longer(cols=c(DoD.Young, DoD.Aged), names_to="DoD.Age", values_to="DoD") %>%
  dplyr::filter((Age=="Y" & DoD.Age=="DoD.Young") | (Age=="A" & DoD.Age=="DoD.Aged")) %>%
  group_by(DoD.Age) %>%
  mutate(DoD_bins = cut_interval(rank(DoD)/n(), length=0.2),
         DoD_group = factor(ifelse(DoD > median(DoD, na.rm=T), "High", "Low"), levels=c("Low","High"))) %>%
  ungroup()
change_age <- c("A"="Aged", "Y"="Young")
merged_mdknn$Age2 <- factor(change_age[as.character(merged_mdknn$Age)], levels=c("Young","Aged"))

ggboxplot(merged_mdknn, x="Type", y="DoD", fill="Age")

# DoD histogram comparision
DoD_test <- merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  group_by(Age) %>%
  mutate(DoD_scaled = DoD/max(DoD, na.rm=T)) %>%
  ungroup() %>%
  wilcox_test(DoD_scaled ~ Age) 
DoD_test <- DoD_test %>% mutate(p.label = "<2.2e-16", y.position = 435, xmin=0.5, xmax=0.6, groups=list(V1=c("Y","A")))

merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  group_by(Age) %>%
  mutate(DoD_scaled = DoD/max(DoD, na.rm=T)) %>%
  arrange(DoD_scaled) %>%
  gghistogram(x="DoD_scaled", fill="Age2", bins=50) +
  stat_pvalue_manual(DoD_test, label = "p.label", label.size = 4) +
  scale_fill_manual(values=pal_nejm()(2)[2:1]) + 
  scale_color_manual(values=pal_nejm()(2)[2:1]) + 
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  labs(x = "Normalized Degree of Disorder", y = "# of TADs", fill = "Age") +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        axis.title = element_text(face="bold"),
        axis.text = element_text(face="bold")) 
ggsave(file.path(projdir, "scaled_DoD_hist.png"), dpi=300, width=4, height=3)

med_DoD <- median(merged_mdknn$DoD, na.rm=T)
table(merged_mdknn$DoD>med_DoD)

boxplot(list(young=young_mdknn$DoD, aged=aged_mdknn$DoD))
t.test(young_mdknn$DoD, aged_mdknn$DoD)

# merged_mdknn %>%
#   mutate(DoD_group = ifelse(DoD > median(DoD, na.rm=T), "High", "Low")) %>%
#   dplyr::filter(!is.na(DoD)) %>%
#   ggplot(aes(x=DoD, y=domain_score)) +
#   geom_point(aes(color=Age))

merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  ggplot(aes(x=DoD_bins, y=domain_score, color=Age)) +
  stat_summary(fun.data="mean_se", geom="pointrange")

merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  ggplot(aes(x=Age, y=domain_score)) +
  geom_violin_pattern(aes(fill=Age, pattern=DoD_group),
                       color="black",
                       pattern_density=0.1,
                       pattern_color="white",
                       pattern_fill="white",
                       position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age, DoD_group)), color="black", width=0.2, outlier.shape=NA, position=position_dodge(0.9)) +
  scale_fill_manual(values = pal_nejm()(2)[2:1]) +
  theme_bw() + labs(x = NULL, y = "Intra-TAD Connectivity") +
  scale_x_discrete(labels = c("Young", "Aged")) +
  guides(fill="none", pattern=guide_legend(title="Degree of Disorder", override.aes = list(fill="gray25"))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face='bold'),
        legend.text = element_text(size=12),
        legend.margin = margin(t=0, r=0, b=0, l=-40))
ggsave(file.path(projdir, "hi_lo_DoD_vs_domain_score.png"), dpi=300, width=4, height=4)

merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  ggplot(aes(x=DoD_bins, y=domain_score, color=Age)) +
  stat_summary(fun="mean", aes(group=Age, color=Age), geom="path", position=position_dodge(0.1)) +
  stat_summary(fun.data="mean_se", geom="pointrange") 
ggsave(file.path(projdir, "hi_lo_DoD_vs_domain_score_line.png"), dpi=300, width=4, height=4)

coregulation %>%
  dplyr::filter(Age=="Aged") %>%
  left_join(aged_mdknn, by=c("seqnames"="ChromID","start"="Start","end"="End")) %>%
  dplyr::select(-gene) %>%
  distinct() %>%
  mutate(log2Y = log2(Y)) %>%
  ggplot(aes(x=DoD, y=log2Y)) + 
  geom_hex(bins=50) +
  stat_smooth(method="lm", se = T) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  scale_fill_distiller(palette = "RdBu")
ggscatter(x="DoD", y="domain_score", add = "reg.line", conf.int = T, add.params = list(color="gray"))

# Gene expression in DoD TADs ---------------------------------------------

merged_mdknn %>%
  left_join(overlap_promoters_df %>% dplyr::select(gene, seqnames, start, end, Age, ID, log2tpm), 
            by=c("ChromID"="seqnames","Start"="start","End"="end","Age2"="Age")) %>%
  dplyr::filter(!is.na(DoD)) %>%
  group_by(ID, Age, DoD_group) %>%
  summarise(mean_log2tpm = mean(log2tpm, na.rm=T)) %>%
  ggplot(aes(x=Age, y=mean_log2tpm)) +
  geom_violin_pattern(aes(fill=Age, pattern=DoD_group),
                      color="black",
                      pattern_density=0.1,
                      pattern_color="white",
                      pattern_fill="white",
                      position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age, DoD_group)), color="black", width=0.2, outlier.shape=NA, position=position_dodge(0.9)) +
  scale_fill_manual(values = pal_nejm()(2)[2:1]) +
  theme_bw() + labs(x = NULL, y = expression(bold(paste("Mean lo","g"["2"],"(TPM) per TAD")))) +
  scale_x_discrete(labels = c("Young", "Aged")) +
  guides(fill="none", pattern=guide_legend(title="Degree of Disorder", override.aes = list(fill="gray25"))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face='bold'),
        legend.text = element_text(size=12),
        legend.margin = margin(t=0, r=0, b=0, l=-40))
ggsave(file.path(projdir, "hi_lo_DoD_log2tpm.png"), dpi=300, width=4, height=4)

merged_mdknn %>%
  dplyr::filter(!is.na(DoD)) %>%
  left_join(overlap_promoters_df %>% dplyr::select(gene, seqnames, start, end, Age, ID, log2tpm), 
            by=c("ChromID"="seqnames","Start"="start","End"="end","Age2"="Age")) %>%
  group_by(ID, Age, DoD_bins) %>%
  summarise(mean_log2tpm = mean(log2tpm)) %>%
  ggplot(aes(x=DoD_bins, y=mean_log2tpm, color=Age)) +
  stat_summary(fun="mean", aes(group=Age, color=Age), geom="path", position=position_dodge(0.1)) +
  stat_summary(fun.data="mean_se", geom="pointrange") +
  scale_color_manual(values = pal_nejm()(2)[2:1]) +
  theme_bw() + labs(x = "DoD Percentile Bins", y = expression(bold(paste("Mean lo","g"["2"],"(TPM) per TAD")))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face='bold'))
ggsave(file.path(projdir, "hi_lo_DoD_log2tpm_line.png"), dpi=300, width=5, height=4)

overlap_promoters_df %>%
  group_by(ID, Age, bins_v2) %>%
  summarise(mean_log2tpm = mean(log2tpm)) %>%
  ggplot(aes(x=bins_v2, y=mean_log2tpm, color=Age)) +
  stat_summary(fun.data="mean_se", geom="pointrange")


# Histone modifications per DoD group -------------------------------------

young_histone <- read.table(file.path(TAD_dir, "TAD_domain_track_coverage", "young_domain_scores_RPKM.tab"))
colnames(young_histone) <- c("chr","start","end","ATAC","H3K4me3")
aged_histone <- read.table(file.path(TAD_dir, "TAD_domain_track_coverage", "aged_domain_scores_RPKM.tab"))
colnames(aged_histone) <- c("chr","start","end","ATAC","H3K4me3")

merged_histone <- rbind(young_histone %>% mutate(Age = "Y"),
                       aged_histone %>% mutate(Age = "A"))

merged_mdknn %>%
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  ggscatter(x="DoD", y="ATAC", color="Age", add="reg.line", alpha=0.5, size = 1, conf.int = T, facet.by = "Age") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) 

merged_mdknn %>%
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  group_by(Age) %>%
  t_test(ATAC ~ bins, comparisons = list(c("[0,0.25]","(0.25,0.5]"),
                                         c("(0.25,0.5]","(0.5,0.75]"),
                                         c("(0.5,0.75]","(0.75,1]")))

merged_mdknn %>%
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  ggplot(aes(x=DoD_bins, y=ATAC, color=Age, group=Age)) + 
  stat_summary(fun.data="mean_se", geom='path') +
  stat_summary(fun.data="mean_se", geom='pointrange')
ggsave(file.path(projdir, "DoD_bin_ATAC_line.png"), dpi=300, width=4, height=3)

merged_mdknn %>%
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  ggplot(aes(x=DoD_bins, y=H3K4me3, color=Age, group=Age)) + 
  stat_summary(fun.data="mean_se", geom='path') +
  stat_summary(fun.data="mean_se", geom='pointrange')
ggsave(file.path(projdir, "DoD_bin_H3K4me3_line.png"), dpi=300, width=4, height=3)

merged_mdknn %>%
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  ggscatter(x="DoD", y="ATAC", color="Age", add = "reg.line", size=1)

merged_mdknn %>%  
  left_join(merged_histone, by=c("ChromID"="chr","Start"="start","End"="end","Age")) %>%
  dplyr::filter(!is.na(DoD)) %>%
  ggviolin(x="Age", y="ATAC", fill="DoD_group", add="boxplot")

# Histone peaks per DoD group -------------------------------------

library(GenomicRanges)

peaksDir <- "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression\\ATAC_overlap"
young_ATAC_peaks <- read.table(file.path(peaksDir, "atac.d0.Young.optimal.narrowPeak.gz"), sep='\t')
aged_ATAC_peaks <- read.table(file.path(peaksDir, "atac.d0.Aged.optimal.narrowPeak.gz"), sep='\t')

young_hits <- findOverlapPairs(GRanges(setNames(young.domain[,1:4], c("chr","start","end","ID"))), 
                               GRanges(setNames(young_ATAC_peaks[,1:4], c("chr","start","end","ID")))) 
young_hits_counts <- first(young_hits) %>% as.data.frame() %>% group_by(seqnames, start, end, ID) %>% tally(name = "ATAC")

aged_hits <- findOverlapPairs(GRanges(setNames(aged.domain[,1:4], c("chr","start","end","ID"))), 
                              GRanges(setNames(aged_ATAC_peaks[,1:4], c("chr","start","end","ID"))))
aged_hits_counts <- first(aged_hits) %>% as.data.frame() %>% group_by(seqnames, start, end, ID) %>% tally(name = "ATAC")

merged_histone_peaks <- rbind(young_hits_counts %>% mutate(Age = "Y"),
                              aged_hits_counts %>% mutate(Age = "A"))

merged_mdknn %>% 
  left_join(merged_histone_peaks, by=c("ChromID"="seqnames","Start"="start","End"="end","Age")) %>%
  filter(!is.na(DoD)) %>%
  ggplot(aes(x=DoD_bins,  y=ATAC, color=Age, group=Age)) +
  stat_summary(fun="mean", geom="line") +
  stat_summary(fun.data="mean_se", geom="pointrange")
  ggviolin(x="bins", y="ATAC", fill="Age")

# Generate bedpe for visualization of significant points ------------------

young_sig_mdknn <- read.table(file.path(projdir, "young_individual", "MDkNN_young_10kb_individual_longrange.txt"), 
                             sep='\t', header=T)
aged_sig_mdknn <- read.table(file.path(projdir, "aged_individual", "MDkNN_aged_10kb_individual_longrange.txt"), 
                             sep='\t', header=T)

young_sig_mdknn <- young_sig_mdknn %>%
  separate(tad_range, into=c("chr","start","end"), sep = "[:-]") %>%
  separate(bin1_range, into=paste("bin1", c("chr","start","end"), sep="_"), sep = "[:-]") %>%
  separate(bin2_range, into=paste("bin2", c("chr","start","end"), sep="_"), sep = "[:-]")
aged_sig_mdknn <- aged_sig_mdknn %>%
  separate(tad_range, into=c("chr","start","end"), sep = "[:-]", ) %>%
  separate(bin1_range, into=paste("bin1", c("chr","start","end"), sep="_"), sep = "[:-]") %>%
  separate(bin2_range, into=paste("bin2", c("chr","start","end"), sep="_"), sep = "[:-]")

young_sig_mdknn <- young_sig_mdknn %>% 
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  left_join(young_mdknn, by=c("chr"="ChromID", "start"="Start", "end"="End")) %>%
  mutate(bins = cut_interval(rank(DoD)/n(), length=0.25))
young_sig_mdknn <- left_join(young_sig_mdknn, 
                             merged.classified.stats %>% 
                               dplyr::filter(Age=="Y") %>% 
                               dplyr::select(chr,start,end,name,Type),
                             by = c("chr","start","end"))
young_sig_mdknn <- young_sig_mdknn %>% 
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  left_join(young_mdknn, by=c("chr"="ChromID", "start"="Start", "end"="End")) %>%
  mutate(bins = cut_interval(rank(DoD)/n(), length=0.25))

aged_sig_mdknn <- aged_sig_mdknn %>% 
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  left_join(young_mdknn, by=c("chr"="ChromID", "start"="Start", "end"="End")) %>%
  mutate(bins = cut_interval(rank(DoD)/n(), length=0.25))
aged_sig_mdknn <- left_join(aged_sig_mdknn, 
                            merged.classified.stats %>% 
                              dplyr::filter(Age=="A") %>% 
                              dplyr::select(chr,start,end,name,Type),
                            by = c("chr","start","end"))

saveRDS(young_sig_mdknn, file.path(projdir, "young_sig_mdknn.RDS"))
saveRDS(aged_sig_mdknn, file.path(projdir, "aged_sig_mdknn.RDS"))

aged_sig_mdknn %>% distinct(name, .keep_all=T) %>% slice_max(order_by=DoD, n=5)
aged_sig_mdknn %>% distinct(name, .keep_all=T) %>% slice_min(order_by=DoD, n=5)

young_sig_mdknn %>%
  dplyr::select(bin1_chr, bin1_start, bin1_end, bin2_chr, bin2_start, bin2_end, name, DoD) %>%
  setNames(c("chrom1","start1","end1","chrom2","start2","end2","name","score")) %>%
  write.table(file.path(projdir, "young_individual", "young_individual_longrange.bedpe"), sep="\t", col.names=T, row.names=F, quote=F)
aged_sig_mdknn %>%
  dplyr::select(bin1_chr, bin1_start, bin1_end, bin2_chr, bin2_start, bin2_end, name, DoD) %>%
  setNames(c("chrom1","start1","end1","chrom2","start2","end2","name","score")) %>%
  write.table(file.path(projdir, "aged_individual", "aged_individual_longrange.bedpe"), sep="\t", col.names=T, row.names=F, quote=F)


# Make DoD plots ----------------------------------------------------------

#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr2", chromstart = 32000000, chromend = 33000000)
#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr18", chromstart = 35000000, chromend = 36400000)
#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr1", chromstart = 191600000, chromend = 191840000)

#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr17", chromstart = 26240000, chromend = 26680000)
#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr3", chromstart = 88280000, chromend = 88560000)
#params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr17", chromstart = 17360000, chromend = 17720000)
params <- pgParams(assembly = "mm10", x = 0.25, width = 4.5, chrom="chr16", chromstart = 38840000, chromend = 39080000)

png(file.path(projdir, "lo_DoD_aged_plotgardener.png"), res=300, units="in", width=5, height=5)
pageCreate(width = 5, height = 5)
hic <- plotHicTriangle(data = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC\\aged.merged.hic",
                       params = params, 
                       matrix = "log2oe",
                       resolution = 10000,
                       bg = "black",
                       norm = "KR",
                       height = 3, y = 0,
                       zrange = c(0,20),
                       palette = colorRampPalette(colors = rev(RColorBrewer::brewer.pal("RdBu", n=11))))
# hic2 <- plotHicSquare(data = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC\\young.merged.hic",
#                      params = params, 
#                      matrix = "log2oe",
#                      resolution = 10000,
#                      norm = "KR",
#                      height = 4.5, y = 0,
#                      half = "bottom",
#                      zrange = c(0,20),
#                      palette = colorRampPalette(colors = rev(RColorBrewer::brewer.pal("RdBu", n=11))))

annoDomains(hic, 
            data = merged.classified.stats %>% ungroup() %>% 
              dplyr::filter(Age=="A") %>% dplyr::select(chr,start,end) %>% distinct(),
            params = params)
annoPixels(hic, 
           data = aged_sig_mdknn %>%
             dplyr::select(bin1_chr, bin1_start, bin1_end, bin2_chr, bin2_start, bin2_end, name, DoD) %>%
             setNames(c("chrom1","start1","end1","chrom2","start2","end2","name","score")) %>%
             mutate(start1 = as.numeric(start1),
                    end1 = as.numeric(end1),
                    start2 = as.numeric(start2),
                    end2 = as.numeric(end2)),
           type = "circle",
           shift = 0,
           col = "green",
           params = params)
# annoPixels(hic2, 
#            data = young_sig_mdknn %>%
#              dplyr::select(bin1_chr, bin1_start, bin1_end, bin2_chr, bin2_start, bin2_end, name, DoD) %>%
#              setNames(c("chrom1","start1","end1","chrom2","start2","end2","name","score")) %>%
#              mutate(start1 = as.numeric(start1),
#                     end1 = as.numeric(end1),
#                     start2 = as.numeric(start2),
#                     end2 = as.numeric(end2)),
#            type = "circle",
#            shift = 0,
#            col = "green",
#            params = params)
# TAD <- plotRanges(data = merged.classified.stats %>% ungroup() %>% 
#              dplyr::filter(Age=="A") %>% dplyr::select(chr,start,end,Type) %>% distinct(),
#            params = params, y = "b0", height = 0.5, 
#            fill = colorby("Type"))
annoGenomeLabel(plot=TAD, params = params, y = "b0")
pageGuideHide()
dev.off()