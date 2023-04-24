library(ggplot2)
library(dplyr)
library(tidyr)
library(HelloRanges)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS"


# Loop domain data --------------------------------------------------------

young=read.table(file.path(projdir,"loop_domain_score","young.merged_loop_domain_observed.csv"), sep=",", header=T)
aged=read.table(file.path(projdir,"loop_domain_score","aged.merged_loop_domain_observed.csv"), sep=",", header=T)

young_bed = GRanges(young)
aged_bed = GRanges(aged)

merged_loops = bind_rows(young, aged) %>%
  mutate(Age = rep(c("Young","Aged"), times=c(nrow(young), nrow(aged))), 
         inter_norm_sum = interRegion_sum/interRegion_count_nnz,
         intra_norm_sum = intraRegion_sum/intraRegion_count_nnz,
         domain_score = intra_norm_sum / (inter_norm_sum + intra_norm_sum)) %>%
  group_by(Age) %>%
  mutate(bins = cut_width(rank(domain_score)/n(), width=0.2, center=0.5)) %>%
  ungroup()

ggplot(merged_loops, aes(x=bins, y=domain_score)) +
  geom_boxplot(aes(fill=Age))


# Loop domain score by compartment ----------------------------------------

ab.switch = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression\\compartmentBed\\100kb",
                                 "ab.switch.bed"), sep="\t", col.names=c('chrom','start','end','group'))
ab.switch = GRanges(ab.switch)

young_loop_list = list()
aged_loop_list = list()
for(g in unique(ab.switch$group)) {
  young_loop_list[[g]] = eval(R_bedtools_intersect(young_bed, ab.switch[ab.switch$group==g, ], wa=T)) %>% 
    unique() %>% 
    as.data.frame() %>%
    mutate(group = g, 
           Age = "Young",
           inter_norm_sum = interRegion_sum / interRegion_count_nnz,
           intra_norm_sum = intraRegion_sum / intraRegion_count_nnz,
           domain_score = intra_norm_sum / (inter_norm_sum + intra_norm_sum),
           bins = cut_width(rank(domain_score)/n(), width=0.2, center=0.5))

  aged_loop_list[[g]] = eval(R_bedtools_intersect(aged_bed, ab.switch[ab.switch$group==g, ], wa=T)) %>% 
    unique() %>%
    as.data.frame() %>%
    mutate(group = g,
           Age = "Aged",
           inter_norm_sum = interRegion_sum / interRegion_count_nnz,
           intra_norm_sum = intraRegion_sum / intraRegion_count_nnz,
           domain_score = intra_norm_sum / (inter_norm_sum + intra_norm_sum),
           bins = cut_width(rank(domain_score)/n(), width=0.2, center=0.5))
}

young_loop_df = bind_rows(young_loop_list)
aged_loop_df = bind_rows(aged_loop_list)
loop_df = bind_rows(young_loop_df, aged_loop_df)
loop_bed = GRanges(loop_df)

ggplot(loop_df, aes(x=group, y=domain_score)) +
  geom_boxplot(aes(fill=Age))

plot(table(loop_df$Age, loop_df$group))

# Overlap with genes and expression ---------------------------------------

tss_slop_1kb = readRDS(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\14_Multiome",
                                 "tss_slop_1kb.RDS"))

loop_bed_genes = eval(R_bedtools_intersect(tss_slop_1kb, loop_bed, wa=T, wb=T))
tpm.data = readRDS(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC", 
                             "tpm.data.RDS"))
avg.tpm = data.frame(A = rowMeans(tpm.data[,1:3]),
                     Y = rowMeans(tpm.data[,5:7]))

loop_expression = cbind(as.data.frame(second(loop_bed_genes)), 
                        as.data.frame(mcols(first(loop_bed_genes))))
loop_expression = cbind(loop_expression, avg.tpm[loop_expression$gene, ])
loop_expression = loop_expression %>% drop_na() %>% distinct()

loop_expression %>%
  #mutate(log2fc = log2(A/Y)) %>%
  pivot_longer(cols=c(A,Y), names_to="TPM_Age", values_to="TPM") %>%
  filter((TPM_Age=="A" & Age=="Aged") | (TPM_Age=="Y" & Age=="Young")) %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  ggplot(aes(x=group, y=log2(TPM))) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size=0.2) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir,"Figures","loop_domain_compartment_expression.png"),dpi=300,width=5,height=3.5)

loop_expression %>%
  pivot_longer(cols=c(A,Y), names_to="TPM_Age", values_to="TPM") %>%
  filter((TPM_Age=="A" & Age=="Aged") | (TPM_Age=="Y" & Age=="Young")) %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  group_by(TPM_Age) %>%
  mutate(TPM_percentile = rank(TPM)/length(TPM)) %>%
  ggplot(aes(x=bins, y=TPM_percentile)) +
  facet_wrap(~Age) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size=0.2) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  labs(x="Loop Domain Score Percentile Bins", y="Gene Expression (TPM) Percentile") +
  theme(axis.text.x = element_text(size=12, color="black", face="bold", angle=35, hjust=1),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        strip.text = element_text(size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "none")
ggsave(file.path(projdir,"Figures","loop_domain_type_expression.png"),dpi=300,width=4.25,height=3.5)

loop_expression %>%
  pivot_longer(cols=c(A,Y), names_to="TPM_Age", values_to="TPM") %>%
  filter((TPM_Age=="A" & Age=="Aged") | (TPM_Age=="Y" & Age=="Young")) %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  ggplot(aes(x=bins, y=log2(TPM))) +
  stat_summary(fun.data="mean_se", geom="path", aes(color=Age, group=Age)) +
  stat_summary(fun.data="mean_se", geom="pointrange", aes(color=Age))

# Classify loops by TADs --------------------------------------------------

merged.classified.stats = readRDS(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                            "hicInterIntraTAD","merged.classified.stats.RDS"))

young_TAD_bed = merged.classified.stats %>% ungroup() %>% filter(Age=="Y") %>% 
  dplyr::select(chr, start, end, name, score, domain_score, bins, Type) %>% GRanges()
aged_TAD_bed = merged.classified.stats %>% ungroup() %>% filter(Age=="A") %>% 
  dplyr::select(chr, start, end, name, score, domain_score, bins, Type) %>% GRanges()

young_loop_anchors = read.table(file.path(projdir,"young_merged_loops_noHeader.bedpe"),
                                col.names = c("chromosome1","x1","x2",
                                              "chromosome2","y1","y2",
                                              "blank1","blank2","blank3","blank4",
                                              "color","observed",
                                              "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                                              "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                                              "number_collapsed","centroid1","centroid2","radius"))
aged_loop_anchors = read.table(file.path(projdir,"aged_merged_loops_noHeader.bedpe"),
                               col.names = c("chromosome1","x1","x2",
                                             "chromosome2","y1","y2",
                                             "blank1","blank2","blank3","blank4",
                                             "color","observed",
                                             "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                                             "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                                             "number_collapsed","centroid1","centroid2","radius"))

young_loop_anchors = inner_join(merged_loops %>% dplyr::filter(Age=="Young"),
                                young_loop_anchors,
                                by=c("chr"="chromosome1","start"="x1", "end"="y2")) %>%
  mutate(loop_ID = paste0("loop_",seq_len(n())))
aged_loop_anchors = inner_join(merged_loops %>% dplyr::filter(Age=="Aged"),
                                aged_loop_anchors,
                                by=c("chr"="chromosome1","start"="x1", "end"="y2")) %>%
  mutate(loop_ID = paste0("loop_",seq_len(n())))

head(young_loop_anchors)
left_young_anchor = young_loop_anchors %>% dplyr::select(chr, start, x2, loop_ID) %>% setNames(c("chr","start","end","loop_ID")) %>% GRanges()
right_young_anchor = young_loop_anchors %>% dplyr::select(chr, y1, end, loop_ID) %>% setNames(c("chr","start","end","loop_ID")) %>% GRanges()
left_aged_anchor = aged_loop_anchors %>% dplyr::select(chr, start, x2, loop_ID) %>% setNames(c("chr","start","end","loop_ID")) %>% GRanges()
right_aged_anchor = aged_loop_anchors %>% dplyr::select(chr, y1, end, loop_ID) %>% setNames(c("chr","start","end","loop_ID")) %>% GRanges()

loop_TAD_olap = function(gr_a, gr_b) {
  # bedtools_intersect(wa=T, wb=T)
  ans = eval({
    pairs = findOverlapPairs(gr_a, gr_b, ignore.strand=T)
    olap = pintersect(pairs, ignore.strand = T)
    mcols(pairs)$overlap_width = width(olap)
    pairs
  })
  
  # convert to data.frame
  ans_df = cbind(as.data.frame(first(ans)) %>% dplyr::select(seqnames, start, end, loop_ID),
                 as.data.frame(second(ans)) %>% dplyr::select(start, end, name, width, score, domain_score, bins, Type)) %>%
    setNames(c("seqnames","start","end","loop_ID","TAD_start","TAD_end","TAD_name","TAD_width","TAD_score","TAD_domain_score","TAD_bins","TAD_Type")) %>%
    mutate(overlap_width = unlist(mcols(ans)))
  
  # Some anchors only overlap TADs by 1 bp. Get unique loop anchors by maximum overlap width
  ans_df = ans_df %>% group_by(loop_ID) %>% dplyr::slice(which.max(overlap_width)) %>% arrange(seqnames, start, loop_ID) %>% ungroup() %>% distinct()
    
  return(ans_df)
}

merge_loop_TAD_olaps = function(gr_a, gr_b) {
  # bring left and right loop anchors back into the same data frame
  # classify loops as asymmetric-TAD if one anchor is in a TAD and the other is outside of any TAD
  merged_olap = full_join(gr_a, gr_b, by=c("seqnames","loop_ID"), suffix=c("_left","_right")) %>%
    mutate(loop_width = end_right - start_left, 
           Type = if_else(!is.na(TAD_name_left) & !is.na(TAD_name_right), 
                          if_else(TAD_name_left==TAD_name_right, 
                                  "Intra-TAD", 
                                  "Inter-TAD"),
                          "Asymmetric-TAD"))
  return(merged_olap)
}

classify_loop_TAD_olaps = function(loop_anchors, merged_loop_TAD_olaps) {
  # if none of the loop anchors overlap a TAD, the loop is outside of a TAD (called Non-TAD) 
  out = left_join(loop_anchors, merged_loop_TAD_olaps, by="loop_ID")
  out$Type[is.na(out$Type)] = "Non-TAD" 
  
  # Classify loops that lie entirely within a single TAD
  intraTAD_idx = which(out$Type == "Intra-TAD")
  # double check that we're actually looking within the same TAD
  stopifnot(all(out$TAD_name_left[intraTAD_idx] == out$TAD_name_right[intraTAD_idx]))
  # If *both* loop anchors are within 20kb of the edge of a TAD domain, it's a loop domain
  loop_domain_IDs = out %>% dplyr::slice(intraTAD_idx) %>% filter(abs(start-TAD_start_left)<=2e4 & abs(end-TAD_end_right)<=2e4) %>% pull(loop_ID)
  # If *only one* loop anchor is within 20kb of the edge of a TAD domain, it's a TAD boundary loop
  boundary_loop_IDs = out %>% dplyr::slice(intraTAD_idx) %>% filter(xor(abs(start-TAD_start_left)<=2e4, abs(end-TAD_end_right)<=2e4)) %>% pull(loop_ID)
  # Make sure loops are uniquely classified
  stopifnot(!any(loop_domain_IDs %in% boundary_loop_IDs))
  out$Type[out$loop_ID %in% loop_domain_IDs] = "Loop Domain"
  out$Type[out$loop_ID %in% boundary_loop_IDs] = "TAD Boundary"
  
  return(out)
}

left_young_anchor_TAD_olap = loop_TAD_olap(left_young_anchor, young_TAD_bed)
right_young_anchor_TAD_olap = loop_TAD_olap(right_young_anchor, young_TAD_bed)
left_aged_anchor_TAD_olap = loop_TAD_olap(left_aged_anchor, aged_TAD_bed)
right_aged_anchor_TAD_olap = loop_TAD_olap(right_aged_anchor, aged_TAD_bed)

merged_young_anchor_TAD_olap = merge_loop_TAD_olaps(left_young_anchor_TAD_olap, right_young_anchor_TAD_olap)
merged_aged_anchor_TAD_olap = merge_loop_TAD_olaps(left_aged_anchor_TAD_olap, right_aged_anchor_TAD_olap)

young_loop_anchors_merged = classify_loop_TAD_olaps(young_loop_anchors, merged_young_anchor_TAD_olap)
aged_loop_anchors_merged = classify_loop_TAD_olaps(aged_loop_anchors, merged_aged_anchor_TAD_olap)

table(young_loop_anchors_merged$Type)
table(young_loop_anchors_merged$TAD_Type_right)
table(aged_loop_anchors_merged$Type)

classified_loop_df = bind_rows(young_loop_anchors_merged %>% mutate(loop_ID = paste("Young",loop_ID,sep="_")), 
                               aged_loop_anchors_merged %>% mutate(loop_ID = paste("Aged",loop_ID,sep="_")))
table(classified_loop_df$Age, classified_loop_df$TAD_Type_left)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plt = classified_loop_df %>%
  group_by(Age, Type) %>%
  summarise(num=n()) %>%
  ungroup() %>%
  group_by(Age) %>%
  mutate(fraction = num/sum(num)) %>%
  arrange(num) %>%
  mutate(Type = factor(Type, levels=unique(Type)),
         Age = factor(Age, levels=c("Young","Aged"))) %>%
  ungroup() %>%
  ggplot(aes(x=Age, y=num)) +
  geom_col(aes(fill=Type), color="black", position=position_stack()) +
  # geom_text_repel(aes(group=Age, label = sprintf("%d (%s)", num, scales::percent(fraction, accuracy=0.1))),
  #                  size=4, position = position_dodge(0.9), direction="x") +
  # geom_text_repel(aes(group=Age, label = sprintf("%s", scales::percent(fraction, accuracy=0.1))),
  #                 size=4, position = position_dodge(0.9), hjust=0.9, direction="x") +
  #scale_x_continuous(expand=expansion(mult=c(0,0.2)), breaks=seq(0,2500,500)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), breaks=seq(0,4000,500)) +
  #scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_fill_manual(values=rev(cbPalette)) +
  guides(fill=guide_legend(nrow=3)) +
  labs(x=NULL, y="# of Loops") +
  #labs(x="# of Loops", y=NULL) +
  theme_bw() +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.ticks.length.x = unit(0,"in"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "right",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"Figures","loops_TAD_overlap_classification.png"),dpi=300,width=5,height=4)

plt_legend = cowplot::get_legend(plt)
png(file.path(projdir,"Figures","loops_TAD_overlap_classification_legend.png"),res=300,units="in",width=5,height=5)
grid::grid.newpage()
grid::grid.draw(plt_legend)
# print(plt_legend)
dev.off()

## loop anchor annotation by TAD type ------------------

sum(table(young_loop_anchors_merged$TAD_Type_left, young_loop_anchors_merged$TAD_Type_right))
sum(table(aged_loop_anchors_merged$TAD_Type_left, aged_loop_anchors_merged$TAD_Type_right))

table(young_loop_anchors_merged$TAD_Type_left, young_loop_anchors_merged$TAD_Type_right) %>% reshape2::melt() %>%
  ggplot(aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value)) +
  geom_text(aes(label=value)) +
  labs(x="Left Anchor", y="Right Anchor") +
  scale_fill_distiller(palette="RdBu") + 
  theme_minimal() + coord_equal() +
  ggtitle("Young") +
  theme(legend.position="none", 
        title = element_text(size=12, face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        axis.text.x=element_text(size=12,color="black",angle=35,hjust=1),
        axis.text.y=element_text(size=12,color="black"))
ggsave(file.path(projdir,"Figures","young_loop_anchor_TAD_type_hmp.png"), dpi=300, width=4, height=3)

table(aged_loop_anchors_merged$TAD_Type_left, aged_loop_anchors_merged$TAD_Type_right) %>% reshape2::melt() %>%
  ggplot(aes(x=Var1, y=Var2)) +
  geom_tile(aes(fill=value)) +
  geom_text(aes(label=value)) +
  labs(x="Left Anchor", y="Right Anchor") +
  scale_fill_distiller(palette="RdBu") + 
  theme_minimal() + coord_equal() +
  ggtitle("Aged") +
  theme(legend.position="none", 
        title = element_text(size=12, face="bold"),
        axis.title=element_text(size=12,face="bold"), 
        axis.text.x=element_text(size=12,color="black",angle=35,hjust=1),
        axis.text.y=element_text(size=12,color="black"))
ggsave(file.path(projdir,"Figures","aged_loop_anchor_TAD_type_hmp.png"), dpi=300, width=4, height=3)

table(young_loop_anchors_merged$TAD_Type_left, young_loop_anchors_merged$TAD_Type_right, young_loop_anchors_merged$Type)

table(young_loop_anchors_merged$TAD_Type_left, young_loop_anchors_merged$TAD_Type_right) %>% as.matrix() %>% pheatmap()
table(aged_loop_anchors_merged$TAD_Type_left, aged_loop_anchors_merged$TAD_Type_right) %>% as.matrix() %>% pheatmap()

# Loop and TAD domain scores by type --------------------------------------

classified_loop_df %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  ggplot(aes(y=Type, x=domain_score))+
  geom_violin(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age, Type)), width=0.2, position=position_dodge(0.9), fill="white", color="black", outlier.shape=NA) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() +
  labs(x="Loop Domain Score", y=NULL) +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        strip.text = element_text(size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir,"Figures","loop_domain_score_by_type.png"),dpi=300,width=4.25,height=4.5)

intra_TAD_domain_test = classified_loop_df %>%
  filter(Type %in% c("Intra-TAD","TAD Boundary","Loop Domain")) %>%
  group_by(Age) %>%
  wilcox_test(TAD_domain_score_left ~ Type)
classified_loop_df %>%
  filter(Type %in% c("Intra-TAD","TAD Boundary","Loop Domain")) %>%
  ggplot(aes(x=Age, y=TAD_domain_score_left))+
  geom_boxplot(aes(fill=Type))

classified_loop_df %>%
  filter(Type == "Inter-TAD") %>%
  pivot_longer(cols=c(TAD_domain_score_left, TAD_domain_score_right), names_to="TAD_domain_score_direction", values_to="TAD_domain_score") %>%
  ggplot(aes(x=Age, y=TAD_domain_score)) +
  geom_boxplot(aes(fill=TAD_domain_score_direction))


# Annotate loops by their anchors ---------------------------------------------------

# Get loop anchor annotations from ChipSeeker
young.peakAnno = readRDS(file.path(projdir,"young_peakAnno.RDS"))
aged.peakAnno = readRDS(file.path(projdir,"aged_peakAnno.RDS"))

# ChipSeeker adds 1 to the start of each region?
young.peakAnno_df = as.data.frame(young.peakAnno@anno) %>% mutate(Age = "Young", start=start-1)
aged.peakAnno_df = as.data.frame(aged.peakAnno@anno) %>% mutate(Age = "Aged", start=start-1)
peakAnno_df = rbind(young.peakAnno_df, aged.peakAnno_df)

# Rename remove gene-specific information for intron and exon annotations
classified_loop_df$left_annotation = left_join(classified_loop_df, peakAnno_df, by=c("chr"="seqnames","start","x2"="end","Age")) %>% pull(annotation)
classified_loop_df$left_annotation = sapply(classified_loop_df$left_annotation, function(x) ifelse(grepl("Exon",x),"Exon",
                                                                                                   ifelse(grepl("Intron",x),"Intron",x)))
classified_loop_df$right_annotation = left_join(classified_loop_df, peakAnno_df, by=c("chr"="seqnames","y1"="start","end","Age")) %>% pull(annotation)
classified_loop_df$right_annotation = sapply(classified_loop_df$right_annotation, function(x) ifelse(grepl("Exon",x),"Exon",
                                                                                                     ifelse(grepl("Intron",x),"Intron",x)))

table(classified_loop_df$left_annotation, classified_loop_df$Age)
table(classified_loop_df$right_annotation, classified_loop_df$Age)

classified_loop_df$anno_type = "Other"
distal_anno_names = c("Distal Intergenic", "3' UTR", "5' UTR", "Downstream (<=300bp)")
classified_loop_df$anno_type[classified_loop_df$left_annotation=="Promoter" & classified_loop_df$right_annotation=="Promoter"] = "Promoter-Promoter"
classified_loop_df$anno_type[(classified_loop_df$left_annotation %in% distal_anno_names) & (classified_loop_df$right_annotation %in% distal_anno_names)] = "Distal-Distal"
classified_loop_df$anno_type[((classified_loop_df$left_annotation=="Promoter") & (classified_loop_df$right_annotation %in% distal_anno_names)) | 
                               ((classified_loop_df$left_annotation %in% distal_anno_names) & (classified_loop_df$right_annotation=="Promoter"))] = "Promoter-Distal"

classified_loop_df = classified_loop_df %>%
  mutate(left_anno_type = ifelse(grepl("Promoter", classified_loop_df$left_annotation), 
                                 "Promoter","Distal"),
         right_anno_type = ifelse(grepl("Promoter", classified_loop_df$right_annotation), 
                                  "Promoter","Distal"),
         anno_type = ifelse(left_anno_type == right_anno_type, 
                            paste(left_anno_type, right_anno_type, sep="-"),
                            paste("Promoter", "Distal", sep="-")))

# classified_loop_df$anno_type = ifelse(grepl("Promoter", classified_loop_df$left_annotation) & grepl("Promoter", classified_loop_df$right_annotation), 
#                                       "Promoter-Promoter", 
#                                       ifelse(xor(grepl("Promoter", classified_loop_df$left_annotation), grepl("Promoter", classified_loop_df$right_annotation)),
#                                              "Promoter-Distal",
#                                              "Distal-Distal"))
table(classified_loop_df$anno_type, classified_loop_df$Age, classified_loop_df$Type)

saveRDS(classified_loop_df, file.path(projdir, "classified_loop_df.RDS"))
classified_loop_df = readRDS(file.path(projdir, "classified_loop_df.RDS"))

# barplots of classifications
young_anno_tbl = with(classified_loop_df %>% dplyr::filter(Age=="Young" & Type != "Asymmetric-TAD"), table(anno_type, Type))
aged_anno_tbl = with(classified_loop_df %>% dplyr::filter(Age=="Aged" & Type != "Asymmetric-TAD"), table(anno_type, Type))
anno_tbl = with(classified_loop_df %>% dplyr::filter(Type != "Asymmetric-TAD"), table(anno_type, Type, Age))
plot(young_anno_tbl)
par(mar = c(5, 5, 6, 5))
barplot(prop.table(young_anno_tbl, margin=2), 
        col=c("#8A4198FF","#5A9599FF","#FF6348FF"),
        legend.text = TRUE, 
        args.legend = list(x = "top",
                           xpd = T,
                           ncol = 1,
                           bty = "n",
                           inset = c(0, -0.3)))
dev.off()

as.data.frame(anno_tbl) %>%
  group_by(Age, Type) %>%
  mutate(Fraction = Freq/sum(Freq)) %>%
  group_by(Type) %>%
  arrange(desc(Freq)) %>%
  mutate(Type = factor(Type, levels=c("Inter-TAD","Intra-TAD","Loop Domain","TAD Boundary","Non-TAD"))) %>%
  ggplot(aes(x=Type, y=Fraction, fill=anno_type)) +
  facet_wrap(~Age) +
  geom_col(color="black") + 
  geom_text(aes(label=scales::percent(Fraction,accuracy=1)), position=position_stack(vjust=0.5)) +
  scale_fill_npg() + guides(fill = guide_legend(nrow=2)) +
  theme_pubclean() + labs(x=NULL, y="Percentage of Loops") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.text.x = element_text(size=12, face="bold", color="black", angle=35, hjust=1),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_blank(),
        legend.text = element_text(size=12))
ggsave(file.path(projdir, "Figures", "loop_annotation_stacked_barchart.png"), dpi=300, width=5, height=4.5)

fisher.test(as.matrix(young_anno_tbl), hybrid=F, simulate.p.value=T)
f.posthoc = list(fisher.test(young_anno_tbl[c(1,2), ], hybrid=F, simulate.p.value=T),
                 fisher.test(young_anno_tbl[c(2,3), ], hybrid=F, simulate.p.value=T),
                 fisher.test(young_anno_tbl[c(1,3), ], hybrid=F, simulate.p.value=T))
f.posthoc.pvalues = lapply(f.posthoc, function(x) x$p.value)
names(f.posthoc.pvalues) = c(paste(rownames(young_anno_tbl)[1], rownames(young_anno_tbl)[2], sep="."),
                             paste(rownames(young_anno_tbl)[1], rownames(young_anno_tbl)[3], sep="."),
                             paste(rownames(young_anno_tbl)[2], rownames(young_anno_tbl)[3], sep="."))
format(p.adjust(f.posthoc.pvalues, method="holm"), scientific = T)

#mosaicplot(prop.table(young_anno_tbl, margin=2), shade=T)

# alluvial plots of classifications
library(ggalluvial)

classified_loop_df %>%
  filter(Type != "Asymmetric-TAD") %>%
  group_by(Age, anno_type, Type) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>%
  mutate(Type = factor(Type, levels=unique(Type))) %>%
  as.data.frame() %>%
  ggplot(aes(y=count, axis1=Age, axis2=anno_type, axis3=Type, fill=Age)) +
  geom_alluvium(width=1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))


# Get genes within classified domains -------------------------------------

domain_genes = eval(R_bedtools_intersect(tss_slop_1kb, GRanges(classified_loop_df[,c("chr","start","end","Age","loop_ID","anno_type")]), wa=T, wb=T))
domain_genes_df = cbind(second(domain_genes) %>% as.data.frame() %>% dplyr::select(seqnames,start,end,Age,loop_ID,anno_type),
                        first(domain_genes) %>% as.data.frame() %>% dplyr::select(gene)) %>%
  distinct()
domain_genes_df = cbind(domain_genes_df, avg.tpm[domain_genes_df$gene, ]) %>%
  drop_na() %>%
  mutate(log2fc = log2(A/Y))

ggviolin(domain_genes_df, x="anno_type", y="log2fc", fill="Age", add="boxplot") + geom_hline(yintercept=0)

domain_genes_df %>% filter(Age=="Aged" & anno_type=="Distal-Distal") %>%
  pull(gene) %>% unique() %>% writeClipboard()

domain_genes_df %>% filter(Age=="Young" & anno_type=="Distal-Distal") %>%
  dplyr::select(gene, log2fc) %>%
  distinct() %>%
  write.table(file="clipboard-1000", sep="\t", col.names=F, row.names=F, quote=F)

#  Look at ATAC/H3K4me3 signal and expression at loop anchors ----------------------------

aged_signals = read.table(file.path(projdir,"track_signal_at_anchors","aged_RPKM_per_aged_loop_anchors.tab"), sep="\t",
                          col.names = c('chr','start','end','ATAC','H3K4me3'))
aged_peaks = read.table(file.path(projdir,"track_signal_at_anchors","aged_peaks_per_aged_loop_anchors.txt"), sep="\t",
                          col.names = c('chr','start','end','ATAC','H3K4me3'))
young_signals = read.table(file.path(projdir,"track_signal_at_anchors","young_RPKM_per_young_loop_anchors.tab"), sep="\t",
                           col.names = c('chr','start','end','ATAC','H3K4me3'))
young_peaks = read.table(file.path(projdir,"track_signal_at_anchors","young_peaks_per_young_loop_anchors.txt"), sep="\t",
                        col.names = c('chr','start','end','ATAC','H3K4me3'))

merged_signals = full_join(aged_signals, young_signals, by=c("chr","start","end"), suffix=c("_Aged","_Young")) %>% as_tibble()
merged_peaks = full_join(aged_peaks, young_peaks, by=c("chr","start","end"), suffix=c("_Aged","_Young")) %>% as_tibble()
merged_signals = full_join(merged_signals, merged_peaks, by=c("chr","start","end"), suffix=c("Signal","Peak")) %>% as_tibble()

## Add H3K4me3 and ATAC signals to loop anchors -----------------------------
classified_anchor_signals = bind_rows(left_join(classified_loop_df[,c("chr","start","x2","observed","expected_donut","domain_score",
                                                                    "Age","Type","left_anno_type","anno_type","loop_ID")], 
                                              merged_signals, 
                                              by=c("chr", "start", "x2"="end")) %>%
                                      dplyr::rename("end"="x2", "anchor_anno_type"="left_anno_type") %>%
                                      mutate(loop_strength = observed / expected_donut),
                                    left_join(classified_loop_df[,c("chr","y1","end","observed","expected_donut","domain_score",
                                                                    "Age","Type","right_anno_type","anno_type","loop_ID")], 
                                              merged_signals, 
                                              by=c("chr", "y1"="start", "end")) %>%
                                      dplyr::rename("start"="y1", "anchor_anno_type"="right_anno_type") %>%
                                      mutate(loop_strength = observed / expected_donut)) %>% 
  distinct()

loop_width = classified_loop_df %>%
  mutate(loop_width = end-start) %>%
  dplyr::select(loop_ID, loop_width)

classified_anchor_signals = left_join(classified_anchor_signals, loop_width, by="loop_ID")

classified_anchor_signals %>%
  #filter(Type %in% c("Intra-TAD","Loop Domain","TAD Boundary")) %>%
  pivot_longer(cols=c(contains("Peak")), names_to=c("Peak","Peak_Age"), names_pattern="(.+)_(.+)", values_to="Peak_num") %>%
  filter((Age=="Young" & grepl("Young",Peak_Age)) | (Age=="Aged" & grepl("Aged",Peak_Age))) %>%
  group_by(Age, anchor_anno_type, Peak) %>%
  tally()

library(ggridges)

classified_anchor_signals %>%
  #filter(Type %in% c("Intra-TAD","Loop Domain","TAD Boundary")) %>%
  pivot_longer(cols=c(contains("Peak")), names_to=c("Peak","Peak_Age"), names_pattern="(.+)_(.+)", values_to="Peak_num") %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  dplyr::filter((Age=="Young" & grepl("Young",Peak_Age)) | (Age=="Aged" & grepl("Aged",Peak_Age))) %>%
  dplyr::filter(Peak=="ATAC") %>%
  ggplot(aes(x=Peak_num)) +
  facet_wrap(~anchor_anno_type) +
  geom_density_ridges(aes(y=Age, fill=Age), scale=0.9,
                      jittered_points = TRUE,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7)

classified_anchor_signals %>%
  #filter(Type %in% c("Intra-TAD","Loop Domain","TAD Boundary")) %>%
  pivot_longer(cols=c(contains("Peak")), names_to=c("Peak","Peak_Age"), names_pattern="(.+)_(.+)", values_to="Peak_num") %>%
  mutate(Age = factor(Age, levels=c("Young","Aged"))) %>%
  dplyr::filter((Age=="Young" & grepl("Young",Peak_Age)) | (Age=="Aged" & grepl("Aged",Peak_Age))) %>%
  dplyr::filter(Peak=="ATAC") %>%
  ggplot(aes(x=Peak_num)) +
  facet_wrap(~anchor_anno_type) +
  #facet_grid(rows=vars(Age), cols=vars(anchor_anno_type), scales="free_y") +
  #facet_grid(rows=vars(Peak), cols=vars(anchor_anno_type), scales="free") +
  geom_histogram(aes(fill=Age), color="black", binwidth=2, position="dodge") +
  #geom_density(aes(fill=Age), alpha=0.75, position="dodge", adjust=2, kernel="gaussian") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_x_continuous(expand=expansion(mult=c(0.01,0))) +
  labs(x="# ATAC Sites / Loop Anchor", y="# ATAC Sites") +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() + 
  theme(legend.position = "top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        strip.text = element_text(size=12, face="bold", color="black"))
ggsave(file.path(projdir, "Figures", "num_ATAC_peaks_per_anchor.png"), dpi=300, width=4.5, height=3)

## Add expression to anchors -------------------------------------------

anchor_genes = eval(R_bedtools_intersect(tss_slop_1kb, GRanges(classified_anchor_signals), wa=T, wb=T))
classified_anchor_expression = cbind(second(anchor_genes) %>% 
                                       as.data.frame() %>% 
                                       dplyr::select("seqnames","start","end","loop_ID","Age",
                                                     "anchor_anno_type", "anno_type", "Type", "domain_score", "observed", "expected_donut"),
                                     first(anchor_genes) %>% 
                                       as.data.frame() %>% 
                                       dplyr::select("gene"))

classified_anchor_expression_df = classified_anchor_expression %>%
  dplyr::filter(gene %in% rownames(avg.tpm)) %>%
  distinct() %>%
  arrange(loop_ID)

classified_anchor_expression_df = cbind(classified_anchor_expression_df, avg.tpm[classified_anchor_expression_df$gene, ])

classified_anchor_expression_df = classified_anchor_expression_df %>%
  mutate(loop_strength = observed / expected_donut,
         Age = factor(Age, levels=c("Young","Aged"))) %>%
  pivot_longer(cols=c(A,Y), names_to="TPM_Age", values_to="TPM") %>%
  dplyr::filter((TPM_Age=="A" & Age=="Aged") | (TPM_Age=="Y" & Age=="Young")) %>%
  dplyr::filter(Type != "Asymmetric-TAD")

saveRDS(classified_anchor_expression_df, file.path(projdir, "classified_anchor_expression_df.RDS"))
classified_anchor_expression_df = readRDS(file.path(projdir, "classified_anchor_expression_df.RDS"))

classified_anchor_expression_test = classified_anchor_expression_df %>%
  filter(anno_type != "Distal-Distal" & anchor_anno_type == "Promoter") %>%
  mutate(log2TPM = log2(TPM)) %>%
  group_by(anno_type) %>%
  wilcox_test(log2TPM ~ Age)
classified_anchor_expression_test = add_xy_position(classified_anchor_expression_test, x="anno_type", group="Age", dodge=0.9)

classified_anchor_expression_df %>%
  filter(anno_type != "Distal-Distal" & anchor_anno_type == "Promoter") %>%
  mutate(log2TPM = log2(TPM)) %>%
  ggplot(aes(x=anno_type, y=log2TPM)) + 
  geom_boxplot(aes(fill=Age), color="black", outlier.size=0.2, position=position_dodge(0.9)) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  stat_pvalue_manual(classified_anchor_expression_test, tip.length=0.01, bracket.nudge.y = 1) +
  theme_bw() + 
  labs(x="Loop Type", y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position="top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir,"Figures","loop_anchor_expression.png"), dpi=300, width=4, height=5)

classified_anchor_expression_df %>%
  group_by(Age, anno_type, anchor_anno_type) %>%
  summarise(num_genes = n_distinct(gene))

# get tables for GSEA
classified_anchor_expression_df %>% 
  filter(anchor_anno_type=="Promoter" & Age=="Young" & anno_type=="Promoter-Distal") %>% 
  distinct() %>% drop_na() %>%
  mutate(log2fc=log2(A/Y)) %>% 
  dplyr::select(gene,log2fc) %>% 
  write.table(file="clipboard-1000",sep="\t",row.names=F,col.names=F,quote=F)

classified_anchor_expression_df %>% 
  filter(anchor_anno_type=="Promoter" & Age=="Aged" & anno_type=="Promoter-Distal") %>%
  pull(gene) %>% unique() %>% writeClipboard()

classified_anchor_expression_df %>% 
  filter(anchor_anno_type=="Promoter" & Age=="Aged" & anno_type=="Promoter-Distal") %>%
  dplyr::select(seqnames,start,end) %>%
  distinct() %>% head() %>% mutate(width=end-start)

## Annotate promoter end of Promoter-Distal loops 
# loop_anchor_genes_df = left_join(classified_anchor_expression_df[,c('loop_ID','gene','anno_type')] %>% 
#                                    filter(anno_type == "Promoter-Distal"), 
#                                  classified_anchor_signals %>% 
#                                    filter(anno_type == "Promoter-Distal") %>%
#                                    dplyr::select(c('chr','start','end','Age','Type','loop_strength','domain_score','anchor_anno_type','loop_ID')), 
#                                  by='loop_ID') %>%
#   filter(anchor_anno_type=="Promoter" & gene %in% rownames(avg.tpm)) %>%
#   distinct()
# loop_anchor_genes_df = cbind(loop_anchor_genes_df, avg.tpm[loop_anchor_genes_df$gene, ])
# loop_anchor_genes_df = loop_anchor_genes_df %>%
#   pivot_longer(cols=c(A,Y), names_to="TPM_Age", values_to="TPM") %>%
#   filter((grepl("Aged",loop_ID) & TPM_Age=="A") | (grepl("Young",loop_ID) & TPM_Age=="Y"))

## Compare track signals at distal end of promoter-distal loops with local gene expression --------------
promoter_distal_genes_df = left_join(classified_anchor_expression_df %>% filter(anno_type == "Promoter-Distal"), 
                                     classified_anchor_signals %>% 
                                       filter(anno_type == "Promoter-Distal" & anchor_anno_type=="Distal") %>%
                                       dplyr::select(loop_ID, contains("ATAC"), contains("H3K4me3")),
                                     by="loop_ID") 

promoter_distal_genes_df %>%
  filter(Type %in% c("Intra-TAD","Loop Domain","TAD Boundary")) %>%
  pivot_longer(cols=c(contains("Signal")), names_to="signal_type_age", values_to="signal") %>%
  filter((Age=="Young" & grepl("Young",signal_type_age)) | (Age=="Aged" & grepl("Aged",signal_type_age))) %>%
  drop_na() %>% distinct() %>%
  mutate(signal_type = sapply(signal_type_age, function(x) ifelse(grepl("ATAC",x), "ATAC", "H3K4me3"))) %>%
  #filter(signal_type == "ATAC") %>%
  arrange(loop_ID) %>%
  mutate(log2signal = log2(signal),
         log2TPM = log2(TPM),
         log2strength = log2(loop_strength),
         TPM_percentile = rank(TPM)/n(),
         Age = factor(Age, levels=c("Young","Aged"))) %>%
  ggscatter(x="log2signal", y="log2TPM", color="Age", facet.by=c("signal_type","Age"), add = "reg.line", cor.coef = T)

promoter_distal_genes_df %>%
  #filter(Type %in% c("Intra-TAD","Loop Domain","TAD Boundary")) %>%
  pivot_longer(cols=contains("Signal"), names_to="signal_type_age", values_to="signal") %>%
  filter((Age=="Young" & grepl("Young",signal_type_age)) | (Age=="Aged" & grepl("Aged",signal_type_age))) %>%
  drop_na() %>% distinct() %>%
  mutate(signal_type = sapply(signal_type_age, function(x) ifelse(grepl("ATAC",x), "ATAC", "H3K4me3"))) %>%
  #filter(signal_type == "ATAC") %>%
  #arrange(loop_ID) %>%
  mutate(log2signal = log2(signal),
         log2TPM = log2(TPM),
         log2strength = log2(loop_strength),
         TPM_percentile = rank(TPM)/n(),
         Age = factor(Age, levels=c("Young","Aged"))) %>%
  group_by(Age,signal_type) %>%
  mutate(log2signal_bins = cut_width(rank(log2signal)/n(), width = 0.25, boundary=0, closed="left")) %>%
  ggplot(aes(x=log2signal_bins, y=log2TPM)) +
  facet_grid(rows=vars(signal_type)) +
  stat_summary(fun.data="mean_se", geom="path", aes(group=interaction(Age, log2signal_bins), color=Age)) +
  stat_summary(fun.data="mean_se", geom="pointrange", aes(color=Age))

  
  ggboxplot(x="log2signal_bins", y="TPM_percentile", fill="Age", facet.by=c("signal_type","Age"), outlier.size=0.2) +
  labs(x="RPKM Percentile Bins at Promoter-Distal Anchors", y="Local Gene Expression (TPM) Percentile") +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold", hjust=1, angle=35),
        axis.text.y = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")
ggsave(file.path(projdir, "Figures", "loop_anchor_track_signalRPKM_1kb.png"), dpi=300, width=4.5, height=4)

## Compare track signals at distal end of promoter-distal loops with local track signals ---------------------

classified_anchor_signals %>%
  filter(anno_type=="Distal-Distal") %>%
  dplyr::select("Age", anno_type, anchor_anno_type, contains("Signal"), "loop_ID") %>%
  arrange(loop_ID) %>%
  mutate(anchor_orientation = rep(c("left","right"), times=length(Age)/2)) %>%
  pivot_longer(cols=contains("Signal"), 
             names_to=c("Signal",".value"), names_pattern="(.+)_(.+)") 

classified_anchor_signals_plt_df = classified_anchor_signals %>%
  filter(anno_type=="Promoter-Distal") %>%
  dplyr::select("Age", anno_type, anchor_anno_type, loop_strength, loop_width, contains("Signal"), "loop_ID") %>%
  pivot_longer(cols=contains("Signal"), 
               names_to=c("Signal",".value"), names_pattern="(.+)_(.+)") %>%
  pivot_wider(names_from="anchor_anno_type", values_from=c("AgedSignal","YoungSignal")) %>%
  pivot_longer(cols=contains(c("YoungSignal","AgedSignal")),
               names_to=c("Signal_Age",".value"), names_pattern="(.+)_(.+)") %>%
  dplyr::filter(sapply(1:nrow(.), function(i) grepl(Age[i], Signal_Age[i]))) %>%
  drop_na() %>% distinct() %>%
  arrange(loop_ID) %>%
  group_by(Age, Signal) %>%
  mutate(promoter_signal_bin = cut_width(rank(Promoter)/n(), width = 0.2, center = 0.5),
         distal_signal_bin = cut_width(rank(Distal)/n(), width = 0.2, center = 0.5),
         log2promoter_signal = log2(Promoter),
         log2distal_signal = log2(Distal),
         promoter_signal_percentile = rank(Promoter)/n(),
         distal_signal_percentile = rank(Distal)/n(),
         Age = factor(Age, levels=c("Young","Aged")))

plot(classified_anchor_signals_plt_df$log2distal_signal[grepl("Aged", classified_anchor_signals_plt_df$Signal_Age) & classified_anchor_signals_plt_df$Signal=="ATAC"],
     classified_anchor_signals_plt_df$log2promoter_signal[grepl("Aged", classified_anchor_signals_plt_df$Signal_Age) & classified_anchor_signals_plt_df$Signal=="H3K4me3"])
plot(classified_anchor_signals_plt_df$log2distal_signal[grepl("Young", classified_anchor_signals_plt_df$Signal_Age) & classified_anchor_signals_plt_df$Signal=="ATAC"],
     classified_anchor_signals_plt_df$log2promoter_signal[grepl("Young", classified_anchor_signals_plt_df$Signal_Age) & classified_anchor_signals_plt_df$Signal=="H3K4me3"])


ggboxplot(classified_anchor_signals_plt_df, 
          x="distal_signal_bin", y="promoter_signal_percentile", fill="Age", facet.by=c("Signal","Age"), outlier.size=0.2) +
  labs(x="RPKM Percentile Bins at Promoter-Distal Anchors", 
       y="RPKM Percentiles at Local Anchors") +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold", hjust=1, angle=35),
        axis.text.y = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")
ggsave(file.path(projdir, "Figures", "loop_anchor_track_RPKM_proximal_distal.png"), dpi=300, width=4.5, height=4)

# pointrange within-age comparison
classified_anchor_signals_test = classified_anchor_signals_plt_df %>%
  group_by(Signal, Age) %>%
  wilcox_test(log2promoter_signal ~ distal_signal_bin, ref="all")
classified_anchor_signals_test = add_xy_position(classified_anchor_signals_test, scales="free_y",
                                                 x="distal_signal_bin", group="Age", dodge=0.25, fun="mean_se")
classified_anchor_signals_test$y.position[classified_anchor_signals_test$Signal=="H3K4me3" & classified_anchor_signals_test$Age=="Aged"] = 
  classified_anchor_signals_test$y.position[classified_anchor_signals_test$Signal=="H3K4me3" & classified_anchor_signals_test$Age=="Aged"] - 0.35

ggplot(classified_anchor_signals_plt_df, aes(x=distal_signal_bin, y=log2promoter_signal)) +
  facet_grid(rows=vars(Signal), scales="free_y") +
  stat_summary(fun.data="mean_se", aes(color=Age, group=interaction(Age,Signal)), geom="path", position=position_dodge(0.25)) +
  stat_summary(fun.data="mean_se", aes(color=Age), geom="pointrange", size=0.5, fatten=3, position=position_dodge(0.25)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.15)), breaks=seq(4,6,0.5)) +
  stat_pvalue_manual(classified_anchor_signals_test, remove.bracket = T) +
  scale_color_manual(values=rev(pal_nejm()(2))) +
  labs(x="RPKM Percentile at Promoter-Distal Anchors", 
       y=expression(bold(paste("lo","g"["2"],"(RPKM) at Local Anchors")))) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold", hjust=1, angle=35),
        axis.text.y = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "loop_anchor_track_RPKM_proximal_distal_points.png"), dpi=300, width=4, height=4)

# pointrange across-age comparison
classified_anchor_signals_age_test = classified_anchor_signals_plt_df %>%
  group_by(Signal, distal_signal_bin) %>%
  wilcox_test(log2promoter_signal ~ Age)
classified_anchor_signals_age_test = add_xy_position(classified_anchor_signals_age_test, scales="free_y", 
                                                     x="distal_signal_bin", dodge=0, fun="mean_se")
classified_anchor_signals_age_test = mutate(classified_anchor_signals_age_test, y.position=y.position+0.05)

ggplot(classified_anchor_signals_plt_df, aes(x=distal_signal_bin, y=log2promoter_signal)) +
  facet_grid(rows=vars(Signal), scales="free_y") +
  stat_summary(fun.data="mean_se", aes(color=Age, group=interaction(Age,Signal)), geom="path") +
  stat_summary(fun.data="mean_se", aes(color=Age), geom="pointrange", size=0.5, fatten=3) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.15)), breaks=seq(3,7,0.4)) +
  stat_pvalue_manual(classified_anchor_signals_age_test, remove.bracket = T, label.size = 3.5) +
  scale_color_manual(values=rev(pal_nejm()(2))) +
  labs(x="RPKM Percentile at Promoter-Distal Anchors", 
       y=expression(bold(paste("lo","g"["2"],"(RPKM) at Local Anchors")))) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text.x = element_text(color="black", size=12, face="bold", hjust=1, angle=35),
        axis.text.y = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "loop_anchor_track_RPKM_proximal_distal_points_withinAge.png"), dpi=300, width=4, height=4)

## Linear regression on promoter signal -----------------------
signal_model = lm(log2promoter_signal ~ log2distal_signal, 
                  subset = Age=="Young" & Signal=="ATAC", data=classified_anchor_signals_plt_df)
summary(signal_model)

ggplot(classified_anchor_signals_plt_df, aes(x=log2distal_signal, y=log2promoter_signal)) +
  facet_grid(cols=vars(Age), rows=vars(Signal), scales="free_y") +
  geom_point(aes(color=Age), alpha=0.5, size=1) +
  scale_x_continuous(limits=c(0,10)) +
  scale_y_continuous(limits=c(0,10)) +
  stat_smooth(method="lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))
ggsave(file.path(projdir, "Figures", "loop_anchor_track_RPKM_proximal_distal_scatter.png"), dpi=300, width=6, height=4)

classified_anchor_signals_plt_df %>%
  group_by(Age,Signal) %>%
  arrange(loop_strength, .by_group=T) %>%
  ggplot(aes(x=log2distal_signal, y=log2promoter_signal)) +
  facet_grid(cols=vars(Age), rows=vars(Signal), scales="free_y") +
  geom_point(aes(color=loop_strength), alpha=0.5, size=1) +
  scale_color_viridis_c(option = "A") +
  #scale_color_distiller(palette="viridis") +
  scale_x_continuous(limits=c(-1,10)) +
  scale_y_continuous(limits=c(-1,10)) +
  stat_smooth(method="lm") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_bw() + guides(color=guide_colorbar(title="Loop Contact Strength")) +
  labs(x=expression(bold(paste("lo","g"["2"],"(RPKM) at Distal Anchors"))),
       y=expression(bold(paste("lo","g"["2"],"(RPKM) at Promoter Anchors")))) +
  theme(axis.text = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "loop_anchor_track_RPKM_proximal_distal_scatter.png"), dpi=300, width=6, height=5)


# Expression of annotated loop anchors ------------------------------------

classified_loop_df %>%
  mutate(loop_strength = observed/expected_bottom_left) %>%
  ggplot(aes(x=Age, y=log10(loop_strength))) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2)

classified_loop_df %>%
  ggplot(aes(x=Age, y=domain_score)) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2)

loop_domain_plt = classified_loop_df %>%
  dplyr::filter(Type!="Asymmetric-TAD") %>% 
  ggplot(aes(x=Age, y=domain_score)) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2) +
  theme_bw() +
  guides(fill=guide_legend(title="Anchor Type")) +
  scale_fill_npg() +
  labs(y="Loop Domain Score", x=NULL) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        strip.text = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        plot.margin = margin(0,0,0,0),
        legend.position = "top")

expression_plt = ggplot(classified_anchor_expression_df, aes(x=Age, y=log2(TPM))) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2) +
  theme_bw() +
  guides(fill=guide_legend(title="Anchor Type")) +
  scale_fill_npg() +
  labs(y=expression(bold(paste("lo","g"["2"],"(TPM)"))), x=NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        strip.text = element_text(size=12, face="bold"),
        plot.margin = margin(0,0,0,0),
        legend.position = "none")

loop_strength_plt = classified_loop_df %>%
  dplyr::filter(Type!="Asymmetric-TAD") %>%
  mutate(loop_strength = observed/expected_bottom_left) %>% 
  ggplot(aes(x=Age, y=log2(loop_strength))) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2) +
  theme_bw() +
  guides(fill=guide_legend(title="Anchor Type")) +
  scale_fill_npg() +
  labs(y=expression(bold(paste("lo","g"["2"],"(Loop Strength)"))), x=NULL) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        strip.text = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        plot.margin = margin(0,0,0,0),
        legend.text = element_text(size=12))

max_dims <- patchwork::get_max_dim(loop_strength_plt, loop_domain_plt, expression_plt)
loop_strength_plt_aligned <- patchwork::set_dim(loop_strength_plt, max_dims)
loop_domain_plt_aligned <- patchwork::set_dim(loop_domain_plt, max_dims)

#gridExtra::grid.arrange(loop_domain_plt, expression_plt, loop_strength_plt)
ggarrange(plotlist=list(expression_plt, loop_strength_plt_aligned, loop_domain_plt_aligned), 
          nrow=3, common.legend=T, legend="top", legend.grob = get_legend(loop_strength_plt)) 
ggsave(file.path(projdir,"Figures","loop_strength_expression_plt.png"), dpi=300, width=5, height=8)

library(ggpmisc)
ggplot(classified_anchor_expression_df, aes(x=log2(loop_strength), y=domain_score)) +
  facet_grid(rows=vars(anno_type)) +
  geom_point(aes(color=Age), alpha=0.5) +
  stat_smooth(aes(group=Age, color=Age), method="lm") +
  stat_poly_eq(formula = y~x, 
               aes(group=Age, color=Age, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  theme_bw()

# Distance between loop anchors -------------------------------------------

classified_loop_df %>%
  group_by(anno_type, Age) %>%
  tally()

loop_width_test = classified_loop_df %>%
  mutate(log10width = log10(end-start)) %>%
  group_by(anno_type) %>%
  wilcox_test(log10width ~ Age)
loop_width_test = add_xy_position(loop_width_test, x="anno_type", group="Age", dodge=0.9)
loop_width_test = loop_width_test %>% mutate(plabel = sapply(p, function(x) ifelse(x<2.2e-16, "<2.2e-16", x)))

classified_loop_df %>%
  mutate(width = end-start,
         Age = factor(Age, levels=c("Young","Aged"))) %>%
  ggplot(aes(x=anno_type, y=log10(width))) +
  geom_violin(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  geom_boxplot(aes(group=interaction(Age,anno_type)), fill="white", width=0.1, color="black", position=position_dodge(0.9), outlier.shape=NA) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1)), breaks=seq(4,8,0.5)) +
  stat_pvalue_manual(loop_width_test, label="plabel", size=4, tip.length=0.01) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme_bw() + 
  labs(x="Loop Type", y=expression(bold(paste("Distance Between Loop Anchors [lo","g"["10"],"(bp)]")))) +
  theme(axis.text = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"Figures","distance_between_loop_anchors_by_group.png"), dpi=300, width=5, height=5)

classified_loop_df %>%
  mutate(width = end-start) %>%
  ggplot(aes(x=Age, y=log10(width))) +
  geom_boxplot(aes(fill=anno_type), color="black", outlier.size=0.2, position=position_dodge(0.9)) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1)), breaks=seq(4,8,0.5)) +
  # stat_pvalue_manual(loop_width_test, label="plabel", size=4, tip.length=0.01,
  #                    step.increase = 0.1, bracket.nudge.y = 0.1, step.group.by = "Age") +
  scale_fill_npg() +
  theme_bw() + guides(fill = guide_legend(title="Anchor\nType", nrow=3)) +
  labs(x=NULL, y=expression(bold(paste("Distance Between Loop Anchors [lo","g"["10"],"(bp)]")))) +
  #labs(x=NULL, y=expression(bold(paste("lo","g"["10"],"[Distance Between Loop Anchors (bp)]")))) +
  theme(axis.text = element_text(color="black", size=12, face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.margin = margin(0,0,0,-40))
ggsave(file.path(projdir,"Figures","distance_between_loop_anchors.png"), dpi=300, width=3.5, height=4.5)

# Test --------------------------------------------------------------------


classified_anchor_signals %>%
  filter(anno_type=="Promoter-Distal") %>%
  pivot_longer(cols=c(ATAC_Aged, ATAC_Young), names_to="ATAC_Age", values_to="ATAC") %>%
  pivot_longer(cols=c(H3K4me3_Aged, H3K4me3_Young), names_to="H3K4me3_Age", values_to="H3K4me3") %>%
  filter((Age=="Young" & grepl("Young",ATAC_Age) & grepl("Young",H3K4me3_Age)) | (Age=="Aged" & grepl("Aged",ATAC_Age) & grepl("Aged",H3K4me3_Age))) %>%
  drop_na() %>% distinct() %>%
  pivot_longer(cols=c(ATAC, H3K4me3), names_to="signal_type", values_to="signal") %>%
  arrange(loop_ID) %>%
  mutate(log2signal = log2(signal)) %>%
  ggboxplot(x="Age", y="log2signal", fill="anchor_anno_type", facet.by="signal_type")

  # pivot_wider(names_from=anchor_anno_type, values_from=signal) %>%
  # ggscatter(x="Promoter", y="Distal", color="Age", facet.by="signal_type")

# Look at ATAC and H3K4me3 at Promoter-Distal loop anchors
classified_anchor_signals_plt_df = classified_anchor_signals %>%
  pivot_longer(cols=c(ATAC_Aged, ATAC_Young), names_to="ATAC_Age", values_to="ATAC") %>%
  pivot_longer(cols=c(H3K4me3_Aged, H3K4me3_Young), names_to="H3K4me3_Age", values_to="H3K4me3") %>%
  filter((Age=="Young" & grepl("Young",ATAC_Age) & grepl("Young",H3K4me3_Age)) | (Age=="Aged" & grepl("Aged",ATAC_Age) & grepl("Aged",H3K4me3_Age))) %>%
  drop_na() %>% distinct() %>%
  pivot_longer(cols=c(ATAC, H3K4me3), names_to="signal_type", values_to="signal") %>%
  mutate(Age = factor(Age, levels=c("Young","Aged")),
         log2signal = log2(signal)) %>%
  ungroup()

classified_anchor_signals_test = classified_anchor_signals_plt_df %>%
  group_by(Age, signal_type) %>%
  wilcox_test(log2signal ~ anno_type) 
classified_anchor_signals_test = classified_anchor_signals_test %>% 
  add_xy_position(x="Age", dodge=0.9) %>%
  mutate(padj_label = sapply(p.adj, function(x) ifelse(x<2.2e-16, "<2.2e-16", x)))

ggboxplot(classified_anchor_signals_plt_df, x="Age", y="log2signal", 
          facet.by = "signal_type", fill="anno_type", outlier.size=0.2) +
  stat_pvalue_manual(classified_anchor_signals_test, tip.length=0.01, label="padj_label", step.increase = 0.05, step.group.by = c("signal_type", "Age")) +
  scale_y_continuous(breaks=seq(-2,12,2), expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_npg() +
  theme_bw() +
  guides(fill = guide_legend(title="Anchor\nType", nrow=2)) +
  labs(y=expression(bold(paste("lo","g"["2"],"(RPKM)"))), x=NULL) +
  theme(legend.position = "top",
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold", color="black"),
        strip.text = element_text(size=12, color="black", face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.box.margin = margin(0,0,0,-40))
ggsave(file.path(projdir, "Figures", "loop_anchor_signals.png"), dpi=300, width=5, height=5)
