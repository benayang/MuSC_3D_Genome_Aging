library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)
library(ggsci)
library(VennDiagram)
library(grid)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS"

young=read.table(file.path(projdir,"young_merged_loops_noHeader.bedpe"),
                 col.names = c("chromosome1","x1","x2",
                               "chromosome2","y1","y2",
                               "blank1","blank2","blank3","blank4",
                               "color","observed",
                               "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                               "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                               "number_collapsed","centroid1","centroid2","radius"))
aged=read.table(file.path(projdir,"aged_merged_loops_noHeader.bedpe"),
                col.names = c("chromosome1","x1","x2",
                              "chromosome2","y1","y2",
                              "blank1","blank2","blank3","blank4",
                              "color","observed",
                              "expected_bottom_left","expected_donut","expected_horizontal","expected_vertical",
                              "fdr_bottom_left","fdr_donut","fdr_horizontal","fdr_vertical",
                              "number_collapsed","centroid1","centroid2","radius"))

young$Age="Y"
aged$Age="A"

dim(aged)
dim(inner_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2")))
dim(anti_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2")))
dim(anti_join(young[,1:6], aged[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2")))

write.table(inner_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))[,c("chromosome1","x1","y2")], 
                       file.path(projdir, "olap_stable_loop_domains.bed"), sep='\t', row.names=F, col.names=F, quote=F)
write.table(anti_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))[,c("chromosome1","x1","y2")], 
            file.path(projdir, "olap_aged_loop_domains.bed"), sep='\t', row.names=F, col.names=F, quote=F)
write.table(anti_join(young[,1:6], aged[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))[,c("chromosome1","x1","y2")], 
            file.path(projdir, "olap_young_loop_domains.bed"), sep='\t', row.names=F, col.names=F, quote=F)

plt.df = rbind(young,aged)
plt.df$Enrichment_bottom_left = plt.df$observed/plt.df$expected_bottom_left
plt.df$Age = factor(plt.df$Age, levels=c("Y","A"))

plt.df %>% group_by(Age) %>% summarise(enrichment=mean(log10(Enrichment_bottom_left)))

score.test = plt.df %>% 
  mutate(log10enrichment = log10(Enrichment_bottom_left)) %>%
  wilcox_test(log10enrichment ~ Age) %>% 
  add_xy_position(x="Age")

ggplot(plt.df, aes(x=Age,y=log10(Enrichment_bottom_left))) +
  geom_violin(aes(fill=Age),color="black") +
  geom_boxplot(width=0.2, fill="white", color="black", outlier.shape=NA) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  theme_bw() +
  #stat_pvalue_manual(score.test, tip.length = 0.01, label.size = 5, bracket.nudge.y = 0.1) +
  labs(x="",y=expression(bold(paste("lo","g"["10"],"(Bottom Left Enrichment)")))) +
  theme(legend.position="none",
        axis.text=element_text(color="black",size=14,face="bold"),
        axis.title=element_text(color="black",size=14,face="bold"))
ggsave(file.path(projdir,"loop_enrichment_violin.png"),dpi=300,width=3,height=4)

venn.plot <- draw.pairwise.venn(6252+65, 6252+38, 6252, c("Young", "Aged"), fill = pal_nejm()(2), margin=0.05,
                                cex = 3, cat.cex = 3, cat.fontfamily = 'Sans', fontfamily = 'Sans');
png(file.path(projdir, "Figures","loop_venn_diagram.png"), res=300, units="in", width=5, height=5)
grid.newpage()
grid.draw(venn.plot)
dev.off()

venn.plot <- draw.pairwise.venn(844+2033, 844+2634, 844, c("Aged", "Young"), fill = pal_nejm()(2), margin=0.05, cat.pos = 180, cat.dist = 0.05,
                                cex = 2, cat.cex = 2, cat.fontfamily = 'Sans', fontfamily = 'Sans');
png(file.path(projdir, "Figures","olap_loop_venn_diagram.png"), res=300, units="in", width=5, height=5)
grid.newpage()
grid.draw(venn.plot)
dev.off()


# Intersect with compartment ----------------------------------------------

library(HelloRanges)

stable = inner_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
unique_aged = anti_join(aged[,1:6], young[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))
unique_young = anti_join(young[,1:6], aged[,1:6], by=c("chromosome1","x1","x2","chromosome2","y1","y2"))

stable_bed = stable[ ,c("chromosome1","x1","y2")] %>% setNames(c("chr","start","end")) %>% GRanges()
unique_aged_bed = unique_aged[ ,c("chromosome1","x1","y2")] %>% setNames(c("chr","start","end")) %>% GRanges()
unique_young_bed = unique_young[ ,c("chromosome1","x1","y2")] %>% setNames(c("chr","start","end")) %>% GRanges()

ab.switch = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression\\compartmentBed\\100kb",
                                 "ab.switch.bed"), sep="\t", col.names=c('chrom','start','end','group'))
ab.switch = GRanges(ab.switch)

aged.ab.bed = GRanges(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\without KR normalization",
                                           "aged.ab_100kb.bed"), sep='\t', col.names=c('chrom','start','end','compartment')))
young.ab.bed = GRanges(read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\without KR normalization",
                                            "young.ab_100kb.bed"), sep='\t', col.names=c('chrom','start','end','compartment')))

R_bedtools_intersect(stable_bed, ab.switch, f=1)

eval({
  pairs <- findOverlapPairs(stable_bed, aged.ab.bed, ignore.strand = TRUE)
  olap = width(first(pairs)) / width(second(pairs))
  keep = olap>=1
  ans <- pintersect(pairs, ignore.strand = TRUE)
  ans$group = second(pairs)$group
  ans = ans[keep,]
})

eval({
  pairs <- findOverlapPairs(unique_aged_bed, aged.ab.bed, ignore.strand = TRUE)
  olap = width(first(pairs)) / width(second(pairs))
  keep = olap>=1
  ans <- pintersect(pairs, ignore.strand = TRUE)
  ans$group = second(pairs)$group
  ans = ans[keep,]
  table(ans$group)
})

prop.table(table(eval({
  pairs <- findOverlapPairs(unique_aged_bed, ab.switch, ignore.strand = TRUE)
  ans <- pintersect(pairs, ignore.strand = TRUE)
  ans$group = second(pairs)$group
})))

prop.table(table(eval({
  pairs <- findOverlapPairs(unique_young_bed, ab.switch, ignore.strand = TRUE)
  ans <- pintersect(pairs, ignore.strand = TRUE)
  ans$group = second(pairs)$group
})))

# Intersect loops with TADs -----------------------------------------------

young_domains = young %>% dplyr::select(chromosome1, x1, y2) %>% distinct() %>% setNames(c("chr","start","end")) %>% GRanges()
aged_domains = aged %>% dplyr::select(chromosome1, x1, y2) %>% distinct() %>% setNames(c("chr","start","end")) %>% GRanges()

aged.domain = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                   "aged.merged","40kb","aged.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))
young.domain = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer",
                                    "young.merged","40kb","young.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))
aged_TADs = aged.domain %>% dplyr::select(1:3) %>% setNames(c("chr","start","end")) %>% GRanges()
young_TADs = young.domain %>% dplyr::select(1:3) %>% setNames(c("chr","start","end")) %>% GRanges()

#young_overlap = R_bedtools_intersect(a=young_domains, b=young_TADs, g="mm10", wo=T, f=1)
young_overlap = eval({
  genome <- Seqinfo(genome = "mm10")
  gr_a <- young_domains
  gr_b <- young_TADs
  pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE, 
                            type = "within")
  olap <- pintersect(pairs, ignore.strand = TRUE)
  keep <- width(olap)/width(first(pairs)) == 1
  pairs <- pairs[keep]
  ans <- pairs
  mcols(ans)$overlap_width <- width(olap)[keep]
  ans
})

aged_overlap = eval({
  genome <- Seqinfo(genome = "mm10")
  gr_a <- aged_domains
  gr_b <- aged_TADs
  pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE, 
                            type = "within")
  olap <- pintersect(pairs, ignore.strand = TRUE)
  keep <- width(olap)/width(first(pairs)) == 1
  pairs <- pairs[keep]
  ans <- pairs
  mcols(ans)$overlap_width <- width(olap)[keep]
  ans
})

length(aged_overlap) / length(aged_domains)
length(young_overlap) / length(young_domains)


# Differential loops ------------------------------------------------------

young_aged_overlap = eval({
  genome <- Seqinfo(genome = "mm10")
  gr_a <- young_domains
  gr_b <- aged_domains
  pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE, 
                            type = "within")
  olap <- pintersect(pairs, ignore.strand = TRUE)
  keep <- width(olap)/width(first(pairs)) == 1
  pairs <- pairs[keep]
  ans <- pairs
  mcols(ans)$overlap_width <- width(olap)[keep]
  ans
})
