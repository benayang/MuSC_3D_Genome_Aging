library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(ggsci)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\hiccups_diff\\loop_domain_track_coverage\\olap_analysis"

files = list.files(file.path(projdir), pattern="loop.tab")

data.list = list()
for(f in files) {
  fname = gsub(".tab","",f,fixed=T)
  print(fname)
  data.list[[fname]] = read.table(file.path(projdir, f), sep='\t', header=F, col.names=c("chr","start","end","ATAC","H3K4me3"))
  data.list[[fname]]$fname = fname
  data.list[[fname]]$score_age = ifelse(grepl("aged_scores",fname,fixed=T), "Aged", "Young")
  data.list[[fname]]$loop_group = ifelse(grepl("stable_loop",fname,fixed=T), "Static", ifelse(grepl("young_loop",fname,fixed=T),"Young","Aged"))
}

data_df = bind_rows(data.list) %>% mutate(group = paste(score_age, loop_group, sep="_")) %>% distinct()
table(data_df$fname, data_df$score_age, data_df$loop_group)
table(data_df$group)

plt_data = data_df %>% 
  pivot_longer(cols=c(ATAC,H3K4me3), names_to="mod") %>%
  pivot_wider(id_cols = c(chr,start,end,loop_group), names_from = c(mod, score_age), values_from = value) %>%
  mutate(H3K4me3 = log2(H3K4me3_Young/H3K4me3_Aged),
         ATAC = log2(ATAC_Young/ATAC_Aged)) %>%
  pivot_longer(cols=c(H3K4me3, ATAC), names_to = "log2fc") 

plt_data_test = plt_data %>% group_by(log2fc) %>% wilcox_test(value ~ loop_group)
plt_data_test = plt_data_test %>% add_xy_position(x="loop_age", dodge=0.9)

ggplot(plt_data, aes(x=loop_group, y=H3K4me3_Young)) +
  geom_boxplot(aes(fill=loop_group))

ggplot(plt_data, aes(x=loop_group, y=value)) +
  facet_wrap(~log2fc) +
  geom_boxplot(aes(fill=loop_group), outlier.size=0.1, color="black", position = position_dodge(0.9)) +
  #stat_pvalue_manual(plt_data_test, tip.length = 0.01) +
  scale_y_continuous(expand=expansion(mult=c(0.1,0.1))) +
  theme_bw() + labs(x="Loop Age", y=expression(bold(paste("lo","g"["2"]," Signal Fold-Change (Young/Aged)")))) +
  theme(legend.position="top",
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        axis.text = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "foldchange_boxplot.png"), dpi=300, width=5, height=5)
