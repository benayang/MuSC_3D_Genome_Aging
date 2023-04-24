library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpattern)
library(rstatix)
library(ggsci)
library(ggpubr)
library(latex2exp)

projdir="C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\track_coverage\\fold_change"

tasks=c("aged_scores_per_A_to_B","aged_scores_per_B_to_A",
        "young_scores_per_A_to_B","young_scores_per_B_to_A",
        "aged_scores_per_aged.A","aged_scores_per_aged.B",
        "young_scores_per_young.A","young_scores_per_young.B",
        "young_scores_per_StaticA","young_scores_per_StaticB",
        "aged_scores_per_StaticA","aged_scores_per_StaticB")

data=list()
for(taskname in tasks){
  tmp = read.table(file.path(projdir,paste0(taskname,".tab")),header=F,col.names=c("chr","start","end","ATAC","H3K4me3"))
  if(grepl("A_to_B",taskname) | grepl("B_to_A",taskname)) {
    splitname=unlist(strsplit(taskname,split="_",fixed=T))
    fname = paste(paste(splitname[1:3],collapse="_"),paste(splitname[4:6],collapse=""),sep = "_")  
  } else {
    fname = taskname
  }
  print(taskname)
  print(fname)
  tmp$Group = fname
  tmp = tmp %>% separate(Group,into=c("Age",NA,NA,"Group"),sep="_",remove=F)
  data[[fname]] = tmp
}

df=do.call(rbind,data) %>% drop_na()


# Look at compartment switch ----------------------------------------------

comp.switch.test = df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  group_by(Group,BW) %>%
  wilcox_test(FC~Age) %>%
  add_xy_position(x="Group",group="Age",dodge=0.9,fun="mean_se")

comp.switch.win.age.test = df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  group_by(Age,BW) %>%
  wilcox_test(FC~Group) %>%
  add_xy_position(x="Group",group="Age",fun="mean_se",dodge=0.9)

df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  group_by(Group,BW,Age) %>%
  summarise(count=n())

df %>% 
  pivot_longer(cols=c("ATAC","H3K4me3"),names_to="BW",values_to="FC") %>%
  filter(Group %in% c("AtoB","BtoA","StaticA","StaticB")) %>%
  mutate(Age=factor(Age,levels=c("young","aged"))) %>%
  ggplot(aes(x=Group,y=FC)) + 
  stat_summary(aes(group=Age), fun.data="mean_se",geom="errorbar", width=0.5, position=position_dodge(0.9)) +
  stat_summary(aes(fill=Age), fun.data="mean_se",color="black",geom="bar",position=position_dodge(0.9)) +
  #geom_boxplot(aes(fill=Age), color="black", outlier.size=0.1) +
  facet_grid(rows=vars(BW), scales="free_y") +
  scale_fill_manual(values=rev(pal_nejm()(2)), labels=c("Y","A")) +
  stat_pvalue_manual(comp.switch.test, tip.length=0.01, label.size=5) +
  #stat_pvalue_manual(comp.switch.win.age.test, tip.length=0.01, step.group.by = "BW", step.increase = 0.05, bracket.nudge.y = 0.5) +
  scale_x_discrete(labels=c(TeX("$A \\rightarrow B$ ",bold=T),TeX("$B \\rightarrow A$ ",bold=T),TeX("Static A ",bold=T),TeX("Static B ",bold=T))) +
  scale_y_continuous(expand = expansion(mult=c(0,0.2))) +
  theme_bw() +
  labs(x=NULL, y="Fold Change") +
       #y=expression(bold(paste("lo","g"["2"],"(RPKM + 1)")))) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"compartment_switch_fold_change.png"),dpi=300,width=4.25,height=5)


# Look at fold change in compartment switch -------------------------------

fc_df = df %>% 
  distinct() %>% 
  mutate(coord=paste(chr,start,end,sep="_")) %>% 
  pivot_wider(id_cols=c(coord,Group), names_from=Age, values_from = c(H4K20me1,ATAC,H3K27me3,H3K4me3)) %>% 
  mutate(H4K20me1=(H4K20me1_aged/H4K20me1_young),
         ATAC=(ATAC_aged/ATAC_young),
         H3K27me3=(H3K27me3_aged/H3K27me3_young),
         H3K4me3=(H3K4me3_aged/H3K4me3_young)) %>%
  drop_na() %>%
  mutate(Group=factor(Group),
         coord=factor(coord))

idx=which(apply(fc_df[,c("H4K20me1","ATAC","H3K27me3","H3K4me3")],1,function(x) all(is.finite(x))))
fc_df = fc_df[idx,] %>%
  mutate(H4K20me1=log2(H4K20me1),
         ATAC=log2(ATAC),
         H3K27me3=log2(H3K27me3),
         H3K4me3=log2(H3K4me3))
# fc_df.test = fc_df %>%
#   filter(Group %in% c("AtoB","BtoA","static")) %>%
#   select(Group,H4K20me1,ATAC,H3K27me3,H3K4me3) %>%
#   pivot_longer(cols=c(H4K20me1,ATAC,H3K27me3,H3K4me3), names_to="Mark", values_to="L2FC") %>%
#   ungroup() %>%
#   mutate(Group=factor(Group,levels=c("static","BtoA","AtoB")),
#          Mark=factor(Mark)) %>%
#   group_by(Mark) %>%
#   wilcox_test(L2FC~Group)

fc_df.test = data.frame(
  H3K27me3=c(wilcox.test(pull(subset(fc_df,Group=="AtoB","H3K27me3")), pull(subset(fc_df,Group=="static","H3K27me3")))$p.value,
             wilcox.test(pull(subset(fc_df,Group=="BtoA","H3K27me3")), pull(subset(fc_df,Group=="static","H3K27me3")))$p.value,
             wilcox.test(pull(subset(fc_df,Group=="AtoB","H3K27me3")), pull(subset(fc_df,Group=="BtoA","H3K27me3")))$p.value),
  ATAC=c(wilcox.test(pull(subset(fc_df,Group=="AtoB","ATAC")), pull(subset(fc_df,Group=="static","ATAC")))$p.value,
         wilcox.test(pull(subset(fc_df,Group=="BtoA","ATAC")), pull(subset(fc_df,Group=="static","ATAC")))$p.value,
         wilcox.test(pull(subset(fc_df,Group=="AtoB","ATAC")), pull(subset(fc_df,Group=="BtoA","ATAC")))$p.value),
  H3K4me3=c(wilcox.test(pull(subset(fc_df,Group=="AtoB","H3K4me3")), pull(subset(fc_df,Group=="static","H3K4me3")))$p.value,
            wilcox.test(pull(subset(fc_df,Group=="BtoA","H3K4me3")), pull(subset(fc_df,Group=="static","H3K4me3")))$p.value,
            wilcox.test(pull(subset(fc_df,Group=="AtoB","H3K4me3")), pull(subset(fc_df,Group=="BtoA","H3K4me3")))$p.value),
  H4K20me1=c(wilcox.test(pull(subset(fc_df,Group=="AtoB","H4K20me1")), pull(subset(fc_df,Group=="static","H4K20me1")))$p.value,
             wilcox.test(pull(subset(fc_df,Group=="BtoA","H4K20me1")), pull(subset(fc_df,Group=="static","H4K20me1")))$p.value,
             wilcox.test(pull(subset(fc_df,Group=="AtoB","H4K20me1")), pull(subset(fc_df,Group=="BtoA","H4K20me1")))$p.value),
  row.names = c("static.AtoB","static.BtoA","AtoB.BtoA"))

plt = fc_df %>%
  filter(Group %in% c("AtoB","BtoA","static")) %>%
  mutate(Group=factor(Group,levels=c("static","BtoA","AtoB"))) %>%
  select(Group,H4K20me1,ATAC,H3K27me3,H3K4me3) %>%
  pivot_longer(cols=c(H4K20me1,ATAC,H3K27me3,H3K4me3), names_to="Mark", values_to="L2FC") %>%
  ggplot(aes(x=Group,y=L2FC)) +
  facet_wrap(~Mark, scales = "free_y") +
  #stat_summary(fun.data="mean_se",geom="errorbar", width=0.5) +
  #stat_summary(aes(fill=Group), fun.data="mean_se",geom="bar") +
  geom_boxplot(aes(fill=Group), outlier.size=0.01, color="black", position=position_dodge(0.9)) +
  scale_y_continuous(expand = expansion(mult=c(0.05,0.1))) +
  theme_bw() + 
  #geom_bracket(xmin=1, xmax=2, y.position=-0.2, label=fc_df.test[1,"ATAC"]) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="ATAC"]),
  #              xmax=with(plt.data.df, x[Group=="BtoA" & Mark=="ATAC"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="AtoB" & Mark=="ATAC"])+0.5,
  #              label=7.63e-8, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="ATAC"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="ATAC"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="ATAC"])+0.5,
  #              label=7.26e-05, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="BtoA" & Mark=="ATAC"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="ATAC"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="ATAC"])+2,
  #              label=4.39e-04, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H3K27me3"]),
  #              xmax=with(plt.data.df, x[Group=="BtoA" & Mark=="H3K27me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="AtoB" & Mark=="H3K27me3"])+0.5,
  #              label=8.41e-8, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H3K27me3"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H3K27me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H3K27me3"])+1.5,
  #              label=1.23e-3, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="BtoA" & Mark=="H3K27me3"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H3K27me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H3K27me3"])+3,
  #              label=3.17e-5, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H3K4me3"]),
  #              xmax=with(plt.data.df, x[Group=="BtoA" & Mark=="H3K4me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="AtoB" & Mark=="H3K4me3"])+0.5,
  #              label=0.210, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H3K4me3"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H3K4me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H3K4me3"])+0.5,
  #              label=0.0529, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="BtoA" & Mark=="H3K4me3"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H3K4me3"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H3K4me3"])+2,
  #              label=0.883, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H4K20me1"]),
  #              xmax=with(plt.data.df, x[Group=="BtoA" & Mark=="H4K20me1"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="AtoB" & Mark=="H4K20me1"])+0.5,
  #              label=0.502, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="AtoB" & Mark=="H4K20me1"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H4K20me1"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H4K20me1"])+0.5,
  #              label=0.346, tip.length = 0.01) +
  # geom_bracket(xmin=with(plt.data.df, x[Group=="BtoA" & Mark=="H4K20me1"]),
  #              xmax=with(plt.data.df, x[Group=="static" & Mark=="H4K20me1"]),
  #              y.position=with(plt.data.df,ymax_final[Group=="static" & Mark=="H4K20me1"])+2,
  #              label=0.921, tip.length = 0.01) +
  #scale_fill_jama(labels=c(TeX("Static"),TeX("$B \\rightarrow A$"),TeX("$A \\rightarrow B$"))) +
  labs(x=NULL,
       y=expression(bold(paste("RPKM lo","g"["2"],"Fold-Change (Aged/Young)")))) +
  guides(fill=guide_legend(title=NULL)) +
  theme(axis.text.y=element_text(size=14,color="black",face="bold"),
        axis.text.x=element_blank(),
        axis.title=element_text(size=14,color="black",face="bold"),
        axis.ticks.x = element_blank(),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14))
plt
#scale_x_discrete(labels=c(TeX("Static",bold=T),TeX("$B \\rightarrow A$", bold=T),TeX("$A \\rightarrow B$", bold=T))) +
ggsave(file.path(projdir,"compartment_switch_foldchange_barplot.png"),dpi=300,width=4,height=5)
#ggsave(file.path(projdir,"compartment_switch_foldchange_boxplots.png"),dpi=300,width=6.5,height=5)

plt.data = ggplot_build(plt)
plt.data.df = as.data.frame(plt.data[[1]][1])
plt.data.df$Group=rep(c("AtoB","BtoA","static"),each=5)
plt.data.df$Mark=rep(c("ATAC","H3K27ac","H3K27me3","H3K4me3","H4K20me1"),times=3)
bracket.coord = as.data.frame(plt.data[[1]][1])$x
ymax.list =  as.data.frame(plt.data[[1]][1])$ymax_final


# Look at compartments only -----------------------------------------------

compartment.df = df %>% 
  pivot_longer(cols=c("H4K20me1","ATAC","H3K27me3","H3K4me3"),names_to="BW",values_to="RPKM") %>%
  filter(Group %in% c("young.A","young.B","aged.A","aged.B")) %>%
  separate(Group, into=c(NA,"Compartment"), sep="\\.", remove=F) %>%
  mutate(Age=factor(Age,levels=c("young","aged")))

comp.test = compartment.df %>% 
  mutate(log2RPKM=log2(RPKM+1)) %>%
  group_by(Age,BW) %>%
  rstatix::wilcox_test(log2RPKM~Compartment) %>%
  add_xy_position(x="Age",group="Compartment",dodge=0.9)
comp.test = comp.test %>% mutate(label="<2.2e-16")

comp.BW.test = compartment.df %>% 
  mutate(log2RPKM=log2(RPKM+1)) %>%
  group_by(BW) %>%
  rstatix::pairwise_wilcox_test(log2RPKM~Group)
comp.BW.test = comp.BW.test %>% add_xy_position(x="Age",group="BW",dodge=0.9)

with(compartment.df, t.test(RPKM[Group=="aged.A" & BW=="H3K27me3"],RPKM[Group=="aged.B" & BW=="H3K27me3"])$p.value)

compartment.df %>%
  group_by(Age,Compartment,BW) %>%
  summarise(avg=mean(log2(RPKM+1)),
            sem=sd(log2(RPKM+1))/n()) %>%
  ggplot(aes(x=Age,y=avg)) +
  facet_wrap(~BW) +
  geom_col_pattern(aes(fill=Age,pattern=Compartment), position=position_dodge(0.9),
                   color = 'black',
                   pattern_density = 0.1, 
                   pattern_fill    = 'white',
                   pattern_colour  = 'black') +
  geom_errorbar(aes(ymin=avg-sem, ymax=avg+sem, group=Compartment), position=position_dodge(0.9)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_pattern_manual(values=c("none","stripe")) +
  guides(fill="none") +
  theme_bw() +
  labs(x="",
       y=expression(bold(paste("lo","g"["2"],"(RPKM + 1)")))) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))

ggplot(compartment.df, aes(x=Age,y=log2(RPKM+1))) + 
  geom_boxplot_pattern(aes(fill=Age, pattern=Compartment), outlier.size=0.1,
                       color = 'black',
                       pattern_density = 0.1, 
                       pattern_fill    = 'white',
                       pattern_colour  = 'black',
                       position=position_dodge(0.9)) +
  facet_grid(cols=vars(BW)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_pattern_manual(values=c("none","stripe")) +
  #stat_pvalue_manual(comp.test, tip.length=0.01, label="label") +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_y_continuous(breaks=seq(0,12,3), expand = expansion(mult=c(0.05,0.1))) +
  guides(fill="none") +
  theme_bw() +
  labs(x="",
       y=expression(bold(paste("lo","g"["2"],"(RPKM + 1)")))) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
ggsave(file.path(projdir,"compartment.png"),dpi=300,width=8,height=5)


# compartment only fold change --------------------------------------------

fc_df_comp = compartment.df %>% 
  distinct() %>% 
  mutate(coord=paste(chr,start,end,sep="_")) %>% 
  pivot_wider(id_cols=c(coord,Age,Compartment), names_from=c(BW,Age), values_from = RPKM) %>% 
  mutate(H4K20me1=(H4K20me1_aged/H4K20me1_young),
         ATAC=(ATAC_aged/ATAC_young),
         H3K27me3=(H3K27me3_aged/H3K27me3_young),
         H3K4me3=(H3K4me3_aged/H3K4me3_young)) %>%
  drop_na()
idx=which(apply(fc_df_comp[,c("H4K20me1","ATAC","H3K27me3","H3K4me3")],1,function(x) all(is.finite(x))))
fc_df_comp = fc_df_comp[idx,] %>%
  mutate(H4K20me1=log2(H4K20me1),
         ATAC=log2(ATAC),
         H3K27me3=log2(H3K27me3),
         H3K4me3=log2(H3K4me3))

fc_df_comp %>%
  mutate(Compartment=factor(Compartment,levels=c("A","B"))) %>%
  select(Compartment,ATAC,H4K20me1,H3K27me3,H3K4me3) %>%
  pivot_longer(cols=c(ATAC,H4K20me1,H3K27me3,H3K4me3), names_to="Mark", values_to="L2FC") %>%
  ggplot(aes(x=Mark,y=L2FC)) +
  geom_boxplot(aes(fill=Compartment))
  geom_boxplot_pattern(aes(fill=Age, pattern=Compartment), outlier.size=0.1,
                       color = 'black',
                       pattern_density = 0.1, 
                       pattern_fill    = 'white',
                       pattern_colour  = 'black',
                       position=position_dodge(0.9)) 
  facet_grid(cols=vars(BW)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_pattern_manual(values=c("none","stripe")) +
  #stat_pvalue_manual(comp.test, tip.length=0.01, label="label") +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_y_continuous(breaks=seq(0,12,3), expand = expansion(mult=c(0.05,0.1))) +
  guides(fill="none") +
  theme_bw() +
  labs(x="",
       y=expression(bold(paste("lo","g"["2"],"(RPKM + 1)")))) +
  theme(axis.text=element_text(size=14,color="black",face="bold"),
        axis.title=element_text(size=14,color="black",face="bold"),
        strip.text=element_text(size=14,color="black",face="bold"),
        legend.position = "top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))