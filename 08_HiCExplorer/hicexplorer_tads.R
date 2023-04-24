library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggsci)
library(Hmisc)
library(reshape2)
library(scales)
library(rstatix)

# Create filename dataframe -----------------------------------------------

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\08_HiCExplorer"

young.domain.files = grep("bedpe",list.files(file.path(projdir, "young.merged"), pattern = "_domains.bed", recursive=T, full.names = T),fixed=T,invert=T,value=T)
aged.domain.files = grep("bedpe",list.files(file.path(projdir, "aged.merged"), pattern = "_domains.bed", recursive=T, full.names = T),fixed=T,invert=T,value=T)

file.df = data.frame(files = c(young.domain.files, aged.domain.files))
file.df = file.df %>%
  separate(files, into=c(NA,"Age","Resolution","Filename"), sep="/", remove=F) %>%
  mutate(FDR = as.numeric(gsub(grep("thresh", 
                                          unlist(strsplit(Filename,"_",fixed=T)), 
                                          value = T),
                                     pattern = "thresh",
                                     replacement = "")),
         Age = sapply(Age, function(x) capitalize(gsub(".merged","",x,fixed=T)))) %>%
  filter(Resolution != "10kb" & FDR  != 0.001)
file.df$Age = recode(file.df$Age, Young="Y", Aged="A") 
file.df$Age = factor(file.df$Age, levels=c("Y","A"))

# plot number and size of TADs per group ----------------------------------

tad.files = lapply(file.df$files, function(x) read.table(x, sep="\t", as.is=T))
names(tad.files) = file.df$Filename
tad.files = do.call(rbind, tad.files)
tad.files = tad.files[,c(1:5,9)]
colnames(tad.files) = c("chr","start","end","ID","score","color")
tad.files$Filename = rownames(tad.files)
tad.files = tad.files %>% 
  mutate(Filename = sapply(Filename, function(x) unlist(strsplit(x, ".bed", fixed=T))[1]))
tad.files$Filename = sapply(tad.files$Filename, function(x) paste0(x,".bed"))
tad.files$Age = factor(file.df[match(tad.files$Filename, file.df$Filename),"Age"])
tad.files$Resolution = factor(file.df[match(tad.files$Filename, file.df$Filename),"Resolution"],
                              levels=c("40kb","100kb","250kb","500kb"))
tad.files$FDR = factor(file.df[match(tad.files$Filename, file.df$Filename),"FDR"])
tad.files$size = with(tad.files, end-start)

tad.files %>%
  group_by(Age, Resolution, FDR) %>%
  summarise(count = n(),
            med.size = median(size)/1e6) %>%
  pivot_longer(cols=c(count,med.size)) %>%
  ggplot(aes(x=Age, y=value, fill=FDR)) +
  facet_grid(cols=vars(Resolution),rows=vars(name),scales="free_y",labeller=labeller(name=c(count="TAD Count",med.size="Median TAD Size (Mb)"))) +
  geom_col(position=position_dodge(0.9)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_fill_npg() +
  theme_bw() + labs(x=NULL,y=NULL) +
  theme(legend.position = "top",
        legend.title = element_text(size=12,face="bold"),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold", color="black"),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "TAD_num_and_size_allRes.png"), dpi=300, width=5, height=5)

  tad.files %>%
  ggplot(aes(x=Age, y=value, fill=FDR)) +
  facet_grid(cols=vars(Resolution),rows=vars(name),scales="free_y") +
  geom_col(position=position_dodge(0.9)) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  scale_fill_npg() +
  theme_bw() + xlab("") +
  theme(legend.position = "top",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold", color="black"),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))

ggplot(tad.files, aes(x=Age, y=log10(size), fill=FDR)) +
  facet_grid(cols=vars(Resolution)) +
  geom_boxplot(color="black", outlier.size = 0.2) +
  scale_y_continuous(breaks=seq(4,8,0.5)) +
  scale_fill_npg() +
  #scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw() + xlab("") +
  theme(legend.position = "top",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold", color="black"),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "median_TAD_size_allRes.png"), dpi=300, width=6, height=5)

tad.files %>%
  filter(FDR==0.01, Resolution=="100kb") %>%
  pivot_wider(id_cols=c("chr","start","end"), names_from="Age", values_from="score") %>%
  ggplot(aes(x=A,y=Y)) +
  #geom_point(alpha=0.6, color="red") +
  geom_hex(bins=60) +
  geom_abline(linetype="dashed") +
  # geom_smooth() +
  # stat_poly_eq(formula = y~x, 
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #              parse = TRUE) +        
  scale_fill_distiller(palette = "RdYlBu") +
  theme_bw() + 
  coord_fixed(ratio = 1, xlim = c(-1.5,1), ylim = c(-1.5,1)) +
  scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  scale_y_continuous(breaks=seq(-1.5,1.5,0.5)) +
  xlab("TAD separation score in Aged") +
  ylab("TAD separation score in Young") +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "TAD_insulation_scores_100kb_0.01FDR.png"), dpi=300, width=5, height=5)

# get insulation scores at TAD boundaries ---------------------------------

# Make boundary file dataframe
young.boundary.files = list.files(file.path(projdir, "young.merged"), pattern = "_boundaries.bed", recursive=T, full.names = T)
aged.boundary.files = list.files(file.path(projdir, "aged.merged"), pattern = "_boundaries.bed", recursive=T, full.names = T)

# create table of boundary files and associated metadata
bound.file.df = data.frame(files = c(young.boundary.files, aged.boundary.files))
bound.file.df = bound.file.df %>%
  separate(files, into=c(NA,"Age","Resolution","Filename"), sep="/", remove=F) %>%
  mutate(FDR = as.numeric(gsub(grep("thresh", 
                                    unlist(strsplit(Filename,"_",fixed=T)), 
                                    value = T),
                               pattern = "thresh",
                               replacement = "")),
         Age = sapply(Age, function(x) capitalize(gsub(".merged","",x,fixed=T)))) %>%
  dplyr::filter(Resolution != "10kb" & FDR  != 0.001)
bound.file.df$Age = recode(bound.file.df$Age, Young="Y", Aged="A") 
bound.file.df$Age = factor(bound.file.df$Age, levels=c("Y","A"))

# read in and process boundary files
boundary.files = lapply(bound.file.df$files, function(x) read.table(x, sep="\t", as.is=T))
names(boundary.files) = bound.file.df$Filename
boundary.files = do.call(rbind, boundary.files)
colnames(boundary.files) = c("chr","start","end","ID","score","strand")
boundary.files$Filename = rownames(boundary.files)
# remove numeric and .bed suffix
boundary.files = boundary.files %>% 
  mutate(Filename = sapply(Filename, function(x) unlist(strsplit(x, ".bed", fixed=T))[1]))
  # add .bed suffix back in
boundary.files$Filename = sapply(boundary.files$Filename, function(x) paste0(x,".bed"))
# get age, resolution, and FDR from file.df
boundary.files$Age = factor(bound.file.df[match(boundary.files$Filename, bound.file.df$Filename),"Age"])
boundary.files$Resolution = factor(bound.file.df[match(boundary.files$Filename, bound.file.df$Filename),"Resolution"],
                              levels=c("40kb","100kb","250kb","500kb"))
boundary.files$FDR = factor(bound.file.df[match(boundary.files$Filename, bound.file.df$Filename),"FDR"])
boundary.files$size = with(boundary.files, end-start)

boundary.files %>%
  group_by(Age, Resolution, FDR) %>%
  summarise(count = n(),
            med.size = median(size))

# select 40kb resolution
boundary.files.df = boundary.files %>%
  filter(FDR==0.01, Resolution=="40kb") %>%
  pivot_wider(id_cols=c("chr","start","end"), names_from="Age", values_from="score") %>%
  mutate(group = factor(apply(.[,c("Y","A")], 1, function(x) if(is.na(x[1])){"A"} else if (is.na(x[2])) {"Y"} else {"Shared"}),
                        levels=c("Y","A","Shared")))
levels(boundary.files.df$group) = list(Young="Y", Aged="A", Shared="Shared")
#boundary.files.df$Y[is.na(boundary.files.df$Y)] = 0
#boundary.files.df$A[is.na(boundary.files.df$A)] = 0

boundary.test = boundary.files %>%
  filter(FDR==0.01, Resolution=="40kb") %>%
  t_test(score~Age) %>%
  add_xy_position(x="Age")

boundary.files %>%
  filter(FDR==0.01, Resolution=="40kb") %>%
  ggplot(aes(x=Age,y=score)) +
  geom_violin(aes(fill=Age), color="black") +
  geom_boxplot(width=0.1, outlier.shape=NA, color="black") +
  stat_pvalue_manual(boundary.test, bracket.nudge.y = 0.1, tip.length = 0.01) +
  scale_y_continuous(breaks=seq(-1.5,1.5,0.5), expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw() +
  ylab("TAD boundary strength") + xlab("") +
  theme(legend.position = "none",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "TAD_separation_scores_40kb_0.01FDR_violin.png"), width=2.5, height=4.5, dpi=300)

# try violin plots for unique and shared TAD boundaries
boundary.files.df.long = boundary.files.df %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
  mutate(group = factor(group),
         Age = factor(Age, levels=c("Y","A")))
levels(boundary.files.df.long$group) = list(Y="Young", A="Aged", Shared="Shared")

table(boundary.files.df$group)

tmp = subset(boundary.files.df.long, Age=="Y" & (group=="Y" | group=="Shared"))
t.test(tmp$score~tmp$group)
tmp = subset(boundary.files.df.long, Age=="A" & (group=="A" | group=="Shared"))
t.test(tmp$score~tmp$group)
tmp = subset(boundary.files.df.long, group=="A" | group=="Y")
t.test(tmp$score~tmp$group)
tmp = subset(boundary.files.df.long, group=="Shared")
t.test(tmp$score~tmp$Age)

boundary.files.df.long %>%
  ggplot(aes(x=Age,y=score)) +
  geom_violin(aes(fill=group), color="black") +
  geom_boxplot(aes(group=interaction(group,Age)), width=0.2, position=position_dodge(0.9), color="black", outlier.shape=NA) +
  geom_bracket(
    xmin = 0.775, xmax = 1.225, y.position = 1.05,
    label = "2.11e-11", tip.length = 0.01, label.size = 5,
  ) +
  geom_bracket(
    xmin = 1.775, xmax = 2.225, y.position = 1.6,
    label = "2.13e-6", tip.length = 0.01, label.size=5,
  ) +
  geom_bracket(
    xmin = 0.775, xmax = 1.775, y.position = 1.9,
    label = "0.303", tip.length = 0.01, label.size = 5,
  ) +
  geom_bracket(
    xmin = 1.225, xmax = 2.225, y.position = 2.2,
    label = "0.927", tip.length = 0.01, label.size=5,
  ) +
  scale_y_continuous(breaks=seq(-2,2.5,0.5),
                     expand=expansion(mult=c(0.05,0.1))) +
  scale_x_discrete(labels=c("Young","Aged")) +
  scale_fill_manual(values=c(pal_nejm()(2)[c(2,1)], "gray")) +
  theme_bw() +
  guides(fill=guide_legend(title="Boundary")) +
  ylab("TAD Separation Score") + xlab(NULL) +
  theme(legend.position="top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold"))
ggsave(file.path(projdir,"Figures","TAD_separation_scores_40kb_0.01FDR_groupedviolin.png"), width=3.5, height=4.5, dpi=300)

# try split violin plots for unique and shared TAD boundaries
vioplot.test = boundary.files.df %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
  #group_by(Age) %>%
  t_test(score~group)

library("vioplot")
png(file.path(projdir,"Figures","TAD_separation_scores_40kb_0.01FDR_splitViolin.png"), units="in", width=3.75, height=5.75, res=300)
boundary.files.df %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
filter(group=="Shared") %>%
  vioplot(score~Age, data=., col = "gray", plotCentre = "line", side = "right",
          ylim=c(-1.5,1.75), xlab="", ylab="TAD Separation Score",
          xaxt="n", yaxt="n", ann=F)
boundary.files.df %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
  filter(group!="Shared") %>%
  vioplot(score~Age, data=., col = pal_nejm()(2)[c(2,1)], plotCentre = "line", side = "left", add=T)
axis(side = 1, at = 1:2, labels = c("Y","A"), font = 2)
axis(side = 2, las = 2, font = 2, seq(-2,2,0.5))
mtext(side=2, text = "TAD separation score", font=2, line=2.75, cex=1.2)
legend("topright", fill = c(pal_nejm()(2)[c(2,1)], "gray"), legend = c("Young Unique", "Aged Unique","Shared"))
dev.off()

boundary.files.df %>%
  group_by(group) %>%
  summarise(Y_med = median(Y),
            A_med = median(A))

boundary.files.df %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
  ggplot(aes(x=Age,y=score)) +
  facet_grid(cols=vars(group)) +
  geom_violin(aes(fill=Age), color="black") +
  geom_boxplot(width=0.1, outlier.shape=NA, color="black") +
  stat_pvalue_manual(boundary.test, bracket.nudge.y = 0.1, tip.length = 0.01) +
  scale_y_continuous(breaks=seq(-1.5,1.5,0.5), expand=expansion(mult=c(0.1,0.1))) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw() +
  ylab("TAD boundary strength") + xlab("") +
  theme(legend.position = "none",
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "Figures", "TAD_insulation_scores_40kb_0.01FDR_violin.png"), width=2.5, height=4.5, dpi=300)

# write shared and unique TAD boundaries
boundary.files.df %>% filter(group=="Shared") %>% dplyr::select(chr, start, end) %>% write.table(file.path(projdir,"shared_TAD_boundary.bed"),quote=F,row.names=F,col.names=F,sep="\t")
boundary.files.df %>% filter(group=="Young") %>% dplyr::select(chr, start, end) %>% write.table(file.path(projdir,"unique_young_TAD_boundary.bed"),quote=F,row.names=F,col.names=F,sep="\t")
boundary.files.df %>% filter(group=="Aged") %>% dplyr::select(chr, start, end) %>% write.table(file.path(projdir,"unique_aged_TAD_boundary.bed"),quote=F,row.names=F,col.names=F,sep="\t")


# Get shared and unique TAD domains ---------------------------------------

aged.domain = read.table(file.path(projdir,"aged.merged","40kb","aged.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))
young.domain = read.table(file.path(projdir,"young.merged","40kb","young.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_domains.bed"))

domain.df = rbind(aged.domain,young.domain)
domain.df$Age = c(rep("Aged",nrow(aged.domain)), rep("Young",nrow(young.domain)))
colnames(domain.df) = c("chr","start","end","ID","score","strand","start2","end2","color","Age")

domain.df = domain.df %>% 
  pivot_wider(id_cols=c("chr","start","end"), names_from="Age", values_from=c("score","ID")) %>% 
  arrange(chr,start) %>%
  mutate(group = factor(apply(.[,c("score_Young","score_Aged")], 1, 
                              function(x) if(is.na(x[1])){"Aged"} else if (is.na(x[2])) {"Young"} else {"Shared"}),
                        levels=c("Young","Aged","Shared")))

domain.venn.data = data.frame(group = factor(names(table(domain.df$group)), levels=c("Shared","Young","Aged")),
                              fraction = matrix(prop.table(table(domain.df$group))),
                              count = matrix(table(domain.df$group)))


ggplot(domain.venn.data, aes(x="", y=fraction, fill=group)) +
  geom_col(width=1) +
  geom_text(aes(label = sprintf("%d\n(%s)", count, percent(fraction))), size=3, position = position_stack(vjust = 0.5)) +
  #facet_wrap(~Age) +
  coord_polar("y", start=0) +
  scale_y_continuous(labels = scales::percent) +
  xlab("") + ylab("") +
  theme_pubclean() +
  scale_fill_brewer(palette="Blues", labels=c("Common", "Unique Young", "Unique Aged")) +
  theme(axis.text = element_text(size=8, color="black"),
        legend.text = element_text(size=9),
        legend.title = element_blank())
ggsave(file.path(projdir,"Figures","unique_common_domains_piechart_40kb_0.01FDR.png"), dpi=300, width=4, height=4)

domain.df %>% filter(group=="Shared") %>% select(chr, start, end) %>% write.table(file.path(projdir,"shared_TAD_domain.bed"),quote=F,row.names=F,col.names=F,sep="\t")
domain.df %>% filter(group=="Young") %>% select(chr, start, end) %>% write.table(file.path(projdir,"unique_young_TAD_domain.bed"),quote=F,row.names=F,col.names=F,sep="\t")
domain.df %>% filter(group=="Aged") %>% select(chr, start, end) %>% write.table(file.path(projdir,"unique_aged_TAD_domain.bed"),quote=F,row.names=F,col.names=F,sep="\t")

# Plot insulation scores flanking boundaries ------------------------------

young.scores = read.table(file.path(projdir,"young.merged","40kb","young.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bedgraph"),
                          col.names = c("chr","start","end","score"))
aged.scores = read.table(file.path(projdir,"aged.merged","40kb","aged.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_score.bedgraph"),
                         col.names = c("chr","start","end","score"))

inner_join(young.scores,aged.scores,by=c("chr","start","end"),suffix=c(".young",".aged")) %>% glimpse()
anti_join(young.scores,aged.scores,by=c("chr","start","end"),suffix=c(".young",".aged")) %>% glimpse()
anti_join(aged.scores,young.scores,by=c("chr","start","end"),suffix=c(".young",".aged")) %>% glimpse()

boundary.files.df %>% filter(group=="Shared") %>% left_join(.,young.scores,by=c("chr","start","end"))

shared.boundary.files.df = boundary.files.df[boundary.files.df$group=="Shared",]
y.scores.df=list()
a.scores.df=list()
for(i in 1:nrow(shared.boundary.files.df)) {
  tmp_chr=shared.boundary.files.df[i,]$chr
  tmp_start=shared.boundary.files.df[i,]$start - 5*4e4
  tmp_end=shared.boundary.files.df[i,]$end + 4*4e4
  start_range=seq(tmp_start,tmp_end,4e4)
  tmp.df=data.frame(chr=tmp_chr,
                    start=start_range,
                    boundary=i)
  y.scores.df[[i]]=inner_join(young.scores,tmp.df,by=c("chr","start"))
  a.scores.df[[i]]=inner_join(aged.scores,tmp.df,by=c("chr","start"))
}
# only keep shared boundaries that have the same number of bins between young and aged
idx=unlist(lapply(y.scores.df,nrow)) == unlist(lapply(a.scores.df,nrow))
y.scores.df=y.scores.df[idx]
a.scores.df=a.scores.df[idx]
# only keep shared boundaries with 11 bins
y.scores.df=y.scores.df[unlist(lapply(y.scores.df,nrow))==11]
a.scores.df=a.scores.df[unlist(lapply(a.scores.df,nrow))==11]

y.scores.df=do.call(rbind,y.scores.df)
y.scores.df$Age="Young"
a.scores.df=do.call(rbind,a.scores.df)
a.scores.df$Age="Aged"

aggregate.scores.df = rbind(y.scores.df, a.scores.df)
aggregate.scores.df = aggregate.scores.df %>% group_by(Age,boundary) %>% mutate(bin=1:n())
aggregate.scores.df %>% 
  group_by(bin,Age) %>%
  summarise(avg=mean(score,na.rm=T),
            sem=sd(score)/sqrt(n())) %>%
  mutate(bin = factor(bin),
         Age=factor(Age,levels=c("Young","Aged"))) %>%
  ggplot(aes(x=bin,y=avg,group=Age)) +
  geom_path(aes(color=Age)) +
  geom_ribbon(aes(ymin=avg-sem, ymax=avg+sem, fill=Age), alpha=0.4, color=NA) +
  scale_x_discrete(labels=c("-200kb",rep("",4),"Boundary",rep("",4),"200kb")) +
  scale_color_manual(values=pal_nejm()(2)[c(2,1)]) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_bw() + 
  labs(x="Binned Distance From Boundary", y="Average TAD Separation Score") + 
  guides(color=guide_legend(ncol=2)) +
  theme(axis.text=element_text(size=14, face="bold", color="black"),
        axis.title=element_text(size=14, face="bold", color="black"),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_blank(),
        legend.background = element_rect(fill='transparent'),
        legend.position=c(0.5,0.95))
ggsave(file.path(projdir,"Figures","TAD_separation_score_flanking_boundary.png"),dpi=300,width=4.5,height=3.5)

# Try bookending boundary technique to classify TAD merging/split ---------

young.domain_bound = young.domain %>%
  dplyr::select(1:5) %>%
  setNames(c("chr","start","end","ID","Score")) %>%
  mutate(left = start+20000,
         right = end-20000) %>%
  mutate(left = paste(chr,start,left,sep="_"),
         right = paste(chr,right,end,sep="_")) %>%
  dplyr::select(-c("chr","start","end")) %>%
  pivot_longer(cols=c(left,right), names_to="boundary_name", values_to="boundary") %>%
  separate(boundary, into=c("chr","start","end"), sep="_", remove=F) %>% 
  mutate(start = as.integer(start), end = as.integer(end))
young.domain_bound[young.domain_bound$boundary_name=="left",]$start = young.domain_bound[young.domain_bound$boundary_name=="left",]$start - 20000
young.domain_bound[young.domain_bound$boundary_name=="right",]$end = young.domain_bound[young.domain_bound$boundary_name=="right",]$end + 20000
young.domain_bound = young.domain_bound %>% unite("full_boundary", c(chr,start,end) ,sep="_")

aged.domain_bound = aged.domain %>%
  dplyr::select(1:5) %>%
  setNames(c("chr","start","end","ID","Score")) %>%
  mutate(left = start+20000,
         right = end-20000) %>%
  mutate(left = paste(chr,start,left,sep="_"),
         right = paste(chr,right,end,sep="_")) %>%
  dplyr::select(-c("chr","start","end")) %>%
  pivot_longer(cols=c(left,right), names_to="boundary_name", values_to="boundary") %>%
  separate(boundary, into=c("chr","start","end"), sep="_", remove=F) %>% 
  mutate(start = as.integer(start), end = as.integer(end))
aged.domain_bound[aged.domain_bound$boundary_name=="left",]$start = aged.domain_bound[aged.domain_bound$boundary_name=="left",]$start - 20000
aged.domain_bound[aged.domain_bound$boundary_name=="right",]$end = aged.domain_bound[aged.domain_bound$boundary_name=="right",]$end + 20000
aged.domain_bound = aged.domain_bound %>% unite("full_boundary", c(chr,start,end) ,sep="_")

# if looking at the boundary scores files
# bookend.bound = dplyr::full_join(young.domain_bound, aged.domain_bound, by=c("full_boundary"), suffix=c(".Young",".Aged")) 
# bookend.bound = boundary.files.df %>% unite("full_boundary",c(chr,start,end)) %>% left_join(bookend.bound, ., by="full_boundary") %>%
#   select(-full_boundary, -boundary.Aged) %>% rename("boundary"=boundary.Young)
bookend.bound = dplyr::full_join(young.domain_bound, aged.domain_bound, by=c("boundary"), suffix=c(".Young",".Aged")) %>% 
  group_by(ID.Young) %>% 
  mutate(complete_loss = all(is.na(ID.Aged)), 
         left_loss = is.na(ID.Aged)[1] & !is.na(ID.Aged)[2], 
         right_loss = !is.na(ID.Aged)[1] & is.na(ID.Aged)[2]) %>%
  ungroup()
saveRDS(bookend.bound, file.path(projdir,"bookend.bound.RDS"))
bookend.bound = readRDS(file.path(projdir,"bookend.bound.RDS"))

# non_overlap_idx = which(is.na(bookend.bound[[ref]]))
# cur_query_idx = which(bookend.bound[[query]]==boundary_ids[b])
# pre_idx_diff = cur_query_idx[1] - non_overlap_idx
# post_idx_diff = non_overlap_idx - cur_query_idx[2]
# non_overlap_idx[sum(pre_idx_diff>=0)]
# non_overlap_idx[which(post_idx_diff>=0)[1]]

classify_TAD = function(bookend.bound, query, ref, shift_num_cap) {
  boundary_ids = unique(bookend.bound[[query]])
  left_loss_row_idx = which(bookend.bound$left_loss)
  right_loss_row_idx = which(bookend.bound$right_loss)
  
  bookend.bound_type_df = bookend.bound %>% dplyr::select(all_of(query)) %>% mutate(Type = NA) %>% as.data.frame()
  bookend.bound_ref_type_df = bookend.bound %>% dplyr::select(all_of(ref)) %>% mutate(Type = NA) %>% as.data.frame()
  
  b = 1
  while (b <= length(boundary_ids)) {
    tmp = bookend.bound[bookend.bound[[query]]==boundary_ids[b], ]
    tmp = tmp[which(!apply(tmp, 1, function(x) all(is.na(x)))), ]
    cur_query_row_idx = which(bookend.bound[[query]]==boundary_ids[b])
    if(nrow(tmp)==2) {
      if(!any(is.na(tmp[[ref]]))) {
        # If shared boundaries flank a single TAD in Young and Aged, it's shared.
        # Otherwise, it's split.
        bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_ids[b], ]$Type = ifelse(length(unique(tmp[[ref]]))==1, "Shared", "Split")

        ref.ID.nums = sapply(tmp$ID.Aged, function(x) as.numeric(unlist(strsplit(x, split="_", fixed=T))[3]))
        ref.ID.range = paste0("ID_0.01_", seq(ref.ID.nums[1], ref.ID.nums[2], 1))
        bookend.bound_ref_type_df[bookend.bound_ref_type_df[[ref]] %in% ref.ID.range, ]$Type = ifelse(length(unique(tmp[[ref]]))==1, "Shared", "Split")
        
        # tmp$Type = ifelse(length(unique(tmp[[ref]]))==1, "Shared", "Split")
      } else if (!is.na(tmp[[ref]][1]) & is.na(tmp[[ref]][2])) {
        # Merged TADs need to have matching outer boundaries of left-most and right-most TADs in Young compared to Aged.
        # At the first sign of a left non-overlapping boundary, scan forward to find the complement non-overlapping boundary.
        left_loss_row_idx_diff = left_loss_row_idx - cur_query_row_idx[2]
        right_loss_row_idx_diff = right_loss_row_idx - cur_query_row_idx[2]
        closest_left_loss_row_idx = left_loss_row_idx[which(left_loss_row_idx_diff>0)[1:2]]
        closest_right_loss_row_idx = right_loss_row_idx[which(right_loss_row_idx_diff>0)[1:2]]
        # There needs to be at least a closest_right_loss_row_idx for merged or shifted
        if (!is.na(closest_left_loss_row_idx[1]) & !is.na(closest_right_loss_row_idx[1])) {
          # The right-most TAD needs to be a left loss.
          if (closest_left_loss_row_idx[1] < closest_right_loss_row_idx[1]) {
            boundary_id_range = unique(bookend.bound[[query]][cur_query_row_idx[2]:closest_left_loss_row_idx[1]])
            if(all(bookend.bound[[ref]][closest_left_loss_row_idx] %in% tmp[[ref]])) {
              bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_id_range, ]$Type = "Merged"
            } else {
              # Need to set a cap on the number of contiguous TADs that can be classified as "Shifted"
              # THIS METHOD ALLOWS TOO MANY INDETERMINATE CLASSIFICATIONS. IF NUMBER OF TADS DOESN'T MATCH BETWEEN YOUNG AND AGED, IT'S INDETERMINATE.
              # bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_id_range, ]$Type = 
              #   ifelse((closest_left_loss_row_idx[2] - cur_query_row_idx[1] + 1)/2 <= shift_num_cap,"Shifted","Indeterminate")
              
              # Need to propagate to reference IDs
              ref.ID.nums = sapply(na.omit(bookend.bound[bookend.bound[[query]] %in% boundary_id_range, ]$ID.Aged), function(x) as.numeric(unlist(strsplit(x, split="_", fixed=T))[3]))
              ref.ID.range = paste0("ID_0.01_", seq(ref.ID.nums[1], ref.ID.nums[2], 1))
              
              bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_id_range, ]$Type = 
                ifelse(length(unique(boundary_id_range)) == length(unique(ref.ID.range)), "Shifted", "Combination")
              bookend.bound_ref_type_df[bookend.bound_ref_type_df[[ref]] %in% ref.ID.range, ]$Type = 
                ifelse(length(unique(boundary_id_range)) == length(unique(ref.ID.range)), "Shifted", "Combination")
              
              # bookend.bound_ref_type_df[bookend.bound_ref_type_df[[ref]] %in% ref.ID.range, ]$Type = 
              #   ifelse((closest_left_loss_row_idx[2] - cur_query_row_idx[1] + 1)/2 <= shift_num_cap,"Shifted","Indeterminate")
            }
            # Now skip the TADs that were classified as Merged or Shifted
            b = b + length(boundary_id_range)
            next
          }
        } else {
          bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_ids[b], ]$Type = "Combination"
        }
      } else {
        bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_ids[b], ]$Type = "Combination"
      }
    } else {
      bookend.bound_type_df[bookend.bound_type_df[[query]] %in% boundary_ids[b], ]$Type = "Huh?"
    }
    b = b + 1
  }
  bookend.bound$Type = bookend.bound_type_df$Type
  bookend.bound = bookend.bound %>% dplyr::select(-contains("loss"))
  bookend.bound = bookend.bound %>% 
    separate(ID.Aged, into=c(NA,NA,"ID.Aged.num"), sep="_", remove=F) %>%
    mutate(ID.Aged.num = as.numeric(ID.Aged.num))

  # get runs of Type
  type_runs = rle(bookend.bound$Type)
  type_end = cumsum(type_runs$lengths)
  type_start = c(1, lag(type_end)[-1] + 1)
  type_runs_df = as_tibble(data.frame(Type=type_runs$values, start=type_start, end=type_end))
  
  # Now look for Huh? and reclassify as flanking classifications
  #huh_idx = type_runs_df$start[which(type_runs_df$Type=="Huh?")] : type_runs_df$end[which(type_runs_df$Type=="Huh?")]
  huh_idx = which(bookend.bound$Type == "Huh?")
  bookend.bound[huh_idx, ]$Type = bookend.bound_ref_type_df[match(bookend.bound[huh_idx, ]$ID.Aged, bookend.bound_ref_type_df$ID.Aged), ]$Type
 
  # Second pass to weed out false "Merged" or "Shifted" classifications since these should be at least 2 (2 consecutive TADs)
  complex_idx = type_runs_df %>% filter(Type %in% c("Merged","Shifted")) %>% mutate(length = end-start) %>% filter(length==1) %>% pull(start)
  bookend.bound[c(complex_idx, complex_idx+1),]$Type = "Combination"
  
  # Classify any remaining NA types as Indeterminate
  bookend.bound[is.na(bookend.bound$Type), ]$Type = "Combination"
  
  return(list(bookend.bound, type_runs_df))
}

bookend.bound_classified = classify_TAD(bookend.bound, "ID.Young", "ID.Aged", shift_num_cap=4)
bookend.bound_df = bookend.bound_classified[[1]]
type_runs_df = bookend.bound_classified[[2]]

# note that this doesn't mean that distinct contiguous "shifted" TADs can't exist
png(file.path(projdir,"Figures","shifted_TAD_histogram_without_clipping.png"),res=300,width=7,height=7,units="in")
type_runs_df %>% filter(Type == "Shifted") %>% mutate(length = end-start) %>% arrange(desc(length)) %>% mutate(num_tad = (length+1)/2) %>% pull(num_tad) %>% 
  hist(breaks=seq(0,20,1), main="# Contiguous Young TADs in 'Shifted' category without clipping", plot=T)
dev.off()

# saveRDS(bookend.bound_df, file.path(projdir,"bookend.bound_df.RDS"))
bookend.bound_df = readRDS(file.path(projdir,"bookend.bound_df.RDS"))
# aged_vs_young_bookend.bound_df = classify_TAD(bookend.bound, "ID.Aged", "ID.Young")

bookend_type_dist = bookend.bound_df %>%
  #dplyr::filter(Type!="Huh?" & !is.na(ID.Young)) %>%
  pivot_longer(cols=c(ID.Young, ID.Aged), names_to="ID_age", values_to="ID") %>%
  dplyr::filter(!is.na(ID)) %>%
  group_by(Type, ID_age) %>% 
  summarise(num=n()/2) %>%
  ungroup() %>%
  group_by(ID_age) %>%
  mutate(fraction = num/sum(num)) %>%
  arrange(ID_age, num) %>%  
  ungroup() %>%
  mutate(Type=factor(Type, levels = c("Combination", "Split", "Merged", "Shifted", "Shared")),
         ID_age = factor(ID_age, levels=c("ID.Young","ID.Aged")))
  
bookend_colors = rev(hcl.colors(5,"Zissou 1"))
names(bookend_colors) = c("Combination", "Split", "Merged", "Shifted", "Shared")

bookend_type_dist$Colors = bookend_colors[bookend_type_dist$Type]

# plot type distributions
# ggplot(bookend_type_dist, aes(x="",y=fraction,fill=Type)) +
#   geom_col(width=1, color="black", position=position_stack()) +
#   geom_text(aes(label = sprintf("%d\n(%s)", num, percent(fraction))), color="white", size=3, position = position_stack(vjust = 0.5)) +
#   coord_polar("y",start=0) +
#   labs(x=NULL, y=NULL) +
#   theme_pubclean() +
#   scale_fill_uchicago() +
#   theme(axis.text = element_blank(),
#         axis.ticks.length = unit(0,"in"),
#         legend.text = element_text(size=11),
#         legend.title = element_blank())
# ggsave(file.path(projdir,"Figures","young_TAD_classification.png"),dpi=300,width=5,height=5)
ggplot(bookend_type_dist, aes(y=ID_age,x=num,fill=Type)) +
  geom_col(color="black", position=position_stack()) +
  # geom_text(aes(label = sprintf("%d (%s)", num, percent(fraction))), color="white", size=4, position = position_stack(vjust = 0.5)) +
  #geom_text(aes(label = scales::percent(fraction,accuracy=0.1)), color="white", size=4, position = position_stack(vjust = 0.5)) +
  scale_y_discrete(labels=c("Young","Aged")) +
  scale_x_continuous(expand=expansion(mult=c(0,0.1)), breaks=seq(0,3000,500)) +
  labs(y=NULL, x="# of TADs") +
  theme_bw() +
  scale_fill_manual(values=rev(hcl.colors(5,"Zissou 1"))) +
  guides(fill=guide_legend(nrow=2)) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.ticks.length.x = unit(0,"in"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.text = element_text(size=11),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"Figures","TAD_classification_bar.png"),dpi=300,height=2.5,width=4.5)
# ggsave(file.path(projdir,"Figures","TAD_classification_bar.png"),dpi=300,width=3.5,height=3.75)

# plot boundary strength distributions
scale_this = function(x) { return((x-mean(x,na.rm=T))/sd(x,na.rm=T)) }

bookend.bound_plt_df = bookend.bound_df %>%
  drop_na() %>%
  dplyr::select(-contains("boundary_name"), -"boundary") %>%
  distinct() %>%
  dplyr::mutate(Type = factor(Type, levels=levels(bookend_type_dist$Type))) %>%
  dplyr::mutate(diff = Score.Aged - Score.Young) 
bookend.bound_plt_df$diff_zscore = scale_this(bookend.bound_plt_df$diff)
bookend.bound_plt_df %>%
  mutate(Enriched_In = ifelse(diff>0,"Aged","Young")) %>%
  mutate(Differential = factor(ifelse(abs(diff_zscore)>2, "Differential", "Non-Differential"), levels=c("Differential","Non-Differential"))) %>%
  ggplot(aes(x=Score.Young,y=Score.Aged)) +
  facet_wrap(~Type) +
  geom_abline(slope=1,intercept=0,lwd=0.5,color="gray",lty=2) +
  geom_point(aes(color=Differential), size=0.1, alpha=0.6) +
  scale_y_continuous(limits=c(-1.75,1.5), breaks=seq(-2,2,1)) +
  scale_x_continuous(limits=c(-1.75,1.5), breaks=seq(-2,2,1)) +
  scale_color_manual(values=c("red","black")) +
  coord_equal() +
  theme_bw() + labs(x="Young TAD Boundary Score", y="Aged TAD Boundary Score") +
  guides(color=guide_legend(override.aes = list(size=3))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=11, face="bold"))
ggsave(file.path(projdir,"Figures","TAD_classification_boundary_scores_scatter.png"),dpi=300, width=5,height=5)

# are genes at changed boundaries or TADs interesting?
# sig_young.domain = bookend.bound[bookend.bound$ID.Young %in% (bookend.bound_plt_df %>% dplyr::filter(abs(diff_zscore)>1.5) %>% pull(ID.Young)), ] %>%
#   separate(full_boundary.Young, into=c("chr","start","end"), sep="_") %>% ungroup() %>% dplyr::select(chr,start,end) %>%
#   setNames(c("chr","start","end")) %>% distinct() %>% GRanges()
sig_young.domain = young.domain[young.domain$V4 %in% (bookend.bound_plt_df %>% dplyr::filter(abs(diff_zscore) > 1.5) %>% pull(ID.Young)), ] %>%
  dplyr::select(1:5) %>% setNames(c("chr","start","end","name","score")) %>% distinct() %>% GRanges()
sig_aged.domain = aged.domain[aged.domain$V4 %in% (bookend.bound_plt_df %>% dplyr::filter(abs(diff_zscore) > 1.5) %>% pull(ID.Aged)), ] %>%
  dplyr::select(1:5) %>% setNames(c("chr","start","end","name","score")) %>% distinct() %>% GRanges()


# Find gene promoters overlapping TAD domains ------------------------------------

library(GenomicRanges)
library(HelloRanges)

ab.switch = read.table(file.path("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\compartmentExpression\\compartmentBed\\100kb",
                                 "ab.switch.bed")) %>% setNames(c("chr","start","end","group")) %>% GRanges()

# Get all gene promoters within aged and young TAD domains
tss_slop_1kb = GRanges(read.table(file.path(projdir,"TAD expression","tss_1kb_slop.bed"), col.names = c("chr","start","end","gene")))
get_young_promoter_code = R_bedtools_intersect(a = tss_slop_1kb, 
                                               b = young.domain %>% dplyr::select(1:4) %>% setNames(c("chr","start","end","ID")) %>% GRanges(),
                                               wa = T, wb = T)
young_overlap_promoters = eval(get_young_promoter_code)
get_aged_promoter_code = R_bedtools_intersect(a = tss_slop_1kb, 
                                              b = aged.domain %>% dplyr::select(1:4) %>% setNames(c("chr","start","end","ID")) %>% GRanges(),
                                              wa = T, wb = T)
aged_overlap_promoters = eval(get_aged_promoter_code)

staticA_genes = eval(R_bedtools_intersect(tss_slop_1kb, ab.switch[ab.switch$group=="staticA", ]))

# need to import merged.classified.stats to get domain score information
merged.classified.stats = readRDS(file.path(projdir,"hicInterIntraTAD","merged.classified.stats.RDS"))

merged.classified.stats = merged.classified.stats %>% 
  group_by(Age) %>%
  mutate(bins_v2 = cut_interval(rank(domain_score)/n(), length=0.25))

table(merged.classified.stats$bins_v2, merged.classified.stats$Age)

ggplot(merged.classified.stats, aes(x=bins_v2, y=domain_score)) + geom_point(aes(color=Age))

# Add rearranged TAD type, domain score, and binned domain scores 
young_overlap_promoters@second$Type = bookend.bound_df[match(young_overlap_promoters@second$ID, bookend.bound_df$ID.Young), ]$Type
young_overlap_promoters@second$domain_score = merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$domain_score
young_overlap_promoters@second$bins = merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$bins
young_overlap_promoters@second$bins_v2 = merged.classified.stats[match(young_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$bins_v2

aged_overlap_promoters@second$Type = bookend.bound_df[match(aged_overlap_promoters@second$ID, bookend.bound_df$ID.Aged), ]$Type
aged_overlap_promoters@second$domain_score = merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$domain_score
aged_overlap_promoters@second$bins = merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$bins
aged_overlap_promoters@second$bins_v2 = merged.classified.stats[match(aged_overlap_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$bins_v2

young_overlap_promoters_type = cbind(gene = young_overlap_promoters@first$gene, 
                                     as.data.frame(young_overlap_promoters@second))
aged_overlap_promoters_type = cbind(gene = aged_overlap_promoters@first$gene, 
                                    as.data.frame(aged_overlap_promoters@second))

# pheatmap(drop_na(tpm.data[unique(overlap_promoters$gene), c(1:3,5:7)]), scale="row", cluster_cols = F, show_rownames = F, color = colorRampPalette(rev(brewer.pal(9,"RdBu")))(200))
avg.tpm.data = data.frame(A=rowMeans(tpm.data[,1:3]), Y=rowMeans(tpm.data[,5:7])) %>% tibble::rownames_to_column("gene")

overlap_promoters_df = rbind(left_join(young_overlap_promoters_type, avg.tpm.data[,c("gene","Y")], by="gene") %>% 
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

bin_promoters_df = inner_join(cbind(young_overlap_promoters_type, tpm.data[young_overlap_promoters_type$gene,5:7]) %>% 
                                drop_na() %>% dplyr::select(gene, bins, d0_Y_Rep3, d0_Y_Rep4, d0_Y_Rep5),
                              cbind(aged_overlap_promoters_type, tpm.data[aged_overlap_promoters_type$gene,1:3]) %>% 
                                drop_na() %>% dplyr::select(gene, bins, d0_A_Rep1, d0_A_Rep2, d0_A_Rep4),
                              by="gene", suffix=c("_Young","_Aged"))
bin_promoters_df %>% 
  dplyr::filter(bins_Young=="(0.8,1]" & bins_Aged=="(0.8,1]") %>%
  mutate(aged_avg = rowMeans(.[,c("d0_A_Rep1", "d0_A_Rep2", "d0_A_Rep4")]),
         young_avg = rowMeans(.[,c("d0_Y_Rep3", "d0_Y_Rep4", "d0_Y_Rep5")]),
         aged_sd = rowSds(as.matrix(.[,c("d0_A_Rep1", "d0_A_Rep2", "d0_A_Rep4")])),
         young_sd = rowSds(as.matrix(.[,c("d0_Y_Rep3", "d0_Y_Rep4", "d0_Y_Rep5")]))) %>%
  mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
  dplyr::filter(young_avg > aged_avg) %>%
  dplyr::filter(gene %in% DE.genes.name) %>%
  pull(gene) %>% unique() %>% sort() %>% writeClipboard()

static_A_s2n = data.frame(aged_avg = rowMeans(static.A.tpm[,1:3]),
                          young_avg = rowMeans(static.A.tpm[,5:7]),
                          aged_sd = rowSds(as.matrix(static.A.tpm[,1:3])),
                          young_sd = rowSds(as.matrix(static.A.tpm[,5:7])))
static_A_s2n %>% 
  mutate(s2n = (aged_avg - young_avg)/(aged_sd + young_sd)) %>%
  dplyr::select(s2n) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(gene %in% DE.genes.name) %>%
  write.table(file="clipboard-1000", sep='\t', col.names=F, row.names=F, quote=F)

## TAD classification expression --------------------

overlap_promoters_df %>%
  group_by(Age,Type) %>%
  tally()

shared_gene = overlap_promoters_df %>% filter(Type=="Shared") %>% pull(gene) %>% unique()
shifted_gene = overlap_promoters_df %>% filter(Type=="Shifted") %>% pull(gene) %>% unique()
merged_gene = overlap_promoters_df %>% filter(Type=="Merged") %>% pull(gene) %>% unique()
split_gene = overlap_promoters_df %>% filter(Type=="Split") %>% pull(gene) %>% unique()

overlap_promoters_df_test = overlap_promoters_df %>% group_by(Age) %>% wilcox_test(log2tpm ~ Type, ref.group="Shared") %>% add_xy_position(x="Age",group="Type",dodge=0.75)
ggplot(overlap_promoters_df, aes(x=Type,y=log2tpm)) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size = 0.1, width=0.65, position=position_dodge(0.75)) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  #theme_bw() + scale_fill_manual(values=hcl.colors(5, "Zissou 1")) +
  scale_y_continuous(breaks=seq(-10,15,2.5)) +
  # stat_pvalue_manual(overlap_promoters_df_test, tip.length=0.01, bracket.nudge.y = 0.5, size=4.5) +
  guides(fill = guide_legend(title="TAD\nType", nrow=2)) +
  labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        legend.margin = margin(t=0, r=0, b=0, l=-50, unit="pt"))
ggsave(file.path(projdir,"Figures","TAD_classification_expression.png"),dpi=300, width=3.5,height=4)  

## TAD classification differential expression -----------------------

DE.data = read.table("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\diff_d0_Y_vs_A.tsv",
                     sep="\t", header=T)
DE.data = tibble::rownames_to_column(DE.data, var="gene") %>%
  mutate(gene_symbol = sapply(gene, function(x) unlist(strsplit(x, split="_", fixed=T))[2])) %>%
  dplyr::filter(!grepl("Rik",gene) & !grepl("Gm",gene,fixed=F)) %>%
  dplyr::filter(abs(logFC)>log2(1.5)  & adj.P.Val<0.05)

overlap_promoters_df %>%
  dplyr::filter(Type=="Merged") %>%
  pull(gene) %>% unique() %>% writeClipboard()

overlap_promoters_df %>% 
  dplyr::filter(gene %in% DE.genes.name) %>%
  group_by(Age, Type) %>%
  tally() 

duplicate_genes = overlap_promoters_df %>% 
  dplyr::filter(gene %in% DE.genes.name) %>%
  dplyr::group_by(gene, Type, Age) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) %>%
  pull(gene) %>% unique()

overlap_promoters_df %>% 
  dplyr::filter(!(gene %in% duplicate_genes)) %>%
  dplyr::filter(Type=="Shifted") %>%
  pivot_wider(id_cols=c(gene,Type),
              names_from=Age, values_from=TPM) %>%
  mutate(log2FC = log2(Young/Aged)) %>%
  pull(gene) %>% writeClipboard()
  # dplyr::select(gene, log2FC) %>%
  # write.table(file='clipboard-1000', sep='\t', quote=F, row.names=F, col.names=F)

overlap_promoters_df %>% 
  dplyr::filter(!(gene %in% duplicate_genes)) %>%
  dplyr::filter(gene %in% DE.genes.name) %>%
  mutate(log2TPM = log2(TPM)) %>%
  ggviolin(x="Type",y="log2TPM",fill="Age",add="boxplot")
  
overlap_promoters_df %>% 
  dplyr::filter(!(gene %in% duplicate_genes)) %>%
  dplyr::filter(gene %in% DE.genes.name) %>%
  pivot_wider(id_cols=c(gene,Type),
              names_from=Age, values_from=TPM) %>%
  mutate(log2FC = log2(Young/Aged)) %>%
  ggviolin(x="Type", y="log2FC", add="jitter")

## log2(TPM) TAD type --------------------
overlap_promoters_df %>% group_by(Type_v2) %>% wilcox_test(log2tpm ~ Age)
wilcox.test(overlap_promoters_df[overlap_promoters_df$Type_v2=="Shared" & overlap_promoters_df$Age=="Young",]$log2tpm,
            overlap_promoters_df[overlap_promoters_df$Type_v2=="Unique" & overlap_promoters_df$Age=="Young",]$log2tpm)
wilcox.test(overlap_promoters_df[overlap_promoters_df$Type_v2=="Shared" & overlap_promoters_df$Age=="Aged",]$log2tpm,
            overlap_promoters_df[overlap_promoters_df$Type_v2=="Unique" & overlap_promoters_df$Age=="Aged",]$log2tpm)

ggplot(overlap_promoters_df, aes(x=Type_v2,y=log2tpm)) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size = 0.1, width=0.65, position=position_dodge(0.75)) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  # theme_bw() + scale_fill_manual(values=hcl.colors(5, "Zissou 1")) +
  scale_y_continuous(breaks=seq(-10,15,2.5)) +
  # stat_pvalue_manual(overlap_promoters_df_test, tip.length=0.01, bracket.nudge.y = 0.5, size=4.5) +
  guides(fill = guide_legend(title="TAD\nType", nrow=2)) +
  labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        legend.margin = margin(t=0, r=0, b=0, l=-50, unit="pt"))

## Differentially expressed genes in high Intra-TAD connectivity ----------------------

DE_genes_top = overlap_promoters_df %>%
  dplyr::filter(bins=="(0.8,1]" & (gene %in% DE.data$gene_symbol)) %>%
  dplyr::filter(Age=="Aged") %>%
  arrange(desc(domain_score)) 
  pull(gene) %>% unique() %>% sort() %>% writeClipboard()

DE_genes_overlap = overlap_promoters_df %>%
  dplyr::filter(gene %in% DE.data$gene_symbol)
DE_genes_overlap$logFC = DE.data$logFC[match(DE_genes_overlap$gene, DE.data$gene_symbol)]

DE_genes_overlap %>%
  group_by(bins, Age) %>%
  summarise(up_young = sum(sign(logFC)>0),
            up_aged = sum(sign(logFC)<0),
            up_young_frac = up_young/n(),
            up_aged_frac = up_aged/n()) %>%
  pivot_longer(cols=c(up_young_frac, up_aged_frac), names_to="direction", values_to="frac") %>%
  ggplot(aes(x=bins, y=frac)) +
  facet_wrap(~Age) +
  geom_col(aes(fill=direction), color="black", position="stack") +
  scale_fill_manual(values=pal_nejm()(2)) +
  theme_bw() + labs(x="Intra-TAD Connectivity Bins", y="Fraction of DE genes") +
  theme(axis.text.x = element_text(face="bold", color="black", size=12, angle=35, hjust=1),
        axis.text.y = element_text(face="bold", color="black", size=12),
        axis.title = element_text(face="bold", color="black", size=12),
        strip.text = element_text(face="bold", color="black", size=12),
        legend.position = "top")
ggsave(file.path(projdir,"Figures","DEgene_direction_TAD_connectivity.png"), dpi=300, width=4, height=4)

hist(sign(DE.data$logFC[DE.data$gene_symbol %in% DE_genes_top]))

# Ddit3, Cdk3n, Acly

DE_genes_up_in_young = overlap_promoters_df %>%
  dplyr::filter(bins=="(0.8,1]" & (gene %in% DE.data$gene_symbol[sign(DE.data$logFC)>0])) %>%
  arrange(desc(domain_score)) %>%
  pull(gene) %>% unique() %>% sort() 
writeClipboard(DE_genes_up_in_young)

overlap_promoters_df %>%
  dplyr::filter(bins=="[0,0.2]" & (gene %in% DE.data$gene_symbol)) %>%
  arrange(desc(domain_score)) %>%
  pull(gene) %>% unique() %>% sort() %>% writeClipboard()

## Domain score vs log2(TPM) Scatter plot --------------------
p = ggplot(overlap_promoters_df, aes(x=domain_score, y=log2tpm)) +
  geom_point(aes(color=Age), size=0.1, alpha=0.5) +
  stat_smooth(method="lm", aes(color=Age)) +
  #scale_x_continuous(breaks=seq(0.7,1,0.1)) +
  scale_color_manual(values=rev(pal_nejm()(2))) +
  stat_cor(aes(color=Age,
               label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  # stat_regline_equation(aes(color=Age,
  #                           label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"))) +
  labs(x="TAD Domain Score", y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme_pubr() + theme(legend.position="none", axis.title=element_text(face="bold"), axis.text=element_text(face="bold"))
png(file.path(projdir, "Figures", "domain_score_log2TPM_scatter.png"), res=300, units="in", width=4.5, height=4)
ggExtra::ggMarginal(p, groupFill = T, type="density")
dev.off()

## Domain score bins vs TPM percentile --------------------
domain_bin_TPM_test = overlap_promoters_df %>%
  group_by(Age) %>%
  # t_test(log2tpm ~ bins_v2, comparisons = list(c("[0,0.25]","(0.25,0.5]"),
  #                                              c("(0.25,0.5]","(0.5,0.75]"),
  #                                              c("(0.5,0.75]","(0.75,1]")))
  t_test(log2tpm ~ bins, comparisons = list(c("[0,0.2]","(0.2,0.4]"),
                                            c("(0.2,0.4]","(0.4,0.6]"),
                                            c("(0.4,0.6]","(0.6,0.8]"),
                                            c("(0.6,0.8]","(0.8,1]")))

overlap_promoters_df %>%
  group_by(Age, bins_v2) %>%
  summarise(num_genes = n_distinct(gene),
            mean_tpm = sd(log2tpm))

summary(lm(log2tpm ~ Age, data = overlap_promoters_df))

ggplot(overlap_promoters_df, aes(x=bins, y=log2tpm)) +
  stat_summary(fun="mean", aes(group=Age, color=Age), geom="path", position=position_dodge(0.1)) +
  stat_summary(fun.data="mean_se", aes(color=Age), geom="pointrange", position=position_dodge(0.1)) +
  scale_y_continuous(breaks=seq(3,4,0.1)) +
  theme_bw(base_size = 12) + scale_color_manual(values=rev(pal_nejm()(2))) +
  labs(x="Intra-TAD Connectivity Percentile Bins", y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  theme(axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.text.x = element_text(size=12, color="black", face="bold", angle=25, hjust=1, vjust=1),
        axis.title = element_text(size=12, face="bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "top")
ggsave(file.path(projdir,"Figures","Expression_percentile_vs_TAD_strength_point.png"),dpi=300, width=4,height=3.5)
  
ggplot(overlap_promoters_df, aes(x=bins, y=log2tpm)) +
  facet_grid(cols=vars(Age)) +
  geom_boxplot(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  # stat_pvalue_manual(overlap_promoters_percentile_test, tip.length=0.01, size=4, bracket.nudge.y = 0.025,
  #                    step.increase = 0.05, step.group.by = "Age") +
  scale_y_continuous(breaks=seq(0,1,0.2), labels=scales::percent_format()) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  labs(x="TAD Strength Percentile Bins", y="Gene Expression (TPM) Percentile") +
  theme(axis.text.x = element_text(size=12, color="black", face="bold", angle=30, hjust = 1, vjust=1),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "none")
ggsave(file.path(projdir,"Figures","Expression_percentile_vs_TAD_strength.png"),dpi=300, width=4.25,height=3.5)

overlap_promoters_df %>%
  dplyr::filter(gene %in% staticA_genes$gene) %>%
  ggplot(aes(x=Age, y=log2tpm)) +
  geom_boxplot(aes(fill=Type))


duplicate_genes = overlap_promoters_df %>%
  dplyr::filter(Type=="Split") %>%
  group_by(gene,Type,Age) %>%
  summarise(n=n(), .groups = "drop") %>%
  filter(n>1) %>% pull(gene)

overlap_promoters_df %>%
  #dplyr::filter(gene %in% staticA_genes$gene) %>%
  dplyr::filter(Type=="Split") %>%
  dplyr::filter(!(gene %in% duplicate_genes)) %>%
  dplyr::select(gene,Age,Type,TPM) %>%
  pivot_wider(names_from=Age, values_from=TPM) %>%
  mutate(log2FC = log2(Young/Aged)) %>%
  distinct() %>%
  dplyr::select(gene, log2FC) %>%
  arrange(log2FC) %>%
  #pull(gene) %>% writeClipboard()
  write.table(file="clipboard-1000", sep='\t', row.names=F, col.names=F, quote=F)
  #pull(gene) %>% unique() %>% writeClipboard()

overlap_promoters_df %>% dplyr::filter(gene=="Huwe1")
  
overlap_promoters_df %>% 
  dplyr::filter(gene %in% DE.genes.name, Type=="S") %>% 
  dplyr::select(gene, TPM, Age) %>% distinct() %>% 
  pull(gene) %>% writeClipboard()
  # pivot_wider(names_from=Age, values_from=TPM) %>% 
  # mutate(log2fc = log2(Aged/Young)) %>% 
  # dplyr::select(gene, log2fc) %>% 
  # write.table(file="clipboard-1000", sep="\t", row.names=F, col.names=F, quote=F)

# heatmap
drop_na(tpm.data[unique(overlap_promoters$gene), c(1:3,5:7)]) %>% 
  tibble::rownames_to_column("gene") %>%
  rowwise() %>%
  mutate(Aged = mean(c(d0_A_Rep1,d0_A_Rep2,d0_A_Rep4)),
         Young = mean(c(d0_Y_Rep3,d0_Y_Rep4,d0_Y_Rep5))) %>%
  ungroup() %>%
  dplyr::select(gene,Aged,Young) %>%
  mutate(foldchange = Aged/Young) %>%
  dplyr::select(gene,foldchange) %>%
  arrange(desc(foldchange)) %>%
  write.table(file="clipboard-1000",sep="\t",row.names = F,col.names = F,quote=F)

# Compare TAD boundary scores across age by class -------------------------

bookend_comparison = rbind(bookend.bound_df[,c("ID.Young","Score.Young","Type")] %>% setNames(c("ID","Score","Type")) %>% mutate(Age="Y") %>% drop_na() %>% distinct(),
                           bookend.bound_df[,c("ID.Aged","Score.Aged","Type")] %>% setNames(c("ID","Score","Type")) %>% mutate(Age="A") %>% drop_na() %>% distinct()) %>%
  mutate(Type = factor(Type, levels=rev(c("Combination", "Split", "Merged", "Shifted", "Shared"))),
         Age = factor(Age, levels=c("Y","A")))

bookend_comparison_test = bookend_comparison %>%
  group_by(Type) %>%
  wilcox_test(Score ~ Age) %>%
  add_xy_position(x="Type", group="Age", dodge=0.75)

bookend_comparison_test = bookend_comparison_test %>%
  mutate(p.adj.label = c("0.996","0.731","2.27e-3","3.20e-3","0.103"))

ggplot(bookend_comparison, aes(x=Type,y=Score)) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size=0.1, width=0.65, position=position_dodge(0.75)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.125)), breaks=seq(-2,2,0.5)) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  stat_pvalue_manual(bookend_comparison_test, label="p.adj.label", tip.length=0.01, size = 4.5) +
  labs(x="TAD Type", y="TAD Separation Score") +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir,"Figures","TAD_type_boundary_score.png"),dpi=300,width=5,height=3.5)

bookend.bound_df %>% drop_na() %>% distinct() %>% mutate(foldchange = log((Score.Aged+1)/(Score.Young+1))) %>% 
  ggplot(aes(x=Type,y=foldchange)) + geom_boxplot(aes(fill=Type))

# split violin plot
library(vioplot)
bookend_comparison %>%
  filter(Age=="Aged") %>%
  vioplot(Score~Type, data=., col = rev(bookend_type_dist$Color), lineCol=pal_nejm()(2)[1], plotCentre = "line", side = "right",
          ylim=c(-1.5,1.75), xlab="", ylab="TAD Separation Score",
          xaxt="n", yaxt="n", ann=F)
bookend_comparison %>%
  filter(Age=="Young") %>%
  vioplot(Score~Type, data=., lineCol = pal_nejm()(2)[2], plotCentre = "line", side = "left", add=T)
axis(side = 1, at = 1:5, labels = levels(bookend_comparison$Type), font = 2)
axis(side = 2, las = 2, font = 2, seq(-2,2,0.5))
mtext(side=2, text = "TAD separation score", font=2, line=2.75, cex=1.2)
legend("topright", fill = c(pal_nejm()(2)[c(2,1)], "gray"), legend = c("Young Unique", "Aged Unique","Shared"))


# export classified TAD domains as BED file with color (for pyG --------

young.domain.bed = young.domain[,1:5]
colnames(young.domain.bed) = c("chr","start","end","name","score")
young.domain.bed$Type = bookend.bound_df[match(young.domain.bed$name, bookend.bound_df$ID.Young), ]$Type
young.domain.bed$strand = "."
young.domain.bed$itemRgb = sapply(bookend_type_dist[bookend_type_dist$ID_age=="ID.Young", ]$Colors[match(young.domain.bed$Type, bookend_type_dist[bookend_type_dist$ID_age=="ID.Young", ]$Type)], 
                                  function(x) paste(col2rgb(x),collapse=","))
young.domain.bed = young.domain.bed[ ,c("chr", "start", "end", "name", "score", "strand", "start", "end", "itemRgb","Type")]
write.table(young.domain.bed, file.path(projdir, "young.TAD.domain.classified.bed"), sep='\t', row.names=F, col.names=F, quote=F)

aged.domain.bed = aged.domain[,1:5]
colnames(aged.domain.bed) = c("chr","start","end","name","score")
aged.domain.bed$Type = bookend.bound_df[match(aged.domain.bed$name, bookend.bound_df$ID.Aged), ]$Type
aged.domain.bed$strand = "."
aged.domain.bed$itemRgb = sapply(bookend_type_dist[bookend_type_dist$ID_age=="ID.Aged", ]$Colors[match(aged.domain.bed$Type, bookend_type_dist[bookend_type_dist$ID_age=="ID.Aged", ]$Type)], 
                                 function(x) paste(col2rgb(x),collapse=","))
aged.domain.bed = aged.domain.bed[ ,c("chr", "start", "end", "name", "score", "strand", "start", "end", "itemRgb","Type")]
write.table(aged.domain.bed, file.path(projdir, "aged.TAD.domain.classified.bed"), sep='\t', row.names=F, col.names=F, quote=F)


# Gene promoters overlapping classified TAD boundaries --------------------

young.boundary = read.table(file.path(projdir, "young.merged", "40kb", "young.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed"))
aged.boundary = read.table(file.path(projdir, "aged.merged", "40kb", "aged.merged_40kb_min120kb_max1200kb_step40kb_thresh0.01_delta_0.01_fdr_boundaries.bed"))

get_young_promoter_boundary_code = R_bedtools_intersect(a = tss_slop_1kb, 
                                               b = young.boundary %>% dplyr::select(1:5) %>% setNames(c("chr","start","end","ID","Score")) %>% GRanges(),
                                               wa = T, wb = T)
young_boundary_promoters = eval(get_young_promoter_boundary_code)
get_aged_promoter_boundary_code = R_bedtools_intersect(a = tss_slop_1kb, 
                                              b = aged.boundary %>% dplyr::select(1:5) %>% setNames(c("chr","start","end","ID","Score")) %>% GRanges(),
                                              wa = T, wb = T)
aged_boundary_promoters = eval(get_aged_promoter_boundary_code)

second(young_boundary_promoters)$boundary = unite(data=as.data.frame(second(young_boundary_promoters)), "boundary", 1:3, sep="_") %>% pull(boundary)
second(aged_boundary_promoters)$boundary = unite(data=as.data.frame(second(aged_boundary_promoters)), "boundary", 1:3, sep="_") %>% pull(boundary)

# need to import merged.classified.stats to get domain score information
merged.classified.stats = readRDS(file.path(projdir,"hicInterIntraTAD","merged.classified.stats.RDS"))

# Add rearranged TAD type, domain score, and binned domain scores 
young_boundary_promoters@second$Type = bookend.bound_df[match(young_boundary_promoters@second$boundary, bookend.bound_df$full_boundary.Young), ]$Type
young_boundary_promoters@second$ID = bookend.bound_df[match(young_boundary_promoters@second$boundary, bookend.bound_df$full_boundary.Young), ]$ID.Young
young_boundary_promoters@second$bins = merged.classified.stats[match(young_boundary_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="Y",]$name), ]$bins
aged_boundary_promoters@second$Type = bookend.bound_df[match(aged_boundary_promoters@second$boundary, bookend.bound_df$full_boundary.Aged), ]$Type
aged_boundary_promoters@second$ID = bookend.bound_df[match(aged_boundary_promoters@second$boundary, bookend.bound_df$full_boundary.Aged), ]$ID.Aged
aged_boundary_promoters@second$bins = merged.classified.stats[match(aged_boundary_promoters@second$ID, merged.classified.stats[merged.classified.stats$Age=="A",]$name), ]$bins

young_boundary_promoters_type = cbind(gene = young_boundary_promoters@first$gene, 
                                     as.data.frame(young_boundary_promoters@second))
aged_boundary_promoters_type = cbind(gene = aged_boundary_promoters@first$gene, 
                                    as.data.frame(aged_boundary_promoters@second))

boundary_promoters_df = rbind(left_join(young_boundary_promoters_type, avg.tpm.data[,c("gene","Y")], by="gene") %>% 
                               dplyr::rename(TPM="Y") %>% mutate(Age="Young") %>% drop_na(),
                             left_join(aged_boundary_promoters_type, avg.tpm.data[,c("gene","A")], by="gene") %>% 
                               dplyr::rename(TPM="A") %>% mutate(Age="Aged") %>% drop_na()) %>%
  mutate(Type = factor(Type, levels=c("Shared","Shifted","Merged","Split","Combination")),
         Type_v2 = if_else(Type=="Shared","Shared","Unique"),
         Age = factor(Age, levels=c("Young","Aged")),
         log2tpm = log2(TPM)) %>%
  group_by(Age) %>%
  mutate(tpm_percentile = rank(TPM)/length(TPM),
         score_percentile = rank(Score)/length(Score),
         score_bins = cut_width(rank(Score)/n(), width=0.2, center=0.5)) %>%
  ungroup()

# is expression correlated with boundary score?
boxplot(tpm_percentile ~ score_bins, data=boundary_promoters_df[boundary_promoters_df$Age=="Young",])

boundary_promoters_df_test = boundary_promoters_df %>% group_by(Type) %>% wilcox_test(log2tpm ~ Age) %>% add_xy_position(x="Type",group="Age",dodge=0.75)
# boundary_promoters_df_test = boundary_promoters_df %>% group_by(Age) %>% wilcox_test(log2tpm ~ Type) %>% add_xy_position(x="Age",group="Type",dodge=0.75)
boundary_promoters_df_test = boundary_promoters_df_test %>%
  mutate(p.label = c("0.0309", "9.23e-3", "0.0712", "0.0708", "4.47e-4"))

ggplot(boundary_promoters_df, aes(x=Type,y=log2tpm)) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size = 0.1, width=0.65, position=position_dodge(0.75)) +
  theme_bw() + scale_fill_manual(values=rev(pal_nejm()(2))) +
  #theme_bw() + scale_fill_manual(values=hcl.colors(5, "Zissou 1")) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-10,15,5)) +
  stat_pvalue_manual(boundary_promoters_df_test, label="p.label", tip.length=0.01, bracket.nudge.y = 0.5, size=4.5) +
  #guides(fill = guide_legend(title="TAD\nType", nrow=2)) +
  labs(x="TAD Type", y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  #labs(x=NULL, y="Gene Expression (TPM) Percentile") +
  theme(axis.text.x = element_text(size=12, color="black", face="bold"),
        axis.text.y = element_text(size=12, color="black", face="bold"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=11, face="bold"),
        legend.text = element_text(size=11),
        legend.margin = margin(t=0, r=0, b=0, l=-50, unit="pt"))
ggsave(file.path(projdir,"Figures","TAD_classification_boundary_expression.png"),dpi=300, width=5,height=3.5)  

# Look at lost and gained TAD boundaries by A/B compartment ----------------------------------

library(ggrepel)
library(latex2exp)

# lost (full) TAD boundaries in Young 
lost_TAD_boundaries = boundary.files.df %>% dplyr::filter(group=="Young") %>% dplyr::select(chr,start,end,Y) %>% distinct() %>% drop_na() %>% 
  mutate(Age="Young") %>% rename("Y"="Score") %>% GRanges()
# gained TAD boundaries in Aged
gained_TAD_boundaries = boundary.files.df %>% dplyr::filter(group=="Aged") %>% dplyr::select(chr,start,end,A) %>% distinct() %>% drop_na() %>% 
  mutate(Age = "Aged") %>% rename("A"="Score") %>% GRanges()
shared_TAD_boundaries = boundary.files.df %>% dplyr::filter(group=="Shared") %>% pivot_longer(cols=c(Y,A), names_to="Age", values_to="Score") %>%
  mutate(Age = sapply(Age, function(x) ifelse(x=="Y" & !is.na(x),"Young","Aged"))) %>% dplyr::select(chr,start,end,Score,Age) %>% 
  distinct() %>% drop_na() %>% arrange(Age) %>% GRanges()

# boundary.files.df %>% dplyr::filter(group=="Young") %>% dplyr::select(chr,start,end,Y) %>% drop_na() %>% distinct() %>%
#   write.table(file.path(projdir,""))

# Find overlap with A/B switching groups
lost_TAD_boundaries_ab_code = R_bedtools_intersect(lost_TAD_boundaries, ab.switch, wao=T)
lost_TAD_boundaries = eval(lost_TAD_boundaries_ab_code)
lost_TAD_boundaries = lost_TAD_boundaries[unlist(mcols(lost_TAD_boundaries)) > 1, ]

lost_TAD_boundaries_df = cbind(as.data.frame(first(lost_TAD_boundaries)),
                               as.data.frame(second(lost_TAD_boundaries))) %>%
  dplyr::select(Age,Score,group) %>%
  distinct() %>%
  group_by(group) %>%
  summarise(count = n())

gained_TAD_boundaries_ab_code = R_bedtools_intersect(gained_TAD_boundaries, ab.switch, wao=T)
gained_TAD_boundaries = eval(gained_TAD_boundaries_ab_code)
gained_TAD_boundaries = gained_TAD_boundaries[unlist(mcols(gained_TAD_boundaries)) > 1, ]

gained_TAD_boundaries_df = cbind(as.data.frame(first(gained_TAD_boundaries)),
                               as.data.frame(second(gained_TAD_boundaries))) %>% 
  dplyr::select(Age,Score,group) %>%
  distinct() %>%
  group_by(group) %>%
  summarise(count = n())

# make new data frame combining gained and lost boundaries to compare separation scores
gained_lost_boundaries = rbind(as.data.frame(first(lost_TAD_boundaries)) %>% 
                                 setNames(c("seqnames","start","end","width","strand","Score","Age")) %>% 
                                 mutate(group="Lost") %>% 
                                 distinct(),
                               as.data.frame(first(gained_TAD_boundaries)) %>% 
                                 setNames(c("seqnames","start","end","width","strand","Score","Age")) %>%
                                 mutate(group="Gained") %>% 
                                 distinct(),
                               as.data.frame(shared_TAD_boundaries) %>%
                                 setNames(c("seqnames","start","end","width","strand","Score","Age")) %>%
                                 mutate(group="Shared") %>% 
                                 distinct()) %>%
  mutate(group = factor(group, levels=c("Shared","Lost","Gained")),
         Age = factor(Age, levels=c("Young","Aged")))

gained_lost_boundaries_test = gained_lost_boundaries %>% group_by(Age) %>% wilcox_test(Score~group) %>% add_xy_position(x="Age", group="group", dodge=0.9)
wilcox.test(gained_lost_boundaries %>% dplyr::filter(Age=="Young", group=="Shared") %>% pull(Score),
            gained_lost_boundaries %>% dplyr::filter(Age=="Aged", group=="Shared") %>% pull(Score))
wilcox.test(gained_lost_boundaries %>% dplyr::filter(group=="Lost") %>% pull(Score),
            gained_lost_boundaries %>% dplyr::filter(group=="Gained") %>% pull(Score))
ggplot(gained_lost_boundaries, aes(x=Age, y=Score)) +
  geom_boxplot(aes(fill=group), color="black", outlier.size=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("gray",pal_nejm()(2))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1))) +
  theme_bw() + 
  guides(fill=guide_legend(title="Boundary", nrow=2)) +
  labs(x=NULL, y="TAD Separation Score") +
  #stat_pvalue_manual(gained_lost_boundaries_test, tip.length=0.01, label="p.adj", size=4) +
  geom_bracket(
    xmin = 1.225, xmax = 2.225, y.position = 1.9,
    label = "0.473", tip.length = 0.01, label.size=4,
  ) +
  geom_bracket(
    xmin = 0.775, xmax = 1.775, y.position = 2.25,
    label = "0.851", tip.length = 0.01, label.size=4,
  ) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face='bold'),
        legend.text = element_text(size=12),
        legend.margin = margin(t=0, r=0, b=0, l=-30))
ggsave(file.path(projdir,"Figures","TAD_separation_score_gained_lost.png"),dpi=300,width=2.5,height=3.25)

# Distance between gained/lost TAD boundaries to gene promoters
gene_anno <- rtracklayer::readGFF("C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\get_tss\\gencode.vM25.basic.annotation.gff3.gz")
gene_anno_subset = gene_anno[gene_anno$type=="gene" & gene_anno$gene_type=="protein_coding", c("seqid","start","end","gene_name")]
gene_anno_subset = gene_anno_subset[!grepl("^Gm", gene_anno_subset$gene_name) & !grepl("Rik", gene_anno_subset$gene_name), ]
gene_anno_subset = GRanges(gene_anno_subset)

quiescence_genes = read.table(file.path(projdir, "Quiescence genes.txt"))
diff_genes = read_xlsx(file.path(projdir,"GO_0042692_muscle_cell_differentiation.xlsx"))
genes_up_in_aged = avg.tpm.data[avg.tpm.data$A > avg.tpm.data$Y,]$gene
genes_down_in_aged = avg.tpm.data[avg.tpm.data$A < avg.tpm.data$Y,]$gene

get_distance_to_genes = function(gene_list, boundary_group) {
  tmp = eval({
    genome <- Seqinfo(genome = "mm10")
    gr_a <- GRanges(boundary.files.df[boundary.files.df$group==boundary_group,])
    gr_b <- resize(tss_slop_1kb, fix = "center", width = 1)[which(mcols(tss_slop_1kb)$gene %in% gene_list), ]
    # gr_b <- gene_anno_subset[unlist(mcols(gene_anno_subset) %in% gene_list), ]
    hits <- nearest(gr_a, gr_b, ignore.strand = TRUE, select = "all")
    ans <- pair(gr_a, gr_b, hits, all.x = TRUE)
    mcols(ans)$distance <- distance(ans)
    mcols(ans)$group = boundary_group
    ans
  })
  return(tmp)
}

closest_random_genes = function(gene_list, boundary_group) {
  tmp = eval({
    genome <- Seqinfo(genome = "mm10")
    gr_a <- sample(GRanges(boundary.files.df), replace = F, size = sum(boundary.files.df$group==boundary_group))
    gr_b <- gene_anno_subset[unlist(mcols(gene_anno_subset) %in% gene_list), ]
    hits <- nearest(gr_a, gr_b, ignore.strand = TRUE, select = "all")
    ans <- pair(gr_a, gr_b, hits, all.x = TRUE)
    mcols(ans)$distance <- distance(ans)
    mcols(ans)$group = boundary_group
    ans
  })
  return(tmp)
}

closest_dist_df = rbind(
  rbind(mcols(get_distance_to_genes(rownames(tpm.data), "Shared")),
        mcols(get_distance_to_genes(rownames(tpm.data), "Y")),
        mcols(get_distance_to_genes(rownames(tpm.data), "A"))) %>%
    as.data.frame() %>%
    mutate(gene_group = "Aged"),
  rbind(mcols(get_distance_to_genes(rownames(tpm.data), "Shared")),
        mcols(get_distance_to_genes(rownames(tpm.data), "Y")),
        mcols(get_distance_to_genes(rownames(tpm.data), "A"))) %>% 
    as.data.frame() %>%
    mutate(gene_group = "Young")
)

# closest_dist_df = rbind(mcols(get_distance_to_genes(unique(diff_genes$Symbol), "Shared")),
#                         mcols(get_distance_to_genes(unique(diff_genes$Symbol), "Young")),
#                         mcols(get_distance_to_genes(unique(diff_genes$Symbol), "Aged")))
# closest_dist_df = rbind(mcols(get_distance_to_genes(unique(avg.tpm.data$gene), "Shared")),
#                         mcols(get_distance_to_genes(unique(avg.tpm.data$gene), "Young")),
#                         mcols(get_distance_to_genes(unique(avg.tpm.data$gene), "Aged")))
# closest_random_dist_df = rbind(mcols(closest_random_genes(unique(avg.tpm.data$gene), "Shared")),
#                                mcols(closest_random_genes(unique(avg.tpm.data$gene), "Young")),
#                                mcols(closest_random_genes(unique(avg.tpm.data$gene), "Aged")))
merged_dist_df = as.data.frame(closest_dist_df)
merged_dist_df$random = as.data.frame(closest_random_dist_df)$distance

# statistical tests for bimodal data
ks.test(closest_dist_df$distance[closest_dist_df$group=="Aged"], closest_dist_df$distance[closest_dist_df$group=="Young"])
ks.test(closest_dist_df$distance[closest_dist_df$group=="Aged"], closest_dist_df$distance[closest_dist_df$group=="Shared"])
ks.test(closest_dist_df$distance[closest_dist_df$group=="Young"], closest_dist_df$distance[closest_dist_df$group=="Shared"])

as.data.frame(closest_dist_df) %>% wilcox_test(distance ~ group)
as.data.frame(closest_random_dist_df) %>% wilcox_test(distance ~ group)

closest_dist_test = as.data.frame(closest_dist_df) %>% 
  mutate(group = replace(group, group == "Aged", "Gained")) %>%
  mutate(group = replace(group, group == "Young", "Lost")) %>%
  wilcox_test(distance ~ group) %>% add_xy_position(x="group")
# 
# ggplot(as.data.frame(closest_dist_df) %>% 
#          mutate(group = replace(group, group == "Aged", "Gained")) %>%
#          mutate(group = replace(group, group == "Young", "Lost")), 
#        aes(x=group,y=log10(distance))) + 
#   geom_violin() +
#   geom_boxplot(width=0.1, outlier.shape=NA) +
#   labs(x="TAD Boundary Group", y="log10(distance to nearest gene promoter)") +
#   stat_pvalue_manual(closest_dist_test, tip.length=0.01, y.position = log10(closest_dist_test$y.position) + 0.5, step.increase = 0.1) +
#   theme_bw() + theme(axis.text = element_text(color="black", size=12), axis.title = element_text(color="black", size=12))
# ggsave(file.path(projdir, "Figures", "TAD_boundary_classification_gene_distance.png"), dpi=300, height=3.5, width=2.5)

## REPLACE WITH MARGINAL COMPARISON OF AT PROMOTER AND AWAY FROM PROMOTER

closest_dist_df %>%
  as.data.frame() %>%
  group_by(group) %>%
  summarise(at_promoter = sum(distance==0),
            not_at_promoter = sum(distance!=0)) %>%
  mutate(at_promoter_frac = at_promoter/(at_promoter + not_at_promoter),
         not_at_promoter_frac = not_at_promoter/(at_promoter + not_at_promoter))

closest_dist_df %>% group_by(gene_group, group) %>% summarise(away = sum(distance>0), at = sum(distance==0)) %>% mutate(away_frac = away/(away+at))
closest_dist_df %>% group_by(group) %>% summarise(away = sum(distance>0), at = sum(distance==0)) %>% mutate(away_frac = away/(away+at))

closest_dist_df %>% 
  group_by(group) %>% 
  summarise(away = sum(distance>3000), at = sum(distance<=3000)) %>% 
  mutate(away_frac = away/(away+at),
         at_frac = at/(away+at)) %>%
  pivot_longer(cols=c(away, at), names_to="frac_name", values_to="frac") %>%
  ggplot(aes(x=group, y=frac)) +
    geom_col(aes(fill=frac_name))

closest_dist_df_chisq = closest_dist_df %>% 
  group_by(group) %>% 
  summarise(away = sum(distance>1000), at = sum(distance<=1000)) %>% 
  as.data.frame() %>%
  tibble::column_to_rownames("group") %>%
  t()
fisher.test(closest_dist_df_chisq[,c("A","Shared")])
fisher.test(closest_dist_df_chisq[,c("Y","Shared")])

# split only by age
ggplot(as.data.frame(closest_dist_df) %>% 
         mutate(group = replace(group, group == "A", "Gained")) %>%
         mutate(group = replace(group, group == "Y", "Lost")), 
       aes(y=group, x=log10(distance+1))) + 
  geom_density_ridges(aes(fill=group), jittered_points = F, point_size=0.01) +
  scale_y_discrete(expand = expansion(mult=c(0,1))) +
  scale_fill_manual(values=c(pal_nejm()(2),"grey")) +
  geom_vline(xintercept = log10(1e3 + 1), lty=2) +
  theme_pubr() + labs(x=expression(bold(paste("lo","g"["10"],"(Distance to nearest TSS + 1)"))), y=NULL) + 
  theme(axis.text = element_text(face="bold", size=12),
        axis.title = element_text(face="bold", size=12),
        legend.position = "none")
ggsave(file.path(projdir, "Figures", "TAD_boundary_classification_gene_distance.png"), dpi=300, height=3, width=4)

# split by differentially expressed genes 
ggplot(as.data.frame(closest_dist_df) %>% 
         mutate(group = replace(group, group == "Aged", "Gained")) %>%
         mutate(group = replace(group, group == "Young", "Lost")), 
       aes(y=group, x=log10(distance+1))) + 
  facet_grid(gene_group~.) +
  geom_density_ridges(aes(fill=group), scale = 3, bandwidth=0.25, jittered_points = T, point_size=0.001, panel_scaling = F) +
  scale_y_discrete(expand = expansion(mult=c(0,0))) +
  #scale_y_discrete(expand = expansion(mult=c(0.1, 0.75))) +
  scale_fill_manual(values=c(pal_nejm()(2),"grey")) +
  geom_vline(xintercept = log10(1e3 + 1), lty=2) +
  theme_bw() + theme(axis.text = element_text(color="black", size=12), axis.title = element_text(color="black", size=12))
ggsave(file.path(projdir, "Figures", "TAD_boundary_classification_gene_distance_color.png"), dpi=300, height=5, width=5)

# compare gained/lost TAD boundaries with compartment switching group
as.data.frame(first(lost_TAD_boundaries))
TAD_ab_switch = inner_join(lost_TAD_boundaries_df %>% group_by(group) %>% summarise(num_boundaries = sum(count)),
                           gained_TAD_boundaries_df %>% group_by(group) %>% summarise(num_boundaries = sum(count)),
                           by="group", suffix=c(".lost",".gained"))

TAD_ab_switch_df = TAD_ab_switch %>% 
  pivot_longer(contains("num_boundaries"), names_to="boundary_type", values_to="number") %>%
  group_by(boundary_type) %>%
  mutate(fraction = number/sum(number)) %>%
  arrange(boundary_type, number) %>%
  mutate(group = factor(group, levels=c("AtoB","BtoA","staticB","staticA")))

ggplot(TAD_ab_switch_df, aes(x=boundary_type, y=number)) +
  geom_col(aes(fill=group), color="black", width=0.75, position=position_stack()) +
  geom_label_repel(aes(group=group, label = scales::percent(fraction,accuracy=0.1)), color="black", force=0.0005, 
            direction="y", size=3.5, position = position_stack(vjust=0.35), label.padding = unit(0.1, "lines")) +
  scale_x_discrete(labels=c("Gained", "Lost")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.1)), breaks=seq(0, 2000, 200)) +
  scale_fill_manual(values = rev(pal_jama()(6)[c(1:2,5:6)]),
                    labels = rev(c(TeX("Static A "),TeX("Static B "),TeX("$B \\rightarrow A$"),TeX("$A \\rightarrow B$")))) +
  theme_bw() + labs(x=NULL, y="# of TAD Boundaries") +
  guides(fill=guide_legend(nrow=2)) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.ticks.length.x = unit(0,"in"),
        axis.title = element_text(size=12, face="bold"),
        legend.justification = c(0,1),
        legend.margin = margin(t = 0, r = 0, b = 0, l = -25, unit = "pt"),
        legend.position = "top",
        legend.text = element_text(size=11, hjust=0),
        legend.title = element_blank())
ggsave(file.path(projdir, "Figures", "TAD_boundary_AB_compartment_switch.png"), dpi=300, height=3.5, width=2.25)

## get gene promoters at shared/lost/gained TAD boundaries -------------------- 
lost_promoters = eval(R_bedtools_intersect(tss_slop_1kb, first(lost_TAD_boundaries), wa=T, wb=T))
gained_promoters = eval(R_bedtools_intersect(tss_slop_1kb, first(gained_TAD_boundaries), wa=T, wb=T))
shared_promoters = eval(R_bedtools_intersect(tss_slop_1kb, shared_TAD_boundaries, wa=T, wb=T))

avg.tpm.data[avg.tpm.data$gene %in% first(shared_promoters)$gene, ] %>%
  drop_na() %>%
  distinct() %>%
  pivot_longer(c(A,Y), names_to="TPM.Age", values_to="TPM") %>%
  ggplot(aes(x=TPM.Age, y=log2(TPM))) +
  geom_boxplot() +
  stat_compare_means()

writeClipboard(unique(first(shared_promoters)$gene)[unique(first(shared_promoters)$gene) %in% rownames(tpm.data)])
writeClipboard(unique(first(lost_promoters)$gene)[unique(first(lost_promoters)$gene) %in% rownames(tpm.data)])
writeClipboard(unique(first(gained_promoters)$gene)[unique(first(gained_promoters)$gene) %in% rownames(tpm.data)])
writeClipboard(unique(first(shared_promoters)$gene))

avg.tpm.data[avg.tpm.data$gene %in% unique(first(gained_promoters)$gene), ] %>%
  drop_na() %>%
  mutate(l2fc = log2(A/Y)) %>%
  dplyr::select(gene,l2fc) %>%
  write.table(file="clipboard-1000", sep="\t", row.names = F, col.names = F, quote=F)

gained_lost_promoters = rbind(as.data.frame(first(lost_promoters)) %>% mutate(group="Lost") %>% distinct(),
                              as.data.frame(first(gained_promoters)) %>% mutate(group="Gained") %>% distinct(),
                              as.data.frame(first(shared_promoters)) %>% mutate(group="Shared") %>% distinct()) %>%
  mutate(group = factor(group, levels=c("Shared","Lost","Gained")))
gained_lost_promoters_df = inner_join(avg.tpm.data, gained_lost_promoters, by="gene")

# gained_lost_promoters_df_test = gained_lost_promoters_df %>%
#   drop_na() %>%
#   pivot_longer(c(A,Y), names_to="TPM.Age", values_to="TPM") %>%
#   mutate(TPM.Age = factor(TPM.Age, levels=c("Y","A")),
#          log2tpm = log2(TPM)) %>%
#   group_by(TPM.Age) %>%
#   wilcox_test(log2tpm ~ group, ref.group="Shared") %>%
#   add_xy_position(x="TPM.Age", group="group", dodge=0.9)
# gained_lost_promoters_df_test$p.adj.label = c("6.74e-5", "0.086", "5.38e-5", "0.118")
gained_lost_promoters_df_test = gained_lost_promoters_df %>%
  drop_na() %>%
  pivot_longer(c(A,Y), names_to="TPM.Age", values_to="TPM") %>%
  mutate(TPM.Age = factor(TPM.Age, levels=c("Y","A")),
         log2tpm = log2(TPM)) %>%
  group_by(TPM.Age) %>%
  wilcox_test(log2tpm ~ group, ref.group="Shared") %>%
  add_xy_position(x="TPM.Age", group="group", dodge=0.9)
# gained_lost_promoters_df_test$p.adj.label = c("6.74e-5", "0.086", "5.38e-5", "0.118")
gained_lost_promoters_df_test$p.adj.label = c("0.086", "6.74e-5", "0.118", "5.38e-5")
  
gained_lost_promoters_df %>%
  drop_na() %>%
  pivot_longer(c(A,Y), names_to="TPM.Age", values_to="TPM") %>%
  mutate(TPM.Age = factor(TPM.Age, levels=c("Y","A")),
         group_v2 = if_else(group=="Shared","Shared","Unique")) %>%
  ggplot(aes(x=TPM.Age, y=log2(TPM))) +
  geom_boxplot_pattern(aes(fill=TPM.Age, pattern=group), outlier.size=0.2, 
                       color="black",
                       pattern_density=0.1,
                       pattern_color="white",
                       pattern_fill="white",
                       position=position_dodge(0.9)) +
  scale_pattern_manual(values=c("none","stripe","circle")) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.1)), breaks=seq(-10,15,5)) +
  scale_x_discrete(labels=c("Young","Aged")) +
  theme_bw() + 
  guides(fill="none", pattern=guide_legend(title="Boundary", nrow=2, override.aes = list(fill="gray25"))) +
  labs(x=NULL, y=expression(bold(paste("lo","g"["2"],"(TPM)")))) +
  stat_pvalue_manual(gained_lost_promoters_df_test, tip.length=0.01, label="p.adj.label", size=4,
                   step.increase = 0.05, step.group.by = "TPM.Age") +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.position = "top",
        legend.title = element_text(size=12, face='bold'),
        legend.text = element_text(size=12),
        legend.margin = margin(t=0, r=0, b=0, l=-40))
ggsave(file.path(projdir,"Figures","TAD_boundary_gained_lost_expression.png"),dpi=300,width=2.75,height=3.25)

# heatmap gained/lost expression
tpm.data.tmp = tpm.data %>% tibble::rownames_to_column("gene")
gained_lost_promoters_mat = inner_join(tpm.data.tmp[,c(1:4,6:8)], gained_lost_promoters, by="gene")
gained_lost_promoters_mat = gained_lost_promoters_mat %>% 
  mutate(gene_idx = 1:n()) %>%
  mutate(row_mean = (d0_A_Rep1+d0_A_Rep2+d0_A_Rep4)/3) %>%
  group_by(group) %>%
  arrange(desc(row_mean), .by_group=T)
# +d0_Y_Rep2+d0_Y_Rep3+d0_Y_Rep4
gained_lost_mat = as.matrix(gained_lost_promoters_mat[,c(5:7,2:4)])
rownames(gained_lost_mat) = gained_lost_promoters_mat$gene_idx
# colnames(gained_lost_mat) = c("A","Y")

gained_lost_rowanno = data.frame(gained_lost_promoters_mat$group)
rownames(gained_lost_rowanno) = gained_lost_promoters_mat$gene_idx
colnames(gained_lost_rowanno) = "Group"

pheatmap(gained_lost_mat, annotation_row = gained_lost_rowanno, 
         width = 4, height=5, treeheight_col = 25, breaks = 3,
         annotation_colors = list(Group=c(Shared="gray", Lost=pal_nejm()(2)[2], Gained=pal_nejm()(2)[1])),
         color = colorRampPalette(colors = brewer.pal(7,"RdBu"))(200), 
         show_rownames = F, cluster_cols = T, cluster_rows=F, scale="row",
         filename = file.path(projdir, "Figures", "gained_lost_boundary_heatmap.png"))

# TAD class rGREAT --------------------------------------------------------

library(rGREAT)

bg = boundary.files.df %>% dplyr::select(chr,start,end) %>% distinct()
bg = bookend.bound_df %>% dplyr::select(full_boundary.Young) %>% 
  separate(full_boundary.Young, into=c("chr","start","end"), sep="_") %>% 
  drop_na() %>% distinct() %>% GRanges()

# Bound_List <- lapply(unique(bookend.bound_df$Type), function(x) {
#   bookend.bound_df %>% dplyr::filter(Type == x) %>% separate(boundary, into=c("chr","start","end"), sep="_") %>%
#     drop_na() %>% 
#     distinct() %>%
#     ungroup() %>%
#     mutate(start = as.integer(start), end = as.integer(end)) %>%
#     mutate(start = ifelse(boundary_name.Young=="left", start-20000, start),
#            end = ifelse(boundary_name.Young=="right", end+20000, end)) %>%
#     mutate(start = as.integer(start), end = as.integer(end)) %>%
#     dplyr::select(chr, start, end) %>%
#     dplyr::filter(chr %in% bg$chr & start %in% bg$start & end %in% bg$end)
# })

Bound_List <- lapply(unique(bookend.bound_df$Type), function(x) {
  bookend.bound_df %>% 
    dplyr::filter(Type == x) %>%
    dplyr::select(full_boundary.Young) %>% 
    separate(full_boundary.Young, into=c("chr","start","end"), sep="_") %>% 
    drop_na() %>% distinct() %>% GRanges()
})

test = submitGreatJob(gr = Bound_List[[2]], bg = bg, species="mm10", rule="basalPlusExt", adv_upstream=5, adv_downstream=1)

Bound_GREAT <- lapply(Bound_List, function(x) {
  getEnrichmentTables(submitGreatJob(x, bg = bg, species = "mm10", rule="basalPlusExt", adv_upstream = 5, adv_downstream = 1))
})

bind_rows(Bound_GREAT[[6]]) %>% filter(Hyper_Fold_Enrichment>2, Hyper_Adjp_BH<0.05)

shifted_GREAT = bookend.bound_df %>% 
  filter(Type=="Shifted") %>% 
  select(boundary) %>% 
  separate(boundary, into=c("chr","start","end"), sep="_") %>% 
  drop_na() %>% 
  distinct() %>%
  ungroup() %>%
  select(-ID.Young) %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  submitGreatJob(., species = "mm10", rule="basalPlusExt", adv_upstream = 5, adv_downstream = 1)
shifted_table = getEnrichmentTables(shifted_GREAT)
bind_rows(shifted_table)

# Plot examples of TAD classes --------------------------------------------

# library(Sushi)
# library(HiCcompare)
library(plotgardener)
library(RColorBrewer)
# library(HiCBricks)
mat_dir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC"
aged_mat_files = list.files(file.path(mat_dir,"aged.merged","DumpedMatrices","40kb"), full.names=T)
young_mat_files = list.files(file.path(mat_dir,"young.merged","DumpedMatrices","40kb","allChr"), full.names=T)
young_mat = read.table(grep("_chr3_",young_mat_files,value=T,fixed=T), sep="\t")
aged_mat = read.table(grep("_chr3_",aged_mat_files,value=T,fixed=T), sep="\t")
# young_mat_full = HiCcompare::sparse2full(young_mat, hic.table=F, column.name=NA)
# aged_mat_full = HiCcompare::sparse2full(aged_mat, hic.table=F, column.name=NA)

young_mat_log1p = young_mat %>% mutate(V3 = log(V3+1))
aged_mat_log1p = aged_mat %>% mutate(V3 = log(V3+1))
#chrom = "chr10"; chromstart = 108e6; chromend = 120e6
# chr3:126,572,095-138,694,044
chrom = "chr3"; chromstart = 126572095; chromend = 138694044

young.domain.triangles = young.domain %>% dplyr::filter(V1==chrom, V2>=chromstart, V3<=chromend)
young.domain.triangles$V9 = young.domain.bed$itemRgb[match(young.domain.triangles$V4, young.domain.bed$name)]
aged.domain.triangles = aged.domain %>% dplyr::filter(V1==chrom, V2>=chromstart, V3<=chromend)

bookend.bound_df %>% 
  separate(boundary,into=c("chr","start","end"),sep="_",remove=F) %>% 
  mutate(start=as.numeric(start), end=as.numeric(end)) %>% 
  dplyr::filter(chr==chrom,start>=chromstart,end<=chromend) %>% 
  dplyr::select(c(contains("Young"),"Type",-"boundary_name.Young")) %>% 
  distinct() %>%
  as.data.frame()

## Define genomic region with `params`
params <- pgParams(chrom = chrom, chromstart = chromstart, chromend = chromend,
                   x = 0.5, width = 3, assembly = "mm10")

png(file.path(projdir,"Figures","plotgardener_TAD_classification_figure.png"),res=300,width=4,height=3,units="in")
## Create a plotgardener page
pageCreate(width = 4, height = 3, default.units = "inches")

## Plot Hi-C data
young_hicPlot <- plotHicTriangle(
  data = young_mat_log1p,
  resolution = 40e3,
  #half="top",
  palette = colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  #zrange = c(0,2),
  zrange = c(0,4),
  bg = NA,
  params = params,
  y = 0.5, 
  height = 0.75,
  just = c("left", "top"), default.units = "inches"
)

young_domainsTriangle = annoDomains(plot=young_hicPlot, data=young.domain.triangles, params=params)
young_domainsPlot = plotRanges(young.domain.bed, linecolor="black",
                               params = params, collapse=T, 
                               fill = sapply(strsplit(young.domain.bed$itemRgb, ","), 
                                             function(x) rgb(x[1], x[2], x[3], maxColorValue=255)),
                               y="0.01b", width=3, height=0.1, just = c("left", "top"),
                               default.units = "inches")
aged_domainsPlot = plotRanges(aged.domain, chrom, linecolor="black",
                              params = params, collapse=T, fill = pal_nejm()(2)[1],
                              y="0b", width=3, height=0.1, just = c("left", "top"),
                              default.units = "inches")

aged_hicPlot <- plotHicTriangle(
  data = aged_mat_log1p,
  resolution = 40e3,
  palette = colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))),
  #zrange = c(0,2.5),
  zrange = c(0,4),
  bg = NA,
  #half="bottom",
  flip = T,
  params = params,
  y = "0.01b", 
  height = 0.75,
  just = c("left", "top"), default.units = "inches"
)

aged_domainsTriangle = annoDomains(plot=aged_hicPlot, data=aged.domain.triangles, params=params) 

## Annotate x-axis genome label
annoGenomeLabel(
  plot = young_hicPlot, scale = "Mb", params=params, y = "0b",
  just = c("left", "top")
)

## Annotate heatmap legend
annoHeatmapLegend(
  plot = young_hicPlot, x = 3.5, y = 0.5,
  width = 0.13, height = 0.6, 
  just = c("right", "top"), fontcolor = "black"
)

pageGuideHide()
dev.off()

#par(mfrow=c(2,1),mar=c(2,1,2,1))
# layout(matrix(c(1,2,3,4), 4, 1, byrow = TRUE), widths=3, heights=c(4,1,1,4))
# phic = plotHic(aged_mat_full, chrom, chromstart, chromend, max_y = 40, zrange=c(0,2), palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))))
# addlegend(phic[[1]], palette=phic[[2]], title="score", side="right",
#           bottominset=0.4, topinset=0, xoffset=-.035, labelside="left",
#           width=0.025, title.offset=0.035)
# plotBed(aged.domain,chrom,chromstart,chromend,height=0.1)
# plotBed(young.domain,chrom,chromstart,chromend,height=0.1)
# phic2 = plotHic(young_mat_full, chrom, chromstart, chromend, max_y = 40, zrange=c(0,2), palette=colorRampPalette(rev(brewer.pal(n = 9, "RdBu"))), flip=T)
# labelgenome(chrom, chromstart, chromend, n=4, scale="Mb", edgeblankfraction=0.20)
# dev.off()

# Use TADCompare to classify domains --------------------------------------

library(TADCompare)

# Running TADCompare with pre-specified TADs
# Import Hi-C matrices in 3 column format per chromosome (start, end, count)
# TADCompare converts to nxn matrix internally
mat_dir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\02_HIC"
aged_mat_files = list.files(file.path(mat_dir,"aged.merged","DumpedMatrices","40kb"), full.names=T)
young_mat_files = list.files(file.path(mat_dir,"young.merged","DumpedMatrices","40kb","allChr"), full.names=T)

young_mat = read.table(young_mat_files[1], sep="\t")
aged_mat = read.table(aged_mat_files[1], sep="\t")

young_domain_coords = young.domain[,1:3] %>% setNames(c("chr","start","end")) %>% filter(chr=="chr1") 
aged_domain_coords = aged.domain[,1:3] %>% setNames(c("chr","start","end")) %>% filter(chr=="chr1") 
# Combining the TAD boundaries for both contact matrices
Combined_Bed = list(young_domain_coords, aged_domain_coords)

TD_Compare = TADCompare(cont_mat1 = young_mat, cont_mat2 = aged_mat, resolution = 40000)
TD_Compare_pre = TADCompare(cont_mat1 = young_mat, cont_mat2 = aged_mat, resolution = 40000, pre_tads = Combined_Bed)
table(TD_Compare$TAD_Frame$Boundary %in% bind_rows(Combined_Bed)$end)
table(TD_Compare$Boundary_Scores$Type)
table(TD_Compare_pre$TAD_Frame$Type)
names(TD_Compare)

TD_Compare_pre$Count_Plot
TD_Compare$Boundary_Scores[TD_Compare$Boundary_Scores$Type=="Shifted",]

# Returning the boundaries
head(TD_Compare$TAD_Frame)
TD_Compare$Boundary_Scores[TD_Compare$Boundary_Scores$Type == "Shifted", ]

bed_coords1 = bind_rows(SpectralTAD(young_mat, chr = "chr1", levels = 3))
bed_coords2 = bind_rows(SpectralTAD(aged_mat, chr = "chr1", levels = 3))

head(TD_Compare$TAD_Frame)
head(bed_coords2)

# Visualize TADCompare results --------------------------------------------

young_mat_total <- sum(young_mat$V3)
aged_mat_total <- sum(aged_mat$V3)
(scaling_factor <- young_mat_total / aged_mat_total)
# Rescale matrices depending on which matrix is smaller
if (young_mat_total > aged_mat_total) {
  aged_mat$V3 <- aged_mat$V3 * scaling_factor
} else {
  young_mat$V3 <- young_mat$V3 * (1 / scaling_factor)
}
# Coordinates of interesting regions
start_coord <- 13000000
end_coord   <- 17000000
# Another interesting region
# start_coord <- 3000000 
# end_coord   <- 10000000
p <- DiffPlot(tad_diff    = TD_Compare_pre, 
              cont_mat1   = young_mat,
              cont_mat2   = aged_mat,
              resolution  = 40000,
              start_coord = start_coord,
              end_coord   = end_coord,
              pre_tad     = Combined_Bed,
              show_types  = TRUE, 
              point_size  = 2.5,
              max_height  = 10,
              rel_heights = c(1.5, 2),
              palette     = "RdYlBu")
p

cont_mat1 = HiCcompare::sparse2full(young_mat)
cont_mat2 = HiCcompare::sparse2full(aged_mat)

# TAD boundary scatterplot ------------------------------------------------

# scatter.df=cbind(left_join(boundary.files.df,young.scores,by=c("chr","start","end")),
#                  left_join(boundary.files.df,aged.scores,by=c("chr","start","end")))
scatter.df=left_join(left_join(boundary.files.df,young.scores,by=c("chr","start","end")), aged.scores,by=c("chr","start","end"),suffix=c(".young",".aged"))

scatter.df %>%
  distinct() %>%
  ggplot(aes(x=score.aged,y=score.young)) +
  geom_abline(linetype="dashed") +
  geom_point(aes(color=group), alpha=0.5, stroke=0.1) +
  scale_color_manual(values=c(pal_nejm()(2)[c(2,1)], "gray50")) +
  theme_bw() + 
  coord_fixed(ratio = 1, xlim = c(-1.75,1.75), ylim = c(-1.75,1.75)) +
  #scale_x_continuous(breaks=seq(-1.5,1.5,0.5)) +
  #scale_y_continuous(breaks=seq(-1.5,1.5,0.5)) +
  xlab("TAD boundary strength in Aged") +
  ylab("TAD boundary strength in Young") +
  guides(colour = guide_legend(title="Boundary",override.aes = list(size=3))) +
  theme(axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.position = "top")
ggsave(file.path(projdir, "Figures", "TAD_insulation_scores_40kb_0.01FDR_scatterplot.png"), width=4, height=4, dpi=300)


# Heatmap without clustering ----------------------------------------------

young.mat = boundary.files.df %>% select(Y,group) %>% drop_na()
aged.mat = boundary.files.df %>% select(A,group) %>% drop_na()

young.plt = Heatmap(as.matrix(young.mat[,1]), split = young.mat$group)
draw(young.plt)
aged.plt = Heatmap(as.matrix(aged.mat[,1]), split = aged.mat$group)
draw(aged.plt)

# Heatmap and clustering of TAD separation scores -------------------------

library(pheatmap)
library(circlize)

hmp.df = read.table(file.path(projdir,"union_young_aged_40kb_TAD_score.txt"), col.names=c("chr","start","end","Y","A"))
hmp.df = hmp.df[apply(hmp.df,1,function(x) !all(is.na(x[4:5]))),]
hmp.df = hmp.df %>% mutate(Y=as.numeric(Y), A=as.numeric(A)) 
hmp.mat = as.matrix(hmp.df[,4:5])
hmp.mat = apply(hmp.mat,2,as.numeric)

col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
young.hmp.plt = Heatmap(hmp.mat[!is.na(hmp.mat)[,1],1], cluster_columns = F, km=4, column_names_rot = 90,
                  na_col="black",
                  heatmap_legend_param = list(direction = "horizontal",
                                              title = "TAD Separation Score"))
draw(young.hmp.plt)
aged.hmp.plt = Heatmap(hmp.mat[!is.na(hmp.mat)[,2],2], cluster_columns = F, km=4 , column_names_rot = 90,
                        na_col="black",
                        heatmap_legend_param = list(direction = "horizontal",
                                                    title = "TAD Separation Score"))
young.hmp.plt+aged.hmp.plt
hmp.plt = draw(young.hmp.plt+aged.hmp.plt, heatmap_legend_side = "bottom")
hmp.dend = row_dend(young.hmp.plt)
hmp.cluster = melt(row_order(hmp.plt))

png(file.path(projdir,"Figures","TAD_separation_score_kmeans_hmp.png"), res=300, units="in", width=2, height=5)
hmp.plt
dev.off()

hmp.bxplt = data.frame(hmp.mat)
hmp.bxplt$cluster = hmp.cluster[order(hmp.cluster$value),"L1"]

hmp.bxplt %>%
  pivot_longer(cols=c("Y","A"), names_to="Age", values_to="score") %>%
  mutate(Age=factor(Age, levels=c("Y","A"))) %>%
  ggplot(aes(x=Age,y=score)) +
  geom_boxplot(aes(fill=Age), color="black", outlier.size=0.1) +
  facet_grid(rows=vars(cluster)) +
  theme_bw() + 
  labs(x="",y="") +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  theme(axis.text.x=element_text(size=14, face="bold", color="black"),
        legend.position="none")
ggsave(file.path(projdir,"Figures","TAD_separation_score_kmeans.png"),dpi=300,width=2,height=5)

# write TAD boundary bed files per cluster
for (i in unique(hmp.cluster$L1)) {
  hmp.df[hmp.cluster[hmp.cluster$L1==i,"value"],1:4] %>% 
    arrange(V1,V2) %>%
    filter(V4!=0) %>%
    write.table(file.path(projdir,"TAD expression","kmeans_boundaries",paste0("young_cluster",i,".bedgraph")), sep="\t", quote=F, col.names=F, row.names=F)
  hmp.df[hmp.cluster[hmp.cluster$L1==i,"value"],c(1:3,5)] %>% 
    arrange(V1,V2) %>%
    filter(V5!=0) %>%
    write.table(file.path(projdir,"TAD expression","kmeans_boundaries",paste0("aged_cluster",i,".bedgraph")), sep="\t", quote=F, col.names=F, row.names=F)
}

# Plot venn diagram of TAD boundaries -------------------------------------


boundary.sets = list(which(!is.na(boundary.files.df$Y)),
                     which(!is.na(boundary.files.df$A)))
boundary.sets.list = list(intersection=intersect(boundary.sets[[1]], boundary.sets[[2]]),
                         unique.y=setdiff(boundary.sets[[1]], boundary.sets[[2]]),
                         unique.a=setdiff(boundary.sets[[2]], boundary.sets[[1]]))
venn.data = melt(data.frame(lapply(boundary.sets.list, function(x) length(x))))
venn.data$fraction = venn.data$value/sum(venn.data$value)
venn.data$variable = factor(venn.data$variable, levels=c("intersection","unique.y","unique.a"))

ggplot(venn.data, aes(x="", y=fraction, fill=variable)) +
  geom_col(width=1) +
  geom_text(aes(label = sprintf("%d\n(%s)", value, percent(fraction))), size=3, position = position_stack(vjust = 0.5)) +
  #facet_wrap(~Age) +
  coord_polar("y", start=0) +
  scale_y_continuous(labels = scales::percent) +
  xlab("") + ylab("") +
  theme_pubclean() + 
  scale_fill_brewer(palette="Blues", labels=c("Common", "Unique Young", "Unique Aged")) +
  theme(axis.text = element_text(size=8, color="black"),
        axis.ticks.length = unit(0,"in"),
        legend.text = element_text(size=9),
        legend.title = element_blank())
ggsave(file.path(projdir,"Figures","unique_common_boundaries_piechart_40kb_0.01FDR.png"), dpi=300, width=4, height=4)


# Plot TAD separation score at promoter / non-promoter boundaries ---------

aged.bound.promoter = read.table(file.path(projdir,"TAD expression","promoterBoundary","aged.merged.TADBoundary.promoter.bed"), sep="\t")
aged.bound.nonpromoter = read.table(file.path(projdir,"TAD expression","promoterBoundary","aged.merged.TADBoundary.nonpromoter.bed"), sep="\t")
young.bound.promoter = read.table(file.path(projdir,"TAD expression","promoterBoundary","young.merged.TADBoundary.promoter.bed"), sep="\t")
young.bound.nonpromoter = read.table(file.path(projdir,"TAD expression","promoterBoundary","young.merged.TADBoundary.nonpromoter.bed"), sep="\t")
promoter.bound.df = rbind(aged.bound.promoter, aged.bound.nonpromoter, young.bound.promoter, young.bound.nonpromoter)
colnames(promoter.bound.df) = c("chr","start","end","name","score","strand")
promoter.bound.df$Age = factor(c(rep("A", each=nrow(aged.bound.promoter)+nrow(aged.bound.nonpromoter)),
                                 rep("Y", each=nrow(young.bound.promoter)+nrow(young.bound.nonpromoter))),
                               levels=c("Y","A"))
promoter.bound.df$Group = c(rep("Promoter", nrow(aged.bound.promoter)),
                            rep("Non-Promoter", nrow(aged.bound.nonpromoter)),
                            rep("Promoter", nrow(young.bound.promoter)),
                            rep("Non-Promoter", nrow(young.bound.nonpromoter)))

promoter.bound.df.test = promoter.bound.df %>%
  group_by(Group) %>%
  t_test(score~Age) %>%
  add_xy_position(x="Group", dodge=0.9)

promoter.bound.df %>%
  group_by(Group,Age) %>%
  summarise(mean=mean(score))

tmp = subset(promoter.bound.df, Group=="Promoter")
t.test(tmp$score~tmp$Age)
tmp = subset(promoter.bound.df, Group=="Non-Promoter")
t.test(tmp$score~tmp$Age)
tmp = subset(promoter.bound.df, Age=="A")
t.test(tmp$score~tmp$Group)
tmp = subset(promoter.bound.df, Age=="Y")
t.test(tmp$score~tmp$Group)

ggplot(promoter.bound.df, aes(x=Group, y=score)) +
  geom_violin(aes(fill=Age), color="black", position = position_dodge(0.9)) + 
  geom_boxplot(aes(group=interaction(Group,Age)), width=0.2, color="black", position = position_dodge(0.9), outlier.shape=NA) +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  scale_y_continuous(breaks=seq(-1.5,2.5,0.5), expand=expansion(mult=c(0.05,0.1))) +
  #stat_pvalue_manual(promoter.bound.df.test, tip.length = 0.01, bracket.nudge.y = 0.1) +
  geom_bracket(
    xmin = 0.775, xmax = 1.225, y.position = 1.6,
    label = "0.569", tip.length = 0.01, label.size = 5,
  ) +
  geom_bracket(
    xmin = 1.775, xmax = 2.225, y.position = 1.5,
    label = "0.986", tip.length = 0.01, label.size=5,
  ) +
  geom_bracket(
    xmin = 0.775, xmax = 1.775, y.position = 2,
    label = "0.955", tip.length = 0.01, label.size = 5,
  ) +
  geom_bracket(
    xmin = 1.225, xmax = 2.225, y.position = 2.4,
    label = "0.521", tip.length = 0.01, label.size=5,
  ) +
  theme_bw() + 
  xlab("") +  ylab("TAD Separation Score") +
  theme(legend.position="top",
        axis.text = element_text(size=14, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,"Figures","TAD_separation_scores_40kb_0.01FDR_promoter.png"), width=4, height=5, dpi=300)

# ChIPseeker --------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ChIPseeker)

young.bed = readPeakFile(unlist(subset(bound.file.df, Resolution=="40kb" & Age=="Y" & FDR==0.01, "files")))
aged.bed = readPeakFile(unlist(subset(bound.file.df, Resolution=="40kb" & Age=="A" & FDR==0.01, "files")))

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
young.tagMatrix = getTagMatrix(young.bed, windows=promoter)
aged.tagMatrix = getTagMatrix(young.bed, windows=promoter)

tagHeatmap(aged.tagMatrix, xlim=c(-5000, 5000), color="red")

# plotAvgProf(young.tagMatrix, xlim=c(-5000, 5000),
#             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency",
#             conf = 0.95, resample = 1000)

genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")
young.matrix_extension <- getTagMatrix(young.bed,windows = genebody, nbin = 800,
                                             upstream = 3000, downstream = 3000)
aged.matrix_extension <- getTagMatrix(aged.bed,windows = genebody, nbin = 800,
                                             upstream = 3000, downstream = 3000)
profile_plot = plotPeakProf(list(Young=young.matrix_extension,
                                 Aged=aged.matrix_extension),
                            conf = 0.95, ylab="TAD boundary count frequency") 
profile_plot +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  scale_color_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme(legend.position = "top",
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=14, face="bold"),
        legend.text = element_text(size=14))
ggsave(file.path(projdir,"Figures","TAD_boundary_frequency_genebodies_v2.png"), width=6, height=4, dpi=300)

young.peakAnno <- annotatePeak(young.bed, tssRegion=c(-1000, 1000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")
aged.peakAnno <- annotatePeak(aged.bed, tssRegion=c(-1000, 1000),
                               TxDb=txdb, annoDb="org.Mm.eg.db")

young.peakAnno.df = data.frame(young.peakAnno@annoStat)
aged.peakAnno.df = data.frame(aged.peakAnno@annoStat)
young.peakAnno.df$Count = 0.01*young.peakAnno.df$Frequency * young.peakAnno@peakNum
aged.peakAnno.df$Count = 0.01*aged.peakAnno.df$Frequency * aged.peakAnno@peakNum
young.peakAnno.df$Age = "Y"
aged.peakAnno.df$Age = "A"
peakAnno.df = rbind(young.peakAnno.df, aged.peakAnno.df)
peakAnno.df$Feature = factor(peakAnno.df$Feature)
peakAnno.df$Feature = factor(peakAnno.df$Feature, levels=rev(levels(peakAnno.df$Feature)))

ggplot(peakAnno.df, aes(y=Age, x=Count)) +
  geom_col(aes(fill=Feature)) +
  scale_x_continuous(expand=expansion(mult=c(0,0)), breaks=seq(0,3000,500)) +
  scale_y_discrete(expand=expansion(mult=c(0,0.1))) +
  scale_fill_manual(values=rev(brewer.pal(length(unique(peakAnno.df$Feature)),"Paired"))) +
  theme_bw() + 
  #xlab("Frequency (%)") +
  xlab("# TAD Boundaries") +
  guides(fill = guide_legend(reverse=TRUE)) +
  theme(axis.text=element_text(size=14, face="bold", color="black"),
        axis.title=element_text(size=14, face="bold", color="black"),
        legend.title=element_text(face="bold"),
        legend.text = element_text(size=12))
ggsave(file.path(projdir, "Figures", "TAD_boundary_count_annotations.png"), dpi=300, width=6, height=3)

png(file.path(projdir,"Figures","young_TAD_boundary_anno_pie.png"), width=6, height=4, res=300, units="in")
plotAnnoPie(young.peakAnno)
dev.off()
png(file.path(projdir,"Figures","aged_TAD_boundary_anno_pie.png"), width=6, height=4, res=300, units="in")
plotAnnoPie(aged.peakAnno)
dev.off()
ChIPseeker::upsetplot(young.peakAnno, vennpie=T)

plotDistToTSS(young.peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")
plotDistToTSS(aged.peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

# annotate boundaries ------------------------------------------------------------------

library(annotatr)
library(circlize)
#bed = circlize::generateRandomBed(nr = 1000, nc = 0)
#str(bed)
# bed[1:2, ]
young.bed = read.table(unlist(subset(bound.file.df, Resolution=="40kb" & Age=="Y" & FDR==0.01, "files")),
                sep = '\t', header=F)
aged.bed = read.table(unlist(subset(bound.file.df, Resolution=="40kb" & Age=="A" & FDR==0.01, "files")),
                       sep = '\t', header=F)
young.bed = young.bed[complete.cases(young.bed),1:4] # isolate BED components
young.bed$Age = "Y"
aged.bed = aged.bed[complete.cases(aged.bed),1:4] # isolate BED components
aged.bed$Age = "A"
bed = rbind(young.bed, aged.bed)
colnames(bed) = c("chr","start","end","bound.name","Age")
# write merged table for read_regions
write.table(bed, file.path(projdir, "TAD expression", "young.aged_40kb_0.01FDR.bed"), sep="\t", quote=F, row.names=F, col.names=F)

dm_regions = read_regions(con = file.path(projdir, "TAD expression", "young.aged_40kb_0.01FDR.bed"), 
                          extraCols = c(bound.name="character", Age="character"), 
                          genome = 'mm10',
                          format = 'bed')

grep("mm10", builtin_annotations(), value=T) # get available annotations for mm10
annots = c("mm10_genes_1to5kb","mm10_genes_promoters","mm10_genes_cds",
           "mm10_genes_5UTRs","mm10_genes_exons","mm10_genes_firstexons","mm10_genes_intergenic",
           "mm10_genes_introns","mm10_enhancers_fantom", "mm10_basicgenes")
#annotations = build_annotations(genome = 'mm10', annotations = annots) # Build the annotations (a single GRanges object)
#saveRDS(annotations, file.path(projdir, "annotatr_mm10_annotations.RDS"))
annotations = readRDS(file.path(projdir, "annotatr_mm10_annotations.RDS"))
# Intersect the regions we read in with the annotations
dm_annotated = annotate_regions(
  regions = dm_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)
# A GRanges object is returned
print(dm_annotated)

# Randomize the input regions
dm_random_regions = randomize_regions(
  regions = dm_regions,
  allow.overlaps = TRUE,
  per.chromosome = TRUE)

# Annotate the random regions using the same annotations as above
# These will be used in later functions
dm_random_annotated = annotate_regions(
  regions = dm_random_regions,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = TRUE)

# Find the number of regions per annotation type
# and the number of random regions per annotation type
dm_annsum_rnd = summarize_annotations(
  annotated_regions = dm_annotated,
  annotated_random = dm_random_annotated,
  quiet = TRUE)
print(dm_annsum_rnd)

# Count the occurrences of classifications in the DM_status
# column across the annotation types.
dm_catsum = summarize_categorical(
  annotated_regions = dm_annotated,
  by = c('annot.type', 'Age'),
  quiet = TRUE) %>%
  pivot_wider(names_from=Age, values_from=n, id_cols=annot.type)
print(dm_catsum)

annots_order = unique(dm_annotated$annot$type)

plot_annotation(annotated_regions = dm_annotated,
  annotated_random = dm_random_annotated,
  annotation_order = annots_order,
  plot_title = 'Dist. of Sites Tested for DM (with rndm.)',
  x_label = 'Annotations',
  y_label = 'Count')

plot_coannotations(
  annotated_regions = dm_annotated,
  annotation_order = annots_order,
  axes_label = 'Annotations',
  plot_title = 'Regions in Pairs of Annotations')

# plot number of enhancers
dm_annotated %>%
  as_tibble() %>%
  filter(annot.type == "mm10_enhancers_fantom") %>%
  group_by(Age) %>%
  summarise(count = n())

# collect genes from promoter overlap
young_promoters_overlap = dm_annotated %>%
  as_tibble() %>%
  filter(annot.type == "mm10_genes_promoters" & Age=="Y") %>%
  pull("annot.symbol") %>%
  na.omit()
aged_promoters_overlap = dm_annotated %>%
  as_tibble() %>%
  filter(annot.type == "mm10_genes_promoters" & Age=="A") %>%
  pull("annot.symbol") %>%
  na.omit()

length(intersect(young_promoters_overlap, aged_promoters_overlap))
length(setdiff(young_promoters_overlap, aged_promoters_overlap))
length(setdiff(aged_promoters_overlap, young_promoters_overlap))

# get promoter and non-promoter boundaries
promoters = dm_annotated %>%
  as_tibble() %>%
  filter(annot.type == "mm10_genes_promoters" & complete.cases(annot.symbol))
list(Y.pro = unique(promoters[promoters$Age=="Y","bound.name"]),
     A.pro = unique(promoters[promoters$Age=="A","bound.name"]))

# Write bedpe files -------------------------------------------------------


# write first 2 lines of contact domain bedpe file
fname = file.path(projdir, sprintf("%s_%s_TAD.40kb.bedpe","tadtool","young"))
writeLines(c(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                     "#chr1","x1","x2","chr2","y1","y2","name","score","strand1","strand2","color","score","uVarScore","lVarScore","upSign","loSign"),
             "# OnTAD"),
           fname, sep = "\n")
# write remaining lines of contact domain bedpe file
bedpe.data = data.frame(cbind(young.tad[,1:3], 
                              young.tad[,1:3],
                              matrix('.',nrow(young.tad),4), 
                              "0,0,1", 
                              matrix(0,nrow(young.tad),5)))
colnames(bedpe.data) = c("chr1","x1","x2","chr2","y1","y2","name","score","strand1","strand2","color","score","uVarScore","lVarScore","upSign","loSign")
write.table(bedpe.data, fname, sep="\t", quote=F, row.names=F, col.names = F, append=T)
# store bedpe data in list to combine all chromosomes into one file


