library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(ggpattern)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\05_Arrowhead_HICCUPS\\HOMER"

files = list.files(projdir, pattern = "CTCF.bed")
random_files = list.files(projdir, pattern = "random.bed")

data=list()
for(f in files) {
  fname = gsub(".bed","",f,fixed=T)
  tmp = read.table(file.path(projdir,f))
  colnames(tmp) = c("chr","start","end","ID","chr2","start2","end2","CTCF","score","strand","overlap")
  tmp$fname = fname
  tmp$class = "Sample"
  tmp = tmp %>% separate(fname, into=c("Age","Direction"), sep="_", remove=F)
  data[[fname]] = tmp
}
random_data=list()
for(f in random_files) {
  fname = gsub(".bed","",f,fixed=T)
  tmp = read.table(file.path(projdir,f))
  colnames(tmp) = c("chr","start","end","ID","width","strand1","chr2","start2","end2","CTCF","score","strand","overlap")
  tmp$fname = fname
  tmp$class = "Random"
  tmp = tmp %>% separate(fname, into=c("Age","Direction",NA,NA), sep="_", remove=F)
  random_data[[fname]] = tmp %>% dplyr::select(chr,start,end,ID,chr2,start2,end2,CTCF,score,strand,overlap,fname,Age,Direction,class)
}

plt.df = rbind(do.call(rbind,data), do.call(rbind,random_data))
plt.df = plt.df %>% mutate(Group = paste(Age,ID,Direction,class,sep="_"),
                           size = end-start)

plt.df.collapsed = list()
for(g in unique(plt.df$Group)){
    tmp = plt.df %>% filter(Group==g)
    max_score = max(tmp$score, na.rm=T)
    idx = which(tmp$score==max_score)
    tmp = tmp %>% dplyr::slice(idx)
    plt.df.collapsed[[g]] = tmp
}
plt.df.collapsed = bind_rows(plt.df.collapsed)
head(plt.df.collapsed)

table(plt.df$strand)

type_list=list()
for(g in plt.df.collapsed$Group) {
  tbl = table(plt.df.collapsed[plt.df.collapsed$Group==g,"strand"])
  if(length(tbl)>1) {
    type = ifelse(unname(tbl[1]==tbl[2]), "+/-", names(sort(tbl,decreasing=TRUE))[1])
  } else {
    type = names(tbl)
  }
  type_list[[g]] = rep(unname(type),sum(tbl))
}
plt.df.collapsed$type = unlist(type_list, use.names = F)

Age.labs <- c("Young", "Aged")
names(Age.labs) <- c("young", "aged")

plt.df.collapsed %>%
  dplyr::select(chr,start,end,Age,Direction,type,class) %>%
  distinct() %>%
  mutate(Direction = factor(Direction, levels=c("start","end")),
         Age = factor(Age, levels=c("young","aged"))) %>%
  group_by(Age,Direction,type,class) %>%
  summarise(count=n()) %>%
  ggplot(aes(x=Direction,y=count)) +
  geom_col_pattern(aes(fill=Age,pattern=type), 
                   position="stack",
                   colour="black",
                   pattern_density = 0.05,
                   pattern_fill    = 'white',
                   pattern_colour  = 'white') +
  facet_grid(cols=vars(Age), rows=vars(class), labeller = labeller(Age=Age.labs)) +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_pattern_manual(values=c("circle","stripe","none")) +
  guides(fill="none",pattern=guide_legend(title="CTCF Direction")) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1)),
                     breaks = seq(0,5000,500)) +
  scale_x_discrete(labels=c("Left","Right")) +
  theme_bw() + 
  labs(x="Loop Anchor Locus",y="# Loop Anchors with CTCF Motif") +
  theme(legend.position="top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=14,face="bold"),
        strip.text=element_text(size=14,face="bold",color="black"))
ggsave(file.path(projdir,"CTCF_direction_barplot.png"),dpi=300,width=4,height=7)

