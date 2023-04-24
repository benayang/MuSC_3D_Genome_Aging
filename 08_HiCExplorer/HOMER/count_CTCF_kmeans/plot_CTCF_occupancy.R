library(dplyr)
library(tidyr)
library(ggplot2)

files = grep(pattern="CTCF.bed", x=list.files(), fixed=T, value=T)

data = list()
for (f in files) {
    fname=paste0(unlist(strsplit(f,split="_",fixed=T))[1:2],collapse="_")
    tmp = read.table(f)
    tmp$fname=fname   
    data[[fname]] = tmp
}
plt.df = do.call(rbind,data)
colnames(plt.df) = c("chr","start","end","score","count","fname")
plt.df = plt.df %>% separate(fname, into=c("Age","Cluster"), sep="_", remove=F)

ggplot(plt.df, aes(x=Age,y=count)) +
geom_boxplot(aes(fill=Age), outlier.size=0.1) +
facet_grid(rows=vars(Cluster)) +
xlab("") +
theme(legend.position="top") +
theme_bw()
ggsave("CTCF_overlap_count.png",dpi=300,width=3,height=7)

plt.df = plt.df %>% mutate(Group = ifelse(count>3,">3",ifelse(count>=1 && count <=3, "1-3", "0")))

ggplot(plt.df, aes(x=Age,y=score)) +
geom_boxplot(aes(fill=Group), outlier.size=0.1) +
facet_grid(rows=vars(Cluster)) +
xlab("") +
theme(legend.position="top") +
theme_bw()
ggsave("CTCF_overlap_count_score.png",dpi=300,width=3,height=7)
