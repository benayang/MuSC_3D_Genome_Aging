library(dplyr)
library(tidyr)
library(ggplot2)

files = grep(pattern="CTCF.bed", x=list.files(), fixed=T, value=T)

data = list()
for (f in files) {
    fname=paste0(unlist(strsplit(f,split="_",fixed=T))[1:2],collapse="_")
    tmp = read.table(f)
    tmp = tmp[,c(1:3,ncol(tmp))]
    colnames(tmp) = c("chr","start","end","count")
    tmp$fname=fname
    data[[fname]] = tmp
}
plt.df = do.call(rbind,data)
plt.df = plt.df %>% separate(fname, into=c("Age","Group"), sep="_", remove=F)

ggplot(plt.df, aes(x=Age,y=count)) +
geom_boxplot(aes(fill=Age), outlier.size=0.1) +
facet_grid(rows=vars(Group)) +
xlab("") +
theme(legend.position="top") +
theme_bw()
ggsave("CTCF_overlap_boundary_vs_domain.png",dpi=300,width=3,height=5)