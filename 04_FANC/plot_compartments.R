library(stringr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap) 
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(ggh4x)
library(dplyr)
library(rstatix)
library(latex2exp)

#projdir = "/nas/homes/benyang/HiC/04_FANC"
projdir = "C:/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/04_FANC"
aged.dir = "aged.merged"
young.dir = "young.merged"

list.files(file.path(projdir, aged.dir))


# Find compartment changes in 250kb data ----------------------------------

young.data = read.table(file.path(projdir, "without KR normalization", "young.merged_100kb_noKR.ev"))
aged.data = read.table(file.path(projdir, "without KR normalization", "aged.merged_100kb_noKR.ev"))

#barplot(height=as.numeric(aged.data$V5), border=rgb(1, 0, 0, .4))
#barplot(height=as.numeric(young.data$V5), border=rgb(0, 1, 0, .4), add=T)

#young.data = read.table(file.path(projdir,"young.merged_250kb.eigenvector"))
young.data$Age = "Y"
#aged.data = read.table(file.path("C:/Users/benjy/Downloads/aged.merged_250kb_noKR.ev"))
aged.data$Age = "A"
data = rbind(young.data, aged.data)
colnames(data) = c("chr","start","end","compartment","eigen","score","Age")
data$chr = factor(data$chr, levels=unique(str_sort(data$chr, numeric=T)))
data$Age = factor(data$Age, levels=c("Y","A"))
data[,c(2,3)] = apply(data[,c(2,3)], 2, as.integer) # convert start and end regions to integers

# find where the first nonzero entry is in eigenvector
firstNonZero = data.frame(matrix(nrow=length(unique(data$chr)), ncol=length(unique(data$Age))))
rownames(firstNonZero) = unique(data$chr)
colnames(firstNonZero) = unique(data$Age)
for (chr in unique(data$chr)) {
  for (a in unique(data$Age)) {
    firstNonZero[chr,a] = which(!data[data$chr==chr & data$Age==a,"eigen"]==0)[1]
  }
}
firstNonZero$both = apply(firstNonZero,1,function(x) max(x['A'],x['Y']))

# trim beginning non-zero entries from eigenvectors
data_trim = list()
for(chrom in rownames(firstNonZero)){
  if(!(chrom %in% c("chrM","chrY"))) {
    y.tmp = data %>%
      dplyr::filter(chr==chrom & Age=="Y") %>%
      slice(as.integer(firstNonZero[chrom,"both"]):n())
    a.tmp = data %>%
      dplyr::filter(chr==chrom & Age=="A") %>%
      slice(as.integer(firstNonZero[chrom,"both"]):n())
    data_trim[[chrom]] = rbind(y.tmp, a.tmp)
  }
}
data_trim = do.call(rbind, data_trim)
data_trim = data_trim %>% 
  mutate(height=1) %>%
  mutate_if(is.character, as.factor)

# data_trim = data %>%
#   dplyr::filter(chr!="chrM" & chr!="chrY") %>%
#   group_by(Age,chr) %>%
#   slice(13:n()) %>% # remove the first 13 rows since they're all zeros 
#   ungroup() %>%
#   mutate(chr = factor(chr, levels=unique(chr)),
#          height = 1) %>%
#   mutate_if(is.character, as.factor) %>%
#   data.frame()

# check that there are the same number of regions per chromosome between aged and young
all(sapply(unique(data_trim$chr), function(x) nrow(subset(data_trim, Age=="Y" & chr==x)) == nrow(subset(data_trim, Age=="A" & chr==x))))


# check average size of compartments --------------------------------------

base.scale = 100
ab.boundaries = data_trim %>%
  group_by(chr, Age) %>%
  summarise(#ab.switch =  list(which(diff(as.numeric(compartment)) != 0)),
    ab.switch.runs = list(rle(as.character(compartment))$lengths), # get lengths of compartment runs
    ab.switch.idx = list(cumsum(rle(as.character(compartment))$lengths)), # get runs of compartments and add to get switch indices
    ab.switch.comp = list(rle(as.character(compartment))$values), # get runs of compartments,
    first.comp = compartment[1],
    chr.start = start[1],
    chr.end = end[length(end)]) %>%
  mutate(num.switch = unlist(lapply(ab.switch.runs, function(x) length(x)-1)))
ab.boundaries$ab.sizes = apply(ab.boundaries, 1, function(x) as.integer(base.scale * 1e3 * x$ab.switch.runs))
# separate ab.sizes by compartment
a.tmp = list()
b.tmp = list()
for (i in 1:nrow(ab.boundaries)){
  tmp.sizes = unlist(ab.boundaries[i,"ab.sizes"], use.names=F)
  if (ab.boundaries[i,"first.comp"]=="A") {
    a.tmp[[i]] = tmp.sizes[seq(1,length(tmp.sizes),2)]
    b.tmp[[i]] = tmp.sizes[seq(2,length(tmp.sizes),2)]
  } else {
    a.tmp[[i]] = tmp.sizes[seq(2,length(tmp.sizes),2)]
    b.tmp[[i]] = tmp.sizes[seq(1,length(tmp.sizes),2)]
  }
}
ab.boundaries$b.sizes = b.tmp
ab.boundaries$a.sizes = a.tmp

# plot number of switches per chr per age
ggplot(ab.boundaries, aes(x=chr, y=num.switch)) +
  geom_col(aes(fill=Age), color="black", position=position_dodge(0.9)) +
  #scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_fill_grey() +
  scale_x_discrete(labels=c(1:19, "X")) +
  scale_y_continuous(breaks = seq(0,200,25), expand = expansion(mult=c(0,0.1))) +
  xlab("Chromosome") + ylab("Number of A/B compartment switches") +
  theme_bw() +
  theme(axis.text = element_text(size=12, color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "top")
ggsave(file.path(projdir,  "compartment size", sprintf("aged_young_%dkb_compartment_num_switches.png",base.scale)), dpi=300, width=7, height=5)

ab.boundaries.sizes = stack(list(Aged_A = unlist(subset(ab.boundaries, Age=="A", select=a.sizes), use.names=F),
                                 Aged_B = unlist(subset(ab.boundaries, Age=="A", select=b.sizes), use.names=F),
                                 Young_A = unlist(subset(ab.boundaries, Age=="Y", select=a.sizes), use.names=F),
                                 Young_B = unlist(subset(ab.boundaries, Age=="Y", select=b.sizes), use.names=F)))
ab.boundaries.sizes = ab.boundaries.sizes %>%
  separate(ind, into=c("Age","compartment"), sep="_", remove=F) %>%
  mutate_if(is.character, as.factor)
levels(ab.boundaries.sizes$Age) = list(A="Aged", Y="Young")
ab.boundaries.sizes$Age = factor(ab.boundaries.sizes$Age, levels=c("Y","A"))

num.bins.test = ab.boundaries.sizes %>%
  mutate(log10size = log10(values)) %>%
  group_by(compartment) %>%
  wilcox_test(log10size~Age) %>%
  add_xy_position(x="compartment", dodge=0.9)

num.bins.test2 = ab.boundaries.sizes %>%
  mutate(log10size = log10(values)) %>%
  group_by(Age) %>%
  wilcox_test(log10size~compartment) %>%
  add_xy_position(x="compartment", group="Age", dodge=0.9)
num.bins.test2 = num.bins.test2 %>% mutate(y.position=c(8,7.5))

ggplot(ab.boundaries.sizes, aes(x=compartment, y=log10(values))) +
  geom_violin(aes(fill=Age), color="black") +
  geom_boxplot(aes(group=interaction(compartment,Age)), position=position_dodge(0.9), width=0.2, color="black", outlier.shape=NA) +
  scale_fill_manual(values=rev(pal_nejm()(2)), labels=c("Y"="Young", "A"="Aged")) +
  stat_pvalue_manual(num.bins.test, tip.length = 0.01, label.size=4, bracket.nudge.y = 0.05) +
  stat_pvalue_manual(num.bins.test2, tip.length = 0.01, label.size=4) +
  scale_y_continuous(breaks=seq(4,8,0.5), expand=expansion(mult=c(0.1,0.1))) +
  theme_bw() + xlab(NULL) + ylab(expression(bold(paste("log"["10"],"(compartment size in bases)")))) +
  theme(legend.position = "top",
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir,  "compartment size", sprintf("aged_young_%dkb_compartment_size.png",base.scale)), dpi=300, width=3, height=4)

# write BED files

young.ab.tmp = list()
aged.ab.tmp = list()
for (i in unique(ab.boundaries$chr)){
  # Young
  tmp.switch.idx = unlist(subset(ab.boundaries, chr==i & Age=="Y", "ab.switch.idx"), use.names=F)
  tmp.ends = c(unlist(subset(ab.boundaries, chr==i & Age=="Y" ,"chr.start")) + 100e3 * tmp.switch.idx[1:length(tmp.switch.idx)-1] - 1,
               unlist(subset(ab.boundaries, chr==i & Age=="Y" ,"chr.end"), use.names=F))
  tmp.starts = c(unlist(subset(ab.boundaries, chr==i & Age=="Y" ,"chr.start"), use.names=F),
                 unlist(subset(ab.boundaries, chr==i & Age=="Y" ,"chr.start"), use.names=F) + 100e3 * tmp.switch.idx[1:length(tmp.switch.idx)-1])
  tmp.comp = unlist(ifelse(subset(ab.boundaries, chr==i & Age=="Y", "first.comp")=="A",
                           list(rep(c("A","B"),times=length(tmp.starts))),
                           list(rep(c("B","A"),times=length(tmp.starts)))))
  young.ab.tmp[[i]] = data.frame(chr=rep(i, length(tmp.starts)),
                                 start=tmp.starts,
                                 end=tmp.ends,
                                 compartment=tmp.comp)
  # Aged
  tmp.switch.idx = unlist(subset(ab.boundaries, chr==i & Age=="A", "ab.switch.idx"), use.names=F)
  tmp.ends = c(unlist(subset(ab.boundaries, chr==i & Age=="A" ,"chr.start")) + 100e3 * tmp.switch.idx[1:length(tmp.switch.idx)-1] - 1,
               unlist(subset(ab.boundaries, chr==i & Age=="A" ,"chr.end"), use.names=F))
  tmp.starts = c(unlist(subset(ab.boundaries, chr==i & Age=="A" ,"chr.start"), use.names=F),
                 unlist(subset(ab.boundaries, chr==i & Age=="A" ,"chr.start"), use.names=F) + 100e3 * tmp.switch.idx[1:length(tmp.switch.idx)-1])
  tmp.comp = unlist(ifelse(subset(ab.boundaries, chr==i & Age=="A", "first.comp")=="A",
                           list(rep(c("A","B"),times=length(tmp.starts))),
                           list(rep(c("B","A"),times=length(tmp.starts)))))
  aged.ab.tmp[[i]] = data.frame(chr=rep(i, length(tmp.starts)),
                                start=tmp.starts,
                                end=tmp.ends,
                                compartment=tmp.comp)
}
young.ab.bed = do.call(rbind, young.ab.tmp)
aged.ab.bed = do.call(rbind, aged.ab.tmp)

young.ab.bed[,c("start","end")] = apply(young.ab.bed[,c("start","end")], 2, as.integer) 
aged.ab.bed[,c("start","end")] = apply(aged.ab.bed[,c("start","end")], 2, as.integer) 

write.table(tibble(young.ab.bed),
            file = file.path(projdir, "without KR normalization", "young.ab_100kb.bed"),
            sep="\t", quote = F, col.names = F, row.names = F)
write.table(tibble(aged.ab.bed),
            file = file.path(projdir, "without KR normalization", "aged.ab_100kb.bed"),
            sep="\t", quote = F, col.names = F, row.names = F)

young.num.bins %>%
  group_by(compartment) %>%
  summarise(summary = quantile(size),
            mean = mean(size))
aged.num.bins %>%
  group_by(compartment) %>%
  summarise(summary = quantile(size),
            mean = mean(size))


# Export compartment bedGraph for each sample -----------------------------------

write.table(subset(data_trim, Age=="Y", c("chr","start","end","eigen")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "young.ab.bedGraph"),
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(subset(data_trim, Age=="A", c("chr","start","end","eigen")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "aged.ab.bedGraph"),
            sep="\t", quote=F, row.names=F, col.names=F)

# Export compartment BED for each group -----------------------------------

write.table(subset(young.ab.bed,compartment=="A", c("chr","start","end")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "young.A.bed"),
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(subset(young.ab.bed,compartment=="B", c("chr","start","end")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "young.B.bed"),
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(subset(aged.ab.bed,compartment=="A", c("chr","start","end")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "aged.A.bed"),
            sep="\t", quote=F, row.names=F, col.names=F)
write.table(subset(aged.ab.bed,compartment=="B", c("chr","start","end")),
            file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "aged.B.bed"),
            sep="\t", quote=F, row.names=F, col.names=F)


# create list of indices for A and B compartments in young sample  --------

young.A.list = tapply(X=with(data_trim, compartment[Age=="Y"]), 
                      INDEX=with(data_trim, chr[Age=="Y"]), 
                      FUN=function(x) which(x=="A"))
young.B.list = tapply(X=with(data_trim, compartment[Age=="Y"]), 
                      INDEX=with(data_trim, chr[Age=="Y"]), 
                      FUN=function(x) which(x=="B"))
young.compartment.idx = data.frame(A = young.A.list, B = young.B.list)

length(unlist(young.compartment.idx["chr2",]))
nrow(data_trim[data_trim$Age=="Y" & data_trim$chr=="chr2",])


# Get compartment switches using merge data frames ------------------------

chrom.sizes = read.table("C:/Users/benjy/Dropbox (University of Michigan)/ENGIN-Lab Notes/Lab Notes/Lab Notes Benjamin/Hi-C/get_tss/sizes.mm10",sep="\t")

data_trim_joined = full_join(data_trim[data_trim$Age=="Y",],data_trim[data_trim$Age=="A",], by=c("chr","start","end"), suffix=c(".Young",".Aged")) %>% drop_na() %>% distinct()
#data_trim_joined = merge(data_trim[data_trim$Age=="Y",],data_trim[data_trim$Age=="A",], by=c("chr","start","end"), all=T, suffixes=c(".Young",".Aged"))
data_trim_joined = data_trim_joined %>% 
  mutate(Category = apply(data_trim_joined,1,function(x) ifelse(is.na(x["compartment.Aged"]) | is.na(x["compartment.Young"]), "NA",
                                                                ifelse(x["compartment.Aged"]=="B" & x["compartment.Young"]=="A", "AtoB",
                                                                       ifelse(x["compartment.Aged"]=="A" & x["compartment.Young"]=="B", "BtoA",
                                                                              ifelse(x["compartment.Aged"]=="A" & x["compartment.Young"]=="A", "staticA", "staticB"))))))

cat.tbl = table(data_trim_joined$Category)

data_trim_joined %>% group_by(chr) %>% summarise(min_start=min(start), max_start=max(start), min_end=min(end), max_end=max(end)) %>%
  mutate(size=sapply(chr,function(x) chrom.sizes[chrom.sizes$V1==x,]$V2),
         over=(max_end>size) | (max_start>size),
         flip=(max_start>max_end))

data_trim_joined %>%
  ggplot(aes(x=eigen.Young,y=eigen.Aged)) +
  geom_point(aes(color=Category), size=0.1) +
  scale_color_uchicago(name=NULL, labels=c(TeX("$A \\rightarrow B$"),TeX("$B \\rightarrow A$"),TeX("Static A"),TeX("Static B"))) +
  scale_x_continuous(breaks=seq(-0.1,0.1,0.02)) +
  scale_y_continuous(breaks=seq(-0.1,0.1,0.02)) +
  theme_bw() +
  coord_fixed() +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  labs(x="Young Eigenvector", y="Aged Eigenvector") +
  theme(axis.text=element_text(size=12,color="black",face="bold"),
        axis.title=element_text(size=12,color="black",face="bold"),
        legend.position="top",
        legend.text=element_text(size=12,color="black"))
ggsave(file.path(projdir,"compartment_score_scatterplot.png"),dpi=300,width=4,height=4)

# data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="AtoB") %>%
#   write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "A_to_B.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
# data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="BtoA") %>%
#   write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "B_to_A.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
# data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticA") %>%
#   write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "StaticA.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
# data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticB") %>%
#   write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "StaticB.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
# data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticB" | values=="staticA") %>%
#   write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "static.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

head(data_trim_joined)
data_trim_joined_rle = data_trim_joined %>% drop_na() %>% group_by(chr) %>% 
  summarise(lengths=rle(Category)$lengths, values=rle(Category)$values) %>% mutate(sumLengths = cumsum(lengths))
data_trim_starts = data_trim %>% group_by(Age,chr) %>% summarise(start=start[1]) %>% pivot_wider(names_from=Age,values_from=start)
apply(data_trim_starts,1,function(x) x['Y']==x['A']) # check which bins have different starts 
data_trim_joined_rle$end = as.integer(apply(data_trim_joined_rle,1,function(x) data_trim_starts[data_trim_starts$chr==x["chr"],]$A + as.numeric(x["sumLengths"]) * 100e3 - 1))
tmp_start_list=list()
for(chr in data_trim_joined_rle$chr) {
  tmp_start_list[[chr]] = c(data_trim_starts[data_trim_starts$chr==chr,]$A, 
                            data_trim_starts[data_trim_starts$chr==chr,]$A + 
                              100e3*pull(data_trim_joined_rle[data_trim_joined_rle$chr==chr,"sumLengths"])[1:sum(data_trim_joined_rle$chr==chr)-1])
}
data_trim_joined_rle$start = as.integer(unlist(tmp_start_list,use.names = F))
# need to replace max end coordinate per chromosome with chromosome size
max_idx = data_trim_joined_rle %>% group_by(chr) %>% summarise(max_end_idx = which(end==max(end)))
for(chr in max_idx$chr) {
  data_trim_joined_rle[data_trim_joined_rle$chr==chr,]$end[max_idx[max_idx$chr==chr,]$max_end_idx] = chrom.sizes[chrom.sizes$V1==chr,]$V2
} 

data_trim_joined_rle %>% group_by(chr) %>% summarise(min_start=min(start), max_start=max(start), min_end=min(end), max_end=max(end)) %>%
  mutate(size=sapply(chr,function(x) chrom.sizes[chrom.sizes$V1==x,]$V2),
         over=(max_end>size) | (max_start>size),
         flip=(max_start>max_end))

data_trim_joined_rle %>% select(chr,start,end,values) %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "ab.switch.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="AtoB") %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "A_to_B.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="BtoA") %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "B_to_A.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticA") %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "StaticA.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticB") %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "StaticB.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
data_trim_joined_rle %>% select(chr,start,end,values) %>% filter(values=="staticB" | values=="staticA") %>%
  write.table(file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "static.bed"), quote = F, sep = "\t", row.names = F, col.names = F)


# create lists of changes in aged compared to young -----------------------


A.to.B = lapply(unique(data_trim$chr), function(x) which(with(data_trim, compartment[Age=="A" & chr==x])[unlist(young.compartment.idx[x,"A"])] == "B"))
B.to.A = lapply(unique(data_trim$chr), function(x) which(with(data_trim, compartment[Age=="A" & chr==x])[unlist(young.compartment.idx[x,"B"])] == "A"))
static.A = lapply(unique(data_trim$chr), function(x) which(with(data_trim, compartment[Age=="A" & chr==x])[unlist(young.compartment.idx[x,"A"])] == "A"))
static.B = lapply(unique(data_trim$chr), function(x) which(with(data_trim, compartment[Age=="A" & chr==x])[unlist(young.compartment.idx[x,"B"])] == "B"))
# combine lists of static regions
static = list()
for(i in seq_along(static.A)) {
  static[[i]] = c(unlist(static.A[i]), unlist(static.B[i]))
}



# create BED files for each group of compartment change -------------------

A.to.B.bed = list()
B.to.A.bed = list()
static.bed = list()
for(i in seq_along(unique(data_trim$chr))) {
  chr = unique(data_trim$chr)[i]
  if(length(A.to.B[[i]])>1) {
    A.to.B.bed[[chr]] = cbind(chr, data_trim[data_trim$chr==chr, c("start","end")][A.to.B[[i]],])
  }
  if(length(B.to.A[[i]])>1) {
    B.to.A.bed[[chr]] = cbind(chr, data_trim[data_trim$chr==chr, c("start","end")][B.to.A[[i]],])
  }
  if(length(static[[i]])>1) {
    static.bed[[chr]] = cbind(chr, data_trim[data_trim$chr==chr, c("start","end")][static[[i]],]) 
  }
}
write.table(do.call(rbind, A.to.B.bed), file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "A_to_B.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(do.call(rbind, B.to.A.bed), file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "B_to_A.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(do.call(rbind, static.bed), file.path(projdir, "compartmentExpression", "compartmentBed", "100kb", "static.bed"), quote = F, sep = "\t", row.names = F, col.names = F)


# plot fraction of compartment switches -----------------------------------

# switch.data = data.frame(Frac = c(length(unlist(A.to.B)) / (0.5*nrow(data_trim)), #0.5 because young and aged merged together in data_trim
#                                   length(unlist(B.to.A)) / (0.5*nrow(data_trim)),
#                                   length(unlist(static)) / (0.5*nrow(data_trim))),
#                          col = rev(pal_jama()(3)),
#                          row.names = c("A.to.B","B.to.A","Static"))
switch.data = data_trim_joined_rle %>% group_by(values) %>% summarise(num=sum(lengths)) %>% filter(values!="NA") %>% 
  mutate(Frac=num/sum(num),
         col=rev(pal_jama()(6)[c(1,2,5,6)])) %>%
  as.data.frame()
switch.data = switch.data[c(1,2,4,3),]

png(file.path(projdir, "compartmentExpression", "compartment_fraction_100kb.png"), res=300, units="in", width=5, height=3)
par(font.lab=2, font.axis=2, mar = c(5, 1, 0, 1))
barplot(switch.data$Frac, xlim=c(0,1), 
        col=switch.data$col,
        xlab="Compartment Fraction", horiz = T, beside = T,axes = T,
        cex.axis = 1.5, cex.lab = 1.5)
legend("bottomright",
       legend=c(TeX("Static A "),TeX("Static B "),TeX("$B \\rightarrow A$ "),TeX("$A \\rightarrow B$ ")),
       col=rev(switch.data$col),
       pch=15, cex=1.5)
#text(0.7,1,"hello")
dev.off()


# Get fraction of genome --------------------------------------------------

#eff_mm10_size = 2494787188
eff_mm10_size=2652783500
sum_list = list()
sum_list[1] = ab.boundaries %>% filter(Age=="A") %>% select(a.sizes) %>% pull() %>% unlist() %>% sum()
sum_list[2] = ab.boundaries %>% filter(Age=="A") %>% select(b.sizes) %>% pull() %>% unlist() %>% sum()
sum_list[3] = ab.boundaries %>% filter(Age=="Y") %>% select(a.sizes) %>% pull() %>% unlist() %>% sum()
sum_list[4] = ab.boundaries %>% filter(Age=="Y") %>% select(b.sizes) %>% pull() %>% unlist() %>% sum()
names(sum_list) = c("Aged_A","Aged_B","Young_A","Young_B")

frag_genome_df = data.frame(comp_size = c(unlist(sum_list), unlist(switch.size)))
frag_genome_df$genome_size = c(sum(frag_genome_df$comp_size[1:2]), sum(frag_genome_df$comp_size[1:2]),
                               sum(frag_genome_df$comp_size[3:4]), sum(frag_genome_df$comp_size[3:4]),
                               sum(frag_genome_df$comp_size[5:7]), sum(frag_genome_df$comp_size[5:7]), sum(frag_genome_df$comp_size[5:7]))
frag_genome_df$frac_genome = frag_genome_df$comp_size/frag_genome_df$genome_size

# plot compartment genome fraction
frag_genome_df %>% slice(1:4) %>% mutate(Age=factor(c("Aged","Aged","Young","Young"),levels=c("Young","Aged")),
                                         Compartment=c("A","B","A","B")) %>%
  ggplot(aes(x=Age,y=frac_genome,fill=Age)) +
  geom_col_pattern(aes(pattern=Compartment),position="stack",
                   colour="black",
                   pattern_density = 0.1,
                   pattern_fill    = 'white',
                   pattern_colour  = 'black') +
  scale_fill_manual(values=rev(pal_nejm()(2))) +
  scale_pattern_manual(values=c("none","stripe")) +
  guides(fill="none") +
  scale_y_continuous(expand = expansion(mult=c(0,0.1)), breaks=seq(0,1,0.25)) +
  theme_bw() + 
  labs(x=NULL,y="Genome Fraction") +
  theme(legend.position="top",
        legend.text=element_text(size=14),
        legend.title=element_text(size=14,face="bold"),
        axis.text=element_text(size=14,face="bold",color="black"),
        axis.title=element_text(size=14,face="bold"))
ggsave(file.path(projdir,"AB_compartment_genome_fraction.png"),dpi=300,width=3.25,height=5)

# plot compartments per chromosome ----------------------------------------

#png(file.path(projdir, "compartmentExpression", "compartment_per_chromosome_100kb.png"), res=200, units="in", width=4, height=10)
ggplot(data_trim, aes(x=start, y=height)) +
  facet_nested(chr+Age~., switch="y") +
  #facet_wrap(facets = vars(chr), nrow=length(unique(data_trim$chr)), ncol=1) +
  geom_bar(aes(fill=compartment), stat="identity") +
  scale_x_continuous(breaks=pretty(seq(min(data_trim[data_trim$Age=="A","start"]),
                                       max(data_trim[data_trim$Age=="A","start"])), 
                                   n=20)) +
  scale_y_continuous(limits=c(-0.5,1.5)) +
  #scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  theme_void() + 
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text.y.left = element_text(size=10),
        legend.position = "top",
        legend.key.height = unit(1, 'mm'),
        legend.key.width = unit(1, 'mm'))
#dev.off()
ggsave(file.path(projdir, "compartment_per_chromosome_100kb.png"), dpi=300, width=4, height=10)

# Plot chromosome level compartment changes -------------------------------

head(data)

#png(file.path(projdir, "young_aged_compartments_100kb.png"), res=200, units="in", width=12, height=4)
data %>%
  #filter(Age=="Y") %>%
  filter(!(chr %in% c("chrM","chrY"))) %>%
  ggplot(aes(x=start, y=eigen)) +
  facet_wrap(vars(chr), nrow = 5, ncol = 5, scales="free_x") +
  geom_bar(stat="identity", aes(fill=Age), alpha=0.6, position="identity") +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)], labels=c("Y"="Young","A"="Aged")) +
  #scale_fill_npg() +
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "top",
        legend.title = element_text(face='bold', size=12),
        legend.text = element_text(size=12),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#dev.off()
ggsave(file.path(projdir, "young_aged_compartments_100kb.png"), dpi=300, width=9, height=4)

library(plotgardener)

ab.switch.bed = read.table(file.path(projdir,"compartmentExpression","compartmentBed","100kb","ab.switch.bed"))
ab.switch.bed.subset = ab.switch.bed[ab.switch.bed$V1=="chr10" & ab.switch.bed$V2>=25e6 & ab.switch.bed$V3<=40e6, ]
ab.switch.bed.subset = ab.switch.bed.subset[!ab.switch.bed.subset$V4 %in% c("staticA","staticB"), ]
params = pgParams(assembly='mm10', chrom="chr10", chromstart=25000000, chromend=40000000, x=1)

png(file.path(projdir, 'young_aged_compartments_chr1_100kb_plotgardener.png'), res=300, width=5, height=3, units='in')
pageCreate(width=5, height=3)
young_plt = plotSignal(setNames(data[data$Age=="Y",c("chr","start","end","eigen")], c("chr","start","end","score")), binCap = F,
                       baseline.color = pal_nejm()(2)[2], fill = pal_nejm()(2)[2], linecolor = pal_nejm()(2)[2], 
                       y=0.5, negData = T, range = c(-0.05,0.05), height=1, width=4, params=params)
aged_plt = plotSignal(setNames(data[data$Age=="A",c("chr","start","end","eigen")], c("chr","start","end","score")), binCap = F,
                       baseline.color = pal_nejm()(2)[1], fill = pal_nejm()(2)[1], linecolor = pal_nejm()(2)[1], 
                       y="b0", negData = T, range = c(-0.05,0.05), height=1, width=4, params=params)
annoYaxis(
  plot = young_plt, at = c(-0.025,0,0.025),
  axisLine = T, fontsize = 10
)
annoYaxis(
  plot = aged_plt, at = c(-0.025,0,0.025),
  axisLine = T, fontsize = 10
)
for(i in 1:nrow(ab.switch.bed.subset)) {
  annoHighlight(
    plot = young_plt, 
    chrom = ab.switch.bed.subset[i,'V1'], chromstart = ab.switch.bed.subset[i,'V2'], chromend = ab.switch.bed.subset[i,'V3'],
    y = young_plt$y, height = 2, fill = "#7ecdbb", alpha=0.4
  ) 
}
annoGenomeLabel(aged_plt, y="b0", fontsize=14, fontface='bold', scale = "Mb", params=params)
legendPlot <- plotLegend(
  legend = c("Young","Aged"),
  fill = pal_nejm()(2)[2:1],
  orientation = 'h', fontsize = 14,
  border = FALSE, params = params, 
  x = 2.25, y = 0, width = 2, height = 0.7
)
pageGuideHide()
dev.off()

#png(file.path(projdir, "young_aged_compartments_100kb.png"), res=200, units="in", width=12, height=4)
data %>%
  filter(chr=="chr10") %>%
  ggplot(aes(x=start, y=eigen)) +
  facet_grid(rows="Age") +
  geom_bar(stat="identity", aes(fill=Age), position="identity") +
  geom_segment(aes(x=data[data$chr=="chr10" & data$eigen!=0,"start"][1], 
                   xend=tail(data[data$chr=="chr10","end"],1), 
                   y=0, 
                   yend=0), 
               color="black") +
  scale_fill_manual(values=pal_nejm()(2)[c(2,1)]) +
  scale_x_continuous(limits=c(2e7,4e7),
                     breaks=seq(0,tail(data[data$chr=="chr10","end"],1),1e7)) +
  theme_void() + 
  xlab("") + ylab("") 
theme(legend.position = "top",
      legend.title = element_text(size=12),
      legend.text = element_text(size=12),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.text = element_blank())
#dev.off()
ggsave(file.path(projdir, "young_aged_compartments_chr1_100kb.png"), dpi=300, width=5, height=3)

data_sum = data_trim %>%
  group_by(chr) %>%
  summarise(eq_frac = sum(compartment[Age=="Y"] == compartment[Age=="A"])/length(compartment[Age=="Y"]),
            diff_frac = sum(compartment[Age=="Y"] != compartment[Age=="A"])/length(compartment[Age=="Y"]),
            equal = sum(compartment[Age=="Y"] == compartment[Age=="A"]),
            diff = sum(compartment[Age=="Y"] != compartment[Age=="A"]),
            num.Y = length(compartment[Age=="Y"]),
            num.A = length(compartment[Age=="A"])) 

data %>%
  summarise(A.young.frac = sum(compartment[Age=="Y"]=="A")/length(compartment[Age=="Y"]),
            B.young.frac = sum(compartment[Age=="Y"]=="B")/length(compartment[Age=="Y"]),
            A.aged.frac = sum(compartment[Age=="A"]=="A")/length(compartment[Age=="A"]),
            B.aged.frac = sum(compartment[Age=="A"]=="B")/length(compartment[Age=="A"]))

ggplot(data_sum, aes(x=chr, y=diff_frac)) +
  geom_bar(stat="identity") +
  theme_bw() + xlab("") + ylab("Fraction of A/B compartment change") +
  scale_y_continuous(expand=expansion(mult=c(0,0.1))) +
  ggtitle("250kb") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

png(file.path(projdir, "fanc_250kb_young_aged_eigenvector.png"), res=200, units="in", width=7, height=5)
data %>%
  dplyr::filter(chr != "chrY" & chr != "chrM") %>%
  ggplot(aes(x=V2, y=V5)) +
  facet_wrap(facets="chr", nrow=5, ncol=5, scales = "free_x") +
  geom_bar(stat="identity", position="identity", alpha=0.5, aes(fill=Age)) +
  #geom_path(alpha=0.4, aes(color=Age)) +
  #scale_y_continuous(limits = c(-0.14, 0.14), breaks=seq(-0.12,0.12,0.08)) +
  scale_fill_nejm() +
  theme_bw() + ylab("") + xlab("") +
  ggtitle("250kb") +
  theme(axis.text.x = element_blank(),
        legend.position = "top")
dev.off()

eigen_files = grep(list.files(projdir), pattern="eigen", fixed=T, value = T)
pearsons_files = grep(list.files(projdir), pattern="pearsons", fixed=T, value = T)

pearsons_list = list()
for (f in pearsons_files) {
  pearsons_list[[f]] = read.table(file.path(projdir, f), sep=" ")
}

eigen_list = list()
for (f in eigen_files) {
  eigen_list[[f]] = read.table(file.path(projdir, f), sep=" ")
}

for (i in 1:length(pearsons_files)) {
  tmp = HeatmapAnnotation(`A/B` = anno_barplot(c(NA,unlist(eigen_list[[i]])),
                                               ylim = c(-0.125,0.125)))
  tmp1 = Heatmap(as.matrix(pearsons_list[[i]]), name="Pearsons", 
                 cluster_rows = F, cluster_columns = F,
                 show_row_names = F, show_column_names = F,
                 bottom_annotation = tmp,
                 col = colorRampPalette(rev(brewer.pal(n=7, "RdBu")))(100))
  png(file.path(projdir,"Figures",paste0(pearsons_files[[i]],".png")),width=5,height=5,units="in",res=200)
  draw(tmp1)
  dev.off()
}

unlist(lapply(eigen_list, function(x) {abs(max(as.numeric(unlist(x)), na.rm=T))}))
