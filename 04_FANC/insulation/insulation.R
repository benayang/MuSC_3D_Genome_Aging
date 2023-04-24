library(dplyr)
library(ggplot2)
library(ggpubr)

projdir = "C:\\Users\\benjy\\Dropbox (University of Michigan)\\ENGIN-Lab Notes\\Lab Notes\\Lab Notes Benjamin\\Hi-C\\04_FANC\\insulation"

young.files = list.files(projdir, pattern = "young.merged.boundaries.")
aged.files = list.files(projdir, pattern = "aged.merged.boundaries.")

young.data = list()
aged.data = list()
for(i in c("120kb","200kb","280kb","400kb","600kb")) {
  young.data[[i]] = read.table(file.path(projdir, sprintf("young.merged.boundaries_%s.bed", i)), 
                               sep="\t", as.is=T, col.names = c("chr","start","end","something","score","strand"))
  aged.data[[i]] = read.table(file.path(projdir, sprintf("aged.merged.boundaries_%s.bed", i)), 
                              sep="\t", as.is=T, col.names = c("chr","start","end","something","score","strand"))
  young.data[[i]]$res = i
  young.data[[i]]$Age = "Y"
  aged.data[[i]]$res = i
  aged.data[[i]]$Age = "A"
}

young.data.df = do.call(rbind, young.data)
aged.data.df = do.call(rbind, aged.data)
all.df = rbind(young.data.df, aged.data.df)

filt.all.df = all.df %>%
  filter(score>0.6 & res=="400kb")

filt.all.df.ranges = filt.all.df %>%
  group_by(chr) %>%
  summarise(chr.ranges = list(seq(min(start),max(start),40000)))

tmp = list()
for(chrom in unique(filt.all.df.ranges$chr)){
  tmp.data = data.frame(chr=chrom,
                        start = unlist(subset(filt.all.df.ranges, chr==chrom, "chr.ranges"), use.names = F),
                        end = unlist(subset(filt.all.df.ranges, chr==chrom, "chr.ranges"), use.names = F) + 39999,
                        young = 0,
                        aged = 0)
  
  y.idx = tmp.data$start %in% with(filt.all.df, start[chr==chrom & Age=="Y"])
  a.idx = tmp.data$start %in% with(filt.all.df, start[chr==chrom & Age=="A"])

  tmp.data$young[y.idx] = with(filt.all.df, score[chr==chrom & Age=="Y"])
  tmp.data$aged[a.idx] = with(filt.all.df, score[chr==chrom & Age=="A"])
  tmp[[chrom]] = tmp.data
}
tmp.all = do.call(rbind, tmp)
tmp.all[,2:3] = apply(tmp.all[,2:3], 2, as.integer)

# Write BED ---------------------------------------------------------------

head(tmp.all)
bed.to.write = cbind(tmp.all[tmp.all$young!=0, 1:3],
                     ".",
                     tmp.all[tmp.all$young!=0, 4])
write.table(bed.to.write, 
            file=file.path(projdir, "young_insulation_boundaries_400kb_filt.bed"),
            sep="\t", col.names = F, row.names = F, quote = F)

# Make scatter plot -------------------------------------------------------

tmp.all %>%
  filter(!(young==0 & aged==0)) %>%
  ggplot(aes(x=aged, y=young)) +
  geom_hex(bins=40) +
  geom_abline(linetype="dashed") +
  scale_fill_viridis_c() +
  coord_fixed(xlim=c(0,8.5), ylim=c(0,8.5)) +
  theme_bw() + 
  xlab("TAD boundary strength in Aged") +
  ylab("TAD boundary strength in Young") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=12, face="bold", color="black"),
        axis.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face="bold"))
ggsave(file.path(projdir, "insulation_scatter.png"), dpi=300, width=5, height=5)
